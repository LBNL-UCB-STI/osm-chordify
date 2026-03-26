import geopandas as gpd
import pandas as pd
import pytest
from shapely.geometry import LineString, Point, Polygon
import json

from osm_chordify.main import build_area_mask_from_counties, intersect_road_network_with_county_zones
from osm_chordify.osm.intersect import (
    intersect_polygons_with_zones,
    intersect_road_network_with_zones,
    intersect_road_polygons_with_zones,
    spatial_left_join_with_zones,
    intersect_zones_with_zones,
)


def test_intersection_includes_fixed_length_and_zone_edge_proportion_columns():
    edges = gpd.GeoDataFrame(
        {
            "osm_id": [1],
            "edge_id": [101],
            "edge_length": [10.0],
            "geometry": [LineString([(0, 0), (10, 0)])],
        },
        geometry="geometry",
        crs="EPSG:3857",
    )

    zones = gpd.GeoDataFrame(
        {
            "zone_id": ["A"],
            "geometry": [Polygon([(0, -1), (5, -1), (5, 1), (0, 1)])],
        },
        geometry="geometry",
        crs="EPSG:3857",
    )

    result = intersect_road_network_with_zones(
        road_network=edges,
        road_network_epsg=3857,
        zones=zones,
    )

    assert len(result) == 1
    assert "zone_edge_proportion" in result.columns
    assert "edge_link_length_m" in result.columns
    assert "zone_link_length_m" in result.columns
    assert result.iloc[0]["edge_link_length_m"] == pytest.approx(10.0, abs=1e-6)
    assert result.iloc[0]["zone_link_length_m"] == pytest.approx(5.0, abs=1e-6)

    assert result.iloc[0]["zone_edge_proportion"] == pytest.approx(0.5, abs=1e-6)


def test_intersection_always_includes_length_columns():
    edges = gpd.GeoDataFrame(
        {
            "osm_id": [1],
            "edge_id": [101],
            "edge_length": [10.0],
            "geometry": [LineString([(0, 0), (10, 0)])],
        },
        geometry="geometry",
        crs="EPSG:3857",
    )

    zones = gpd.GeoDataFrame(
        {
            "zone_id": ["A"],
            "geometry": [Polygon([(0, -1), (5, -1), (5, 1), (0, 1)])],
        },
        geometry="geometry",
        crs="EPSG:3857",
    )

    result = intersect_road_network_with_zones(
        road_network=edges,
        road_network_epsg=3857,
        zones=zones,
    )

    assert len(result) == 1
    assert "zone_edge_proportion" in result.columns
    assert "edge_link_length_m" in result.columns
    assert "zone_link_length_m" in result.columns

    assert result.iloc[0]["edge_link_length_m"] == pytest.approx(10.0, abs=1e-6)
    assert result.iloc[0]["zone_link_length_m"] == pytest.approx(5.0, abs=1e-6)


def test_intersection_rejects_geographic_working_crs():
    edges = gpd.GeoDataFrame(
        {
            "osm_id": [1],
            "edge_id": [101],
            "edge_length": [10.0],
            "geometry": [LineString([(-122.0, 37.0), (-121.99, 37.0)])],
        },
        geometry="geometry",
        crs="EPSG:4326",
    )
    zones = gpd.GeoDataFrame(
        {
            "zone_id": ["A"],
            "geometry": [Polygon([(-122.0, 36.99), (-121.995, 36.99), (-121.995, 37.01), (-122.0, 37.01)])],
        },
        geometry="geometry",
        crs="EPSG:4326",
    )

    with pytest.raises(ValueError, match="geographic"):
        intersect_road_network_with_zones(
            road_network=edges,
            road_network_epsg=4326,
            zones=zones,
        )


def test_intersection_can_prefilter_zones_to_network_bbox(monkeypatch):
    edges = gpd.GeoDataFrame(
        {
            "osm_id": [1],
            "edge_id": [101],
            "edge_length": [10.0],
            "geometry": [LineString([(0, 0), (10, 0)])],
        },
        geometry="geometry",
        crs="EPSG:3857",
    )
    zones = gpd.GeoDataFrame(
        {
            "zone_id": ["A"],
            "geometry": [Polygon([(0, -1), (5, -1), (5, 1), (0, 1)])],
        },
        geometry="geometry",
        crs="EPSG:3857",
    )

    calls = {}
    original = __import__("osm_chordify.osm.intersect", fromlist=["_prefilter_zones_to_network_bbox"])._prefilter_zones_to_network_bbox

    def _wrapped_prefilter(polys_proj, edges_proj):
        calls["count"] = calls.get("count", 0) + 1
        return original(polys_proj, edges_proj)

    monkeypatch.setattr(
        "osm_chordify.osm.intersect._prefilter_zones_to_network_bbox",
        _wrapped_prefilter,
    )

    result = intersect_road_network_with_zones(
        road_network=edges,
        road_network_epsg=3857,
        zones=zones,
        prefilter_zones_to_network_bbox=True,
    )

    assert len(result) == 1
    assert calls["count"] == 1


def test_prefilter_drops_far_zones_and_keeps_nearby_zones():
    edges = gpd.GeoDataFrame(
        {
            "osm_id": [1],
            "edge_id": [101],
            "edge_length": [10.0],
            "geometry": [LineString([(0, 0), (10, 0)])],
        },
        geometry="geometry",
        crs="EPSG:3857",
    )
    zones = gpd.GeoDataFrame(
        {
            "zone_id": ["near", "far"],
            "geometry": [
                Polygon([(0, -2), (5, -2), (5, 2), (0, 2)]),
                Polygon([(1000, 1000), (1010, 1000), (1010, 1010), (1000, 1010)]),
            ],
        },
        geometry="geometry",
        crs="EPSG:3857",
    )

    result = intersect_road_network_with_zones(
        road_network=edges,
        road_network_epsg=3857,
        zones=zones,
        prefilter_zones_to_network_bbox=True,
    )

    assert len(result) == 1
    assert result.iloc[0]["zone_zone_id"] == "near"


def test_bbox_prefilter_keeps_voided_rows_for_bbox_cells_without_link_pieces():
    edges = gpd.GeoDataFrame(
        {
            "osm_id": [1],
            "edge_id": [101],
            "edge_length": [10.0],
            "geometry": [LineString([(0, 0), (10, 0)])],
        },
        geometry="geometry",
        crs="EPSG:3857",
    )
    zones = gpd.GeoDataFrame(
        {
            "zone_id": ["hit", "voided", "outside"],
            "geometry": [
                Polygon([(0, -1), (5, -1), (5, 1), (0, 1)]),
                Polygon([(10, -1), (12, -1), (12, 1), (10, 1)]),
                Polygon([(1000, 1000), (1010, 1000), (1010, 1010), (1000, 1010)]),
            ],
        },
        geometry="geometry",
        crs="EPSG:3857",
    )

    result = intersect_road_network_with_zones(
        road_network=edges,
        road_network_epsg=3857,
        zones=zones,
        prefilter_zones_to_network_bbox=True,
    )

    assert set(result["zone_zone_id"]) == {"hit", "voided"}
    voided = result[result["zone_zone_id"] == "voided"].iloc[0]
    assert pd.isna(voided["zone_edge_proportion"])
    assert pd.isna(voided["edge_link_length_m"])
    assert pd.isna(voided["zone_link_length_m"])


def _make_edges():
    return gpd.GeoDataFrame(
        {
            "osm_id": [1],
            "edge_id": [101],
            "edge_length": [10.0],
            "vmt": [100.0],
            "geometry": [LineString([(0, 0), (10, 0)])],
        },
        geometry="geometry",
        crs="EPSG:3857",
    )


def _make_non_overlapping_zone():
    """A zone far from the edge so no intersections occur."""
    return gpd.GeoDataFrame(
        {
            "zone_id": ["Z"],
            "geometry": [Polygon([(1000, 1000), (2000, 1000), (2000, 2000), (1000, 2000)])],
        },
        geometry="geometry",
        crs="EPSG:3857",
    )


class TestEmptyIntersectionResult:
    """
    Regression (M2): the empty-result path used list.insert() with wrong index
    values, producing a different column order than the non-empty path.
    """

    def test_empty_result_returns_geodataframe(self):
        result = intersect_road_network_with_zones(
            road_network=_make_edges(),
            road_network_epsg=3857,
            zones=_make_non_overlapping_zone(),
        )
        assert isinstance(result, gpd.GeoDataFrame)
        assert len(result) == 0

    def test_empty_result_has_zone_edge_proportion_column(self):
        result = intersect_road_network_with_zones(
            road_network=_make_edges(),
            road_network_epsg=3857,
            zones=_make_non_overlapping_zone(),
        )
        assert "zone_edge_proportion" in result.columns

    def test_empty_result_column_order_matches_nonempty(self):
        """
        Regression: empty result with the always-present link-length columns
        must have the same column order as a non-empty result.
        """
        zone_overlapping = gpd.GeoDataFrame(
            {
                "zone_id": ["A"],
                "geometry": [Polygon([(0, -1), (5, -1), (5, 1), (0, 1)])],
            },
            geometry="geometry",
            crs="EPSG:3857",
        )
        nonempty = intersect_road_network_with_zones(
            road_network=_make_edges(),
            road_network_epsg=3857,
            zones=zone_overlapping,
        )
        empty = intersect_road_network_with_zones(
            road_network=_make_edges(),
            road_network_epsg=3857,
            zones=_make_non_overlapping_zone(),
        )
        # The shared columns must appear in the same order in both results
        shared = [c for c in nonempty.columns if c in empty.columns]
        nonempty_order = [c for c in nonempty.columns if c in shared]
        empty_order = [c for c in empty.columns if c in shared]
        assert nonempty_order == empty_order, (
            f"Column order mismatch between empty and non-empty results.\n"
            f"  non-empty: {nonempty_order}\n"
            f"  empty:     {empty_order}"
        )

    def test_empty_result_edge_link_length_m_before_zone_link_length_m(self):
        """edge_link_length_m must come before zone_link_length_m in the empty result."""
        result = intersect_road_network_with_zones(
            road_network=_make_edges(),
            road_network_epsg=3857,
            zones=_make_non_overlapping_zone(),
        )
        cols = list(result.columns)
        assert cols.index("edge_link_length_m") < cols.index("zone_link_length_m")


def test_intersect_road_network_with_county_zones_uses_fips_boundary_fetch(monkeypatch):
    edges = gpd.GeoDataFrame(
        {
            "osm_id": [1],
            "edge_id": [101],
            "edge_length": [10.0],
            "geometry": [LineString([(0, 0), (10, 0)])],
        },
        geometry="geometry",
        crs="EPSG:3857",
    )
    counties = gpd.GeoDataFrame(
        {
            "GEOID": ["06001"],
            "NAME": ["Alameda"],
            "geometry": [Polygon([(0, -1), (5, -1), (5, 1), (0, 1)])],
        },
        geometry="geometry",
        crs="EPSG:3857",
    )

    calls = {}

    def _fake_collect_geographic_boundaries(**kwargs):
        calls["kwargs"] = kwargs
        return counties

    monkeypatch.setattr(
        "osm_chordify.utils.data_collection.collect_geographic_boundaries",
        _fake_collect_geographic_boundaries,
    )

    result = intersect_road_network_with_county_zones(
        road_network=edges,
        road_network_epsg=3857,
        state_fips_code="06",
        county_fips_codes=["001"],
        year=2020,
        work_dir="/tmp/osm-chordify-test",
        area_name="county-test",
    )

    assert calls["kwargs"]["state_fips_code"] == "06"
    assert calls["kwargs"]["county_fips_codes"] == ["001"]
    assert calls["kwargs"]["geo_level"] == "county"
    assert calls["kwargs"]["target_epsg"] == 3857
    assert calls["kwargs"]["cartographic"] is True
    assert len(result) == 1
    assert result.iloc[0]["zone_GEOID"] == "06001"


def test_build_area_mask_from_counties_land_only_preserves_gap(monkeypatch):
    counties = gpd.GeoDataFrame(
        {
            "GEOID": ["06001", "06013"],
            "geometry": [
                Polygon([(0, 0), (1, 0), (1, 1), (0, 1)]),
                Polygon([(3, 0), (4, 0), (4, 1), (3, 1)]),
            ],
        },
        geometry="geometry",
        crs="EPSG:3857",
    )

    calls = {}

    def _fake_collect_geographic_boundaries(**kwargs):
        calls["kwargs"] = kwargs
        return counties

    monkeypatch.setattr(
        "osm_chordify.utils.data_collection.collect_geographic_boundaries",
        _fake_collect_geographic_boundaries,
    )

    mask = build_area_mask_from_counties(
        state_fips_code="06",
        county_fips_codes=["001", "013"],
        year=2020,
        work_dir="/tmp/osm-chordify-test",
        output_epsg=3857,
        include_water=False,
    )

    assert len(mask) == 1
    assert calls["kwargs"]["cartographic"] is True
    assert not mask.iloc[0].geometry.contains(Point(2, 0.5))


def test_build_area_mask_from_counties_whole_area_fills_gap(monkeypatch):
    counties = gpd.GeoDataFrame(
        {
            "GEOID": ["06001", "06013"],
            "geometry": [
                Polygon([(0, 0), (1, 0), (1, 1), (0, 1)]),
                Polygon([(3, 0), (4, 0), (4, 1), (3, 1)]),
            ],
        },
        geometry="geometry",
        crs="EPSG:3857",
    )

    calls = {}

    def _fake_collect_geographic_boundaries(**kwargs):
        calls["kwargs"] = kwargs
        return counties

    monkeypatch.setattr(
        "osm_chordify.utils.data_collection.collect_geographic_boundaries",
        _fake_collect_geographic_boundaries,
    )

    mask = build_area_mask_from_counties(
        state_fips_code="06",
        county_fips_codes=["001", "013"],
        year=2020,
        work_dir="/tmp/osm-chordify-test",
        output_epsg=3857,
        include_water=True,
    )

    assert len(mask) == 1
    assert calls["kwargs"]["cartographic"] is False
    assert mask.iloc[0].geometry.contains(Point(2, 0.5))


def test_intersection_drops_point_only_boundary_touches():
    edges = gpd.GeoDataFrame(
        {
            "osm_id": [1],
            "edge_id": [101],
            "edge_length": [10.0],
            "geometry": [LineString([(0, 0), (10, 0)])],
        },
        geometry="geometry",
        crs="EPSG:3857",
    )
    zones = gpd.GeoDataFrame(
        {
            "zone_id": ["A"],
            "geometry": [Polygon([(10, 0), (11, 0), (11, 1), (10, 1)])],
        },
        geometry="geometry",
        crs="EPSG:3857",
    )

    result = intersect_road_network_with_zones(
        road_network=edges,
        road_network_epsg=3857,
        zones=zones,
    )

    assert len(result) == 0


def test_intersect_road_polygons_with_zones_uses_area_proportion_to_derive_length():
    road_polygons = gpd.GeoDataFrame(
        {
            "edge_id": [101],
            "edge_length": [10.0],
            "geometry": [Polygon([(0, 0), (10, 0), (10, 2), (0, 2)])],
        },
        geometry="geometry",
        crs="EPSG:3857",
    )
    zones = gpd.GeoDataFrame(
        {
            "zone_id": ["left-half"],
            "geometry": [Polygon([(0, 0), (5, 0), (5, 2), (0, 2)])],
        },
        geometry="geometry",
        crs="EPSG:3857",
    )

    result = intersect_road_polygons_with_zones(
        road_network=road_polygons,
        road_network_epsg=3857,
        zones=zones,
        zone_label="aermod",
    )

    assert len(result) == 1
    assert result.iloc[0]["aermod_zone_edge_proportion"] == pytest.approx(0.5, abs=1e-6)
    assert result.iloc[0]["aermod_edge_link_length_m"] == pytest.approx(10.0, abs=1e-6)
    assert result.iloc[0]["aermod_zone_link_length_m"] == pytest.approx(5.0, abs=1e-6)
    assert result.iloc[0]["aermod_edge_surface_m2"] == pytest.approx(20.0, abs=1e-6)
    assert result.iloc[0]["aermod_zone_surface_m2"] == pytest.approx(10.0, abs=1e-6)
    assert result.iloc[0].geometry.geom_type in {"Polygon", "MultiPolygon"}
    assert result.iloc[0]["aermod_zone_id"] == "left-half"


def test_intersect_road_polygons_with_zones_drops_boundary_only_touches():
    road_polygons = gpd.GeoDataFrame(
        {
            "edge_id": [101],
            "edge_length": [10.0],
            "geometry": [Polygon([(0, 0), (10, 0), (10, 2), (0, 2)])],
        },
        geometry="geometry",
        crs="EPSG:3857",
    )
    zones = gpd.GeoDataFrame(
        {
            "zone_id": ["touch-only"],
            "geometry": [Polygon([(10, 0), (12, 0), (12, 2), (10, 2)])],
        },
        geometry="geometry",
        crs="EPSG:3857",
    )

    result = intersect_road_polygons_with_zones(
        road_network=road_polygons,
        road_network_epsg=3857,
        zones=zones,
    )

    assert len(result) == 0


def test_intersect_road_polygons_with_zones_prefilter_keeps_void_rows():
    road_polygons = gpd.GeoDataFrame(
        {
            "edge_id": [101],
            "edge_length": [10.0],
            "geometry": [Polygon([(0, 0), (10, 0), (10, 2), (0, 2)])],
        },
        geometry="geometry",
        crs="EPSG:3857",
    )
    zones = gpd.GeoDataFrame(
        {
            "zone_id": ["hit", "voided", "outside"],
            "geometry": [
                Polygon([(0, 0), (5, 0), (5, 2), (0, 2)]),
                Polygon([(0, 2), (5, 2), (5, 4), (0, 4)]),
                Polygon([(100, 100), (105, 100), (105, 105), (100, 105)]),
            ],
        },
        geometry="geometry",
        crs="EPSG:3857",
    )

    result = intersect_road_polygons_with_zones(
        road_network=road_polygons,
        road_network_epsg=3857,
        zones=zones,
        prefilter_zones_to_network_bbox=True,
    )

    assert set(result["zone_zone_id"]) == {"hit", "voided"}
    voided = result[result["zone_zone_id"] == "voided"].iloc[0]
    assert pd.isna(voided["zone_edge_proportion"])
    assert pd.isna(voided["edge_link_length_m"])
    assert pd.isna(voided["zone_link_length_m"])


def test_intersect_road_polygons_with_zones_prefilters_to_bbox_by_default(monkeypatch):
    road_polygons = gpd.GeoDataFrame(
        {
            "edge_id": [101],
            "edge_length": [10.0],
            "geometry": [Polygon([(0, 0), (10, 0), (10, 2), (0, 2)])],
        },
        geometry="geometry",
        crs="EPSG:3857",
    )
    zones = gpd.GeoDataFrame(
        {
            "zone_id": ["hit"],
            "geometry": [Polygon([(0, 0), (5, 0), (5, 2), (0, 2)])],
        },
        geometry="geometry",
        crs="EPSG:3857",
    )

    calls = {}
    original = __import__("osm_chordify.osm.intersect", fromlist=["_prefilter_zones_to_network_bbox"])._prefilter_zones_to_network_bbox

    def _wrapped_prefilter(polys_proj, edges_proj):
        calls["count"] = calls.get("count", 0) + 1
        return original(polys_proj, edges_proj)

    monkeypatch.setattr(
        "osm_chordify.osm.intersect._prefilter_zones_to_network_bbox",
        _wrapped_prefilter,
    )

    result = intersect_road_polygons_with_zones(
        road_network=road_polygons,
        road_network_epsg=3857,
        zones=zones,
    )

    assert len(result) == 1
    assert calls["count"] == 1


def test_intersect_polygons_with_zones_preserves_existing_columns_and_recomputes_piece_metrics():
    polygons = gpd.GeoDataFrame(
        {
            "edge_osm_id": [1],
            "zone_NAME": ["inmap-cell-1"],
            "zone_link_length_m": [10.0],
            "geometry": [Polygon([(0, 0), (10, 0), (10, 2), (0, 2)])],
        },
        geometry="geometry",
        crs="EPSG:3857",
    )
    zones = gpd.GeoDataFrame(
        {
            "aermod_id": ["a1"],
            "geometry": [Polygon([(0, 0), (5, 0), (5, 2), (0, 2)])],
        },
        geometry="geometry",
        crs="EPSG:3857",
    )

    result = intersect_polygons_with_zones(
        polygons=polygons,
        polygons_epsg=3857,
        zones=zones,
        zone_label="inmap",
    )

    assert len(result) == 1
    assert result.iloc[0]["edge_osm_id"] == 1
    assert result.iloc[0]["zone_NAME"] == "inmap-cell-1"
    assert result.iloc[0]["inmap_aermod_id"] == "a1"
    assert result.iloc[0]["inmap_piece_link_length_m"] == pytest.approx(10.0, abs=1e-6)
    assert result.iloc[0]["inmap_zone_piece_proportion"] == pytest.approx(0.5, abs=1e-6)
    assert result.iloc[0]["inmap_zone_piece_length_m"] == pytest.approx(5.0, abs=1e-6)
    assert result.iloc[0]["inmap_piece_surface_m2"] == pytest.approx(20.0, abs=1e-6)
    assert result.iloc[0]["inmap_zone_surface_m2"] == pytest.approx(10.0, abs=1e-6)


def test_spatial_left_join_with_zones_keeps_unmatched_rows_with_null_zone_fields():
    pieces = gpd.GeoDataFrame(
        {
            "piece_id": [1, 2],
            "geometry": [
                Polygon([(0, 0), (2, 0), (2, 2), (0, 2)]),
                Polygon([(10, 10), (12, 10), (12, 12), (10, 12)]),
            ],
        },
        geometry="geometry",
        crs="EPSG:3857",
    )
    counties = gpd.GeoDataFrame(
        {
            "COUNTYFP": ["001"],
            "geometry": [Polygon([(0, 0), (3, 0), (3, 3), (0, 3)])],
        },
        geometry="geometry",
        crs="EPSG:3857",
    )

    result = spatial_left_join_with_zones(
        gdf=pieces,
        gdf_epsg=3857,
        zones=counties,
        zone_label="county",
    )

    assert len(result) == 2
    matched = result[result["piece_id"] == 1].iloc[0]
    unmatched = result[result["piece_id"] == 2].iloc[0]
    assert matched["county_COUNTYFP"] == "001"
    assert pd.isna(unmatched["county_COUNTYFP"])


def test_intersection_preserves_prior_prefixed_columns_without_edge_stacking():
    edges = gpd.GeoDataFrame(
        {
            "edge_osm_id": [1],
            "zone_NAME": ["old_zone"],
            "zone_edge_proportion": [0.6],
            "edge_link_length_m": [10.0],
            "zone_link_length_m": [6.0],
            "geometry": [LineString([(0, 0), (10, 0)])],
        },
        geometry="geometry",
        crs="EPSG:3857",
    )
    zones = gpd.GeoDataFrame(
        {
            "NAME": ["new_zone"],
            "geometry": [Polygon([(0, -1), (5, -1), (5, 1), (0, 1)])],
        },
        geometry="geometry",
        crs="EPSG:3857",
    )

    result = intersect_road_network_with_zones(
        road_network=edges,
        road_network_epsg=3857,
        zones=zones,
    )

    assert "edge_edge_osm_id" not in result.columns
    assert "zone_zone_NAME" not in result.columns
    assert "edge_osm_id" in result.columns
    assert "zone_NAME" in result.columns
    assert "zone2_NAME" in result.columns
    assert result.iloc[0]["zone_NAME"] == "old_zone"
    assert result.iloc[0]["zone2_NAME"] == "new_zone"


def test_intersection_loads_parquet_inputs(tmp_path):
    edges = gpd.GeoDataFrame(
        {
            "osm_id": [1],
            "edge_id": [101],
            "edge_length": [10.0],
            "geometry": [LineString([(0, 0), (10, 0)])],
        },
        geometry="geometry",
        crs="EPSG:3857",
    )
    zones = gpd.GeoDataFrame(
        {
            "zone_id": ["A"],
            "geometry": [Polygon([(0, -1), (5, -1), (5, 1), (0, 1)])],
        },
        geometry="geometry",
        crs="EPSG:3857",
    )

    try:
        edges_path = tmp_path / "edges.parquet"
        zones_path = tmp_path / "zones.parquet"
        edges.to_parquet(edges_path, index=False)
        zones.to_parquet(zones_path, index=False)
    except ImportError:
        pytest.skip("pyarrow is required for parquet intersection input test")

    result = intersect_road_network_with_zones(
        road_network=edges_path,
        road_network_epsg=3857,
        zones=zones_path,
    )

    assert len(result) == 1
    assert result.iloc[0]["zone_edge_proportion"] == pytest.approx(0.5, abs=1e-6)


def test_intersection_reuses_cached_output(monkeypatch, tmp_path):
    edges = gpd.GeoDataFrame(
        {
            "osm_id": [1],
            "edge_id": [101],
            "edge_length": [10.0],
            "geometry": [LineString([(0, 0), (10, 0)])],
        },
        geometry="geometry",
        crs="EPSG:3857",
    )
    zones = gpd.GeoDataFrame(
        {
            "zone_id": ["A"],
            "geometry": [Polygon([(0, -1), (5, -1), (5, 1), (0, 1)])],
        },
        geometry="geometry",
        crs="EPSG:3857",
    )
    output_path = tmp_path / "cached_intersection.geojson"

    first = intersect_road_network_with_zones(
        road_network=edges,
        road_network_epsg=3857,
        zones=zones,
        output_path=output_path,
    )
    assert len(first) == 1
    assert output_path.exists()
    assert (tmp_path / "cached_intersection.geojson.cache.json").exists()

    def _fail_sjoin(*args, **kwargs):
        raise AssertionError("spatial join should not run when cache is reused")

    monkeypatch.setattr(gpd, "sjoin", _fail_sjoin)

    second = intersect_road_network_with_zones(
        road_network=edges,
        road_network_epsg=3857,
        zones=zones,
        output_path=output_path,
    )
    assert len(second) == 1
    assert second.iloc[0]["zone_edge_proportion"] == pytest.approx(0.5, abs=1e-6)


def test_intersect_zones_with_zones_returns_prefixed_polygon_overlaps():
    zones_a = gpd.GeoDataFrame(
        {
            "county_id": ["001"],
            "geometry": [Polygon([(0, 0), (4, 0), (4, 4), (0, 4)])],
        },
        geometry="geometry",
        crs="EPSG:3857",
    )
    zones_b = gpd.GeoDataFrame(
        {
            "grid_id": ["cell-1"],
            "geometry": [Polygon([(2, 2), (6, 2), (6, 6), (2, 6)])],
        },
        geometry="geometry",
        crs="EPSG:3857",
    )

    result = intersect_zones_with_zones(
        zones_a=zones_a,
        zones_a_epsg=3857,
        zones_b=zones_b,
        zones_b_epsg=3857,
    )

    assert len(result) == 1
    assert "zone_a_county_id" in result.columns
    assert "zone_b_grid_id" in result.columns
    assert result.iloc[0]["zone_a_county_id"] == "001"
    assert result.iloc[0]["zone_b_grid_id"] == "cell-1"
    assert result.iloc[0].geometry.geom_type in {"Polygon", "MultiPolygon"}
    assert result.iloc[0].geometry.area == pytest.approx(4.0, abs=1e-6)


def test_intersect_zones_with_zones_drops_line_or_point_only_touches():
    zones_a = gpd.GeoDataFrame(
        {
            "county_id": ["001"],
            "geometry": [Polygon([(0, 0), (2, 0), (2, 2), (0, 2)])],
        },
        geometry="geometry",
        crs="EPSG:3857",
    )
    zones_b = gpd.GeoDataFrame(
        {
            "grid_id": ["cell-1"],
            "geometry": [Polygon([(2, 0), (4, 0), (4, 2), (2, 2)])],
        },
        geometry="geometry",
        crs="EPSG:3857",
    )

    result = intersect_zones_with_zones(
        zones_a=zones_a,
        zones_a_epsg=3857,
        zones_b=zones_b,
        zones_b_epsg=3857,
    )

    assert len(result) == 0


def test_intersect_zones_with_zones_saves_output(tmp_path):
    zones_a = gpd.GeoDataFrame(
        {
            "county_id": ["001"],
            "geometry": [Polygon([(0, 0), (4, 0), (4, 4), (0, 4)])],
        },
        geometry="geometry",
        crs="EPSG:3857",
    )
    zones_b = gpd.GeoDataFrame(
        {
            "grid_id": ["cell-1"],
            "geometry": [Polygon([(2, 2), (6, 2), (6, 6), (2, 6)])],
        },
        geometry="geometry",
        crs="EPSG:3857",
    )
    output_path = tmp_path / "zone_overlap.geojson"

    result = intersect_zones_with_zones(
        zones_a=zones_a,
        zones_a_epsg=3857,
        zones_b=zones_b,
        zones_b_epsg=3857,
        output_path=output_path,
    )

    loaded = gpd.read_file(output_path)
    assert output_path.exists()
    assert len(result) == 1
    assert len(loaded) == 1
    assert loaded.iloc[0]["zone_a_county_id"] == "001"
    assert loaded.iloc[0]["zone_b_grid_id"] == "cell-1"


def test_intersect_zones_with_zones_loads_parquet_inputs(tmp_path):
    zones_a = gpd.GeoDataFrame(
        {
            "county_id": ["001"],
            "geometry": [Polygon([(0, 0), (4, 0), (4, 4), (0, 4)])],
        },
        geometry="geometry",
        crs="EPSG:3857",
    )
    zones_b = gpd.GeoDataFrame(
        {
            "grid_id": ["cell-1"],
            "geometry": [Polygon([(2, 2), (6, 2), (6, 6), (2, 6)])],
        },
        geometry="geometry",
        crs="EPSG:3857",
    )

    try:
        zones_a_path = tmp_path / "zones_a.parquet"
        zones_b_path = tmp_path / "zones_b.parquet"
        zones_a.to_parquet(zones_a_path, index=False)
        zones_b.to_parquet(zones_b_path, index=False)
    except ImportError:
        pytest.skip("pyarrow is required for parquet zone intersection input test")

    result = intersect_zones_with_zones(
        zones_a=zones_a_path,
        zones_a_epsg=3857,
        zones_b=zones_b_path,
        zones_b_epsg=3857,
    )

    assert len(result) == 1
    assert result.iloc[0]["zone_a_county_id"] == "001"
    assert result.iloc[0]["zone_b_grid_id"] == "cell-1"


def test_intersect_zones_with_zones_distinguishes_land_and_whole_area_masks(monkeypatch):
    counties = gpd.GeoDataFrame(
        {
            "GEOID": ["06001", "06013"],
            "geometry": [
                Polygon([(0, 0), (1, 0), (1, 1), (0, 1)]),
                Polygon([(3, 0), (4, 0), (4, 1), (3, 1)]),
            ],
        },
        geometry="geometry",
        crs="EPSG:3857",
    )

    def _fake_collect_geographic_boundaries(**kwargs):
        return counties

    monkeypatch.setattr(
        "osm_chordify.utils.data_collection.collect_geographic_boundaries",
        _fake_collect_geographic_boundaries,
    )

    land_mask = build_area_mask_from_counties(
        state_fips_code="06",
        county_fips_codes=["001", "013"],
        year=2020,
        work_dir="/tmp/osm-chordify-test",
        output_epsg=3857,
        include_water=False,
    )
    whole_area_mask = build_area_mask_from_counties(
        state_fips_code="06",
        county_fips_codes=["001", "013"],
        year=2020,
        work_dir="/tmp/osm-chordify-test",
        output_epsg=3857,
        include_water=True,
    )
    gap_cell = gpd.GeoDataFrame(
        {
            "grid_id": ["gap-cell"],
            "geometry": [Polygon([(1.5, 0.25), (2.5, 0.25), (2.5, 0.75), (1.5, 0.75)])],
        },
        geometry="geometry",
        crs="EPSG:3857",
    )

    land_overlap = intersect_zones_with_zones(
        zones_a=land_mask,
        zones_a_epsg=3857,
        zones_b=gap_cell,
        zones_b_epsg=3857,
    )
    whole_area_overlap = intersect_zones_with_zones(
        zones_a=whole_area_mask,
        zones_a_epsg=3857,
        zones_b=gap_cell,
        zones_b_epsg=3857,
    )

    assert len(land_overlap) == 0
    assert len(whole_area_overlap) == 1
    assert whole_area_overlap.iloc[0]["zone_b_grid_id"] == "gap-cell"
