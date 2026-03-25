import geopandas as gpd
import pytest
from shapely.geometry import LineString, Point, Polygon

from osm_chordify.main import intersect_road_network_with_county_zones
from osm_chordify.osm.intersect import intersect_road_network_with_zones


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


def test_intersection_can_prefilter_zones_with_buffer(monkeypatch):
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
    original = __import__("osm_chordify.osm.intersect", fromlist=["_prefilter_zones_against_buffer"])._prefilter_zones_against_buffer

    def _wrapped_prefilter(polys_proj, edges_proj, road_buffer_filter_m):
        calls["count"] = calls.get("count", 0) + 1
        calls["road_buffer_filter_m"] = road_buffer_filter_m
        return original(polys_proj, edges_proj, road_buffer_filter_m)

    monkeypatch.setattr(
        "osm_chordify.osm.intersect._prefilter_zones_against_buffer",
        _wrapped_prefilter,
    )

    result = intersect_road_network_with_zones(
        road_network=edges,
        road_network_epsg=3857,
        zones=zones,
        road_buffer_filter_m=25.0,
    )

    assert len(result) == 1
    assert calls["count"] == 1
    assert calls["road_buffer_filter_m"] == 25.0


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
        road_buffer_filter_m=25.0,
    )

    assert len(result) == 1
    assert result.iloc[0]["zone_zone_id"] == "near"


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
    assert len(result) == 1
    assert result.iloc[0]["zone_GEOID"] == "06001"


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
