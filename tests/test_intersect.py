import geopandas as gpd
import pytest
from shapely.geometry import LineString, Point, Polygon

from osm_chordify.osm.intersect import intersect_road_network_with_zones


def test_intersection_includes_proportion_value_columns():
    edges = gpd.GeoDataFrame(
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
        proportional_cols=["vmt"],
    )

    assert len(result) == 1
    assert "proportion" in result.columns
    assert "edge_length_m" not in result.columns
    assert "proportional_length_m" not in result.columns
    assert "proportional_vmt" in result.columns

    assert result.iloc[0]["proportion"] == pytest.approx(0.5, abs=1e-6)
    assert result.iloc[0]["proportional_vmt"] == pytest.approx(50.0, abs=1e-6)


def test_intersection_includes_length_columns_when_requested():
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
        proportional_cols=["edge_length_m", "proportional_length_m"],
    )

    assert len(result) == 1
    assert "proportion" in result.columns
    assert "edge_length_m" in result.columns
    assert "proportional_length_m" in result.columns

    assert result.iloc[0]["edge_length_m"] == pytest.approx(10.0, abs=1e-6)
    assert result.iloc[0]["proportional_length_m"] == pytest.approx(5.0, abs=1e-6)


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

    def test_empty_result_has_proportion_column(self):
        result = intersect_road_network_with_zones(
            road_network=_make_edges(),
            road_network_epsg=3857,
            zones=_make_non_overlapping_zone(),
        )
        assert "proportion" in result.columns

    def test_empty_result_column_order_matches_nonempty_with_length_cols(self):
        """
        Regression: empty result with both edge_length_m and proportional_length_m
        must have the same column order as a non-empty result.

        Old code: insert(0, "edge_length_m") put it BEFORE "proportion".
        Fixed code: append() puts it AFTER "proportion" and "geometry",
        matching the record dict order in the non-empty path.
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
            proportional_cols=["edge_length_m", "proportional_length_m"],
        )
        empty = intersect_road_network_with_zones(
            road_network=_make_edges(),
            road_network_epsg=3857,
            zones=_make_non_overlapping_zone(),
            proportional_cols=["edge_length_m", "proportional_length_m"],
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

    def test_empty_result_edge_length_m_before_proportional_length_m(self):
        """edge_length_m must come before proportional_length_m in the empty result."""
        result = intersect_road_network_with_zones(
            road_network=_make_edges(),
            road_network_epsg=3857,
            zones=_make_non_overlapping_zone(),
            proportional_cols=["edge_length_m", "proportional_length_m"],
        )
        cols = list(result.columns)
        assert cols.index("edge_length_m") < cols.index("proportional_length_m")
