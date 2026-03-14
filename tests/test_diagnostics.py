"""
Regression tests for diagnostics.py.

Covers the data-loss bug in check_duplicate_edge_ids where the function
was overwriting the duplicate_info DataFrame with a tuple before returning it.
"""

import geopandas as gpd
import pandas as pd
import pytest
from shapely.geometry import LineString

from osm_chordify.osm.diagnostics import check_duplicate_edge_ids


def _make_edges(edge_ids, osmids=None):
    """Build a minimal edges GeoDataFrame with the given edge_id values."""
    if osmids is None:
        osmids = list(range(len(edge_ids)))
    return gpd.GeoDataFrame(
        {
            "edge_id": edge_ids,
            "osmid": osmids,
            "geometry": [LineString([(i, 0), (i + 1, 0)]) for i in range(len(edge_ids))],
        },
        geometry="geometry",
        crs="EPSG:4326",
    )


class TestCheckDuplicateEdgeIds:
    """
    Regression: check_duplicate_edge_ids used to overwrite the `duplicate_info`
    DataFrame with a tuple (duplicate_info, examples_df) before returning it.
    Any caller working with the result as a DataFrame would crash.
    """

    def test_no_duplicates_returns_false_and_none(self):
        edges = _make_edges(["a", "b", "c"])
        has_dupes, info = check_duplicate_edge_ids(edges)
        assert has_dupes is False
        assert info is None

    def test_duplicates_found_returns_true(self):
        edges = _make_edges(["a", "a", "b"])
        has_dupes, info = check_duplicate_edge_ids(edges)
        assert has_dupes is True

    # --- regression: return type must be DataFrame, not tuple ---------------

    def test_duplicate_info_is_dataframe_not_tuple(self):
        """
        Regression: the function overwrote duplicate_info with a tuple before
        returning.  The return type must be a DataFrame so callers can use
        .columns, .iloc, etc.
        """
        edges = _make_edges(["dup", "dup", "unique"])
        has_dupes, info = check_duplicate_edge_ids(edges)
        assert has_dupes is True
        assert isinstance(info, pd.DataFrame), (
            f"Expected DataFrame but got {type(info).__name__} — "
            "check_duplicate_edge_ids tuple overwrite bug regressed"
        )

    def test_duplicate_info_has_expected_columns(self):
        edges = _make_edges(["x", "x", "y"])
        _, info = check_duplicate_edge_ids(edges)
        assert "edge_id" in info.columns
        assert "count" in info.columns

    def test_duplicate_info_count_is_correct(self):
        edges = _make_edges(["dup", "dup", "dup", "ok"])
        _, info = check_duplicate_edge_ids(edges)
        row = info[info["edge_id"] == "dup"].iloc[0]
        assert row["count"] == 3

    def test_multiple_duplicate_groups(self):
        edges = _make_edges(["a", "a", "b", "b", "b", "c"])
        has_dupes, info = check_duplicate_edge_ids(edges)
        assert has_dupes is True
        assert len(info) == 2  # two duplicate groups: a and b

    def test_custom_id_column(self):
        edges = _make_edges(["x", "x", "y"])
        edges = edges.rename(columns={"edge_id": "my_id"})
        has_dupes, info = check_duplicate_edge_ids(edges, id_column="my_id")
        assert has_dupes is True
        assert isinstance(info, pd.DataFrame)