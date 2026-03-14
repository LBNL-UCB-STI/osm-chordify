"""
Unit tests for graph.process_tags().

process_tags() is the most business-critical function in the pipeline: it
standardizes oneway/maxspeed/access/motor_vehicle and computes hgv/mdv boolean
flags from weight restrictions.  It had zero test coverage before these tests.

All tests use synthetic osmnx-compatible graphs (no network access).
"""

import geopandas as gpd
import networkx as nx
import osmnx as ox
import pytest
from shapely.geometry import LineString, Point

from osm_chordify.osm.graph import process_tags

# ---------------------------------------------------------------------------
# Config helpers
# ---------------------------------------------------------------------------

_DEFAULT_CONFIG = {
    "weight_limits": {
        "unit": "lbs",
        "mdv_max": 26000,
        "hdv_max": 80000,
    }
}

_KG_CONFIG = {
    "weight_limits": {
        "unit": "kg",
        "mdv_max": 11793,   # 26000 lbs in kg
        "hdv_max": 36287,   # 80000 lbs in kg
    }
}


def _build_graph(node_coords, edges_data):
    """
    Build a minimal osmnx-compatible WGS-84 MultiDiGraph.

    Parameters
    ----------
    node_coords : dict  {node_id: (lon, lat)}
    edges_data  : list of dicts with keys: u, v, and any edge attributes
    """
    node_rows = [
        {"osmid": nid, "x": lon, "y": lat, "geometry": Point(lon, lat)}
        for nid, (lon, lat) in node_coords.items()
    ]
    nodes_gdf = (
        gpd.GeoDataFrame(node_rows, geometry="geometry", crs="EPSG:4326")
        .set_index("osmid")
    )

    edge_rows = []
    for ed in edges_data:
        u, v = ed["u"], ed["v"]
        lon1, lat1 = node_coords[u]
        lon2, lat2 = node_coords[v]
        row = {
            "u": u, "v": v, "key": 0,
            "osmid": ed.get("osmid", 1),
            "length": ed.get("length", 100.0),
            "highway": ed.get("highway", "residential"),
            "geometry": LineString([(lon1, lat1), (lon2, lat2)]),
        }
        # Pass through any additional columns
        for k, val in ed.items():
            if k not in ("u", "v"):
                row[k] = val
        edge_rows.append(row)

    edges_gdf = (
        gpd.GeoDataFrame(edge_rows, geometry="geometry", crs="EPSG:4326")
        .set_index(["u", "v", "key"])
    )
    return ox.graph_from_gdfs(nodes_gdf, edges_gdf)


_NODES = {1: (-122.0, 37.0), 2: (-122.001, 37.0)}

# process_tags() calls .apply() on these columns unconditionally, so they must
# always be present.  osmnx drops all-None columns on graph_from_gdfs →
# graph_to_gdfs roundtrips, so use neutral non-None sentinel values.
_BASE_EDGE = {
    "u": 1, "v": 2, "osmid": 100,
    "oneway": "no",
    "motor_vehicle": "yes",
    "access": "yes",
    "maxspeed": "50",
}


def _single_edge(**attrs):
    ed = {**_BASE_EDGE, **attrs}
    return _build_graph(_NODES, [ed])


def _get_edge_attr(g, attr):
    """Return the attribute value for the single edge in the graph."""
    _, edges = ox.graph_to_gdfs(g)
    return edges.iloc[0][attr]


# ---------------------------------------------------------------------------
# oneway standardization
# ---------------------------------------------------------------------------

class TestProcessTagsOneway:
    def test_oneway_yes_preserved(self):
        g = _single_edge(oneway="yes")
        result = process_tags(g, _DEFAULT_CONFIG)
        assert _get_edge_attr(result, "oneway") == "yes"

    def test_oneway_true_normalized_to_yes(self):
        g = _single_edge(oneway=True)
        result = process_tags(g, _DEFAULT_CONFIG)
        assert _get_edge_attr(result, "oneway") == "yes"

    def test_oneway_false_normalized_to_no(self):
        g = _single_edge(oneway=False)
        result = process_tags(g, _DEFAULT_CONFIG)
        assert _get_edge_attr(result, "oneway") == "no"

    def test_oneway_minus1_preserved(self):
        g = _single_edge(oneway="-1")
        result = process_tags(g, _DEFAULT_CONFIG)
        assert _get_edge_attr(result, "oneway") == "-1"

    def test_oneway_reverse_becomes_minus1(self):
        """BEAM/MATSim require "-1", not "reverse"."""
        g = _single_edge(oneway="reverse")
        result = process_tags(g, _DEFAULT_CONFIG)
        assert _get_edge_attr(result, "oneway") == "-1"


# ---------------------------------------------------------------------------
# access / motor_vehicle standardization
# ---------------------------------------------------------------------------

class TestProcessTagsAccess:
    def test_access_yes_preserved(self):
        g = _single_edge(access="yes")
        result = process_tags(g, _DEFAULT_CONFIG)
        assert _get_edge_attr(result, "access") == "yes"

    def test_access_private_is_no(self):
        g = _single_edge(access="private")
        result = process_tags(g, _DEFAULT_CONFIG)
        assert _get_edge_attr(result, "access") == "no"

    def test_motor_vehicle_private_is_no(self):
        g = _single_edge(motor_vehicle="private")
        result = process_tags(g, _DEFAULT_CONFIG)
        assert _get_edge_attr(result, "motor_vehicle") == "no"


# ---------------------------------------------------------------------------
# hgv / mdv boolean flags from weight limits
# ---------------------------------------------------------------------------

class TestProcessTagsWeightFlags:
    def test_no_maxweight_hgv_and_mdv_default_true(self):
        """Without weight tags, hgv and mdv must both default to True (allowed)."""
        g = _single_edge()
        result = process_tags(g, _DEFAULT_CONFIG)
        assert _get_edge_attr(result, "hgv") == True
        assert _get_edge_attr(result, "mdv") == True

    def test_weight_below_mdv_max_restricts_mdv_and_hgv(self):
        """
        A maxweight below mdv_max (26000 lbs) must set both mdv=False and hgv=False.
        Weight "5t" = 5 metric tons = 11023 lbs < 26000 lbs.
        """
        g = _single_edge(maxweight="5t")
        result = process_tags(g, _DEFAULT_CONFIG)
        assert _get_edge_attr(result, "mdv") == False
        assert _get_edge_attr(result, "hgv") == False

    def test_weight_between_mdv_and_hdv_max_restricts_only_hgv(self):
        """
        A maxweight above mdv_max but below hdv_max restricts hgv but not mdv.
        "20t" = 44092 lbs — above 26000 (mdv_max) but below 80000 (hdv_max).
        """
        g = _single_edge(maxweight="20t")
        result = process_tags(g, _DEFAULT_CONFIG)
        assert _get_edge_attr(result, "mdv") == True
        assert _get_edge_attr(result, "hgv") == False

    def test_weight_above_hdv_max_allows_both(self):
        """
        A maxweight above hdv_max (80000 lbs) means both hgv and mdv are allowed.
        "50t" = 110231 lbs > 80000 lbs.
        """
        g = _single_edge(maxweight="50t")
        result = process_tags(g, _DEFAULT_CONFIG)
        assert _get_edge_attr(result, "hgv") == True
        assert _get_edge_attr(result, "mdv") == True

    def test_hgv_tag_false_restricts_hgv(self):
        """
        An explicit hgv=False tag must produce hgv=False regardless of maxweight.
        Regression: old ``if not value`` guard treated bool False as missing → True.
        """
        g = _single_edge(hgv=False)
        result = process_tags(g, _DEFAULT_CONFIG)
        assert _get_edge_attr(result, "hgv") == False

    def test_hgv_tag_true_allows_hgv(self):
        g = _single_edge(hgv=True)
        result = process_tags(g, _DEFAULT_CONFIG)
        assert _get_edge_attr(result, "hgv") == True

    def test_hgv_and_mdv_output_are_bool_dtype(self):
        """hgv and mdv columns must be bool, not object or int."""
        g = _single_edge(maxweight="5t")
        result = process_tags(g, _DEFAULT_CONFIG)
        _, edges = ox.graph_to_gdfs(result)
        assert edges["hgv"].dtype == "bool", f"hgv dtype is {edges['hgv'].dtype}"
        assert edges["mdv"].dtype == "bool", f"mdv dtype is {edges['mdv'].dtype}"

    def test_maxlength_restricts_hgv(self):
        """Any non-null maxlength tag must set hgv=False."""
        g = _single_edge(maxlength="10")
        result = process_tags(g, _DEFAULT_CONFIG)
        assert _get_edge_attr(result, "hgv") == False

    def test_maxweight_hgv_tag_overrides_maxweight(self):
        """
        maxweight:hgv tag must override maxweight for hgv classification.
        maxweight="50t" (above hdv_max) but maxweight:hgv="5t" (below mdv_max)
        → hgv must be False.
        """
        g = _single_edge(**{"maxweight": "50t", "maxweight:hgv": "5t"})
        result = process_tags(g, _DEFAULT_CONFIG)
        assert _get_edge_attr(result, "hgv") == False

    def test_bare_numeric_maxweight_treated_as_metric_tons(self):
        """
        OSM bare numeric weight "5" means 5 metric tons (11023 lbs).
        That is below mdv_max=26000 lbs, so mdv and hgv must both be False.
        Regression: bare numerics were previously treated as raw lbs (5 lbs).
        """
        g = _single_edge(maxweight="5")
        result = process_tags(g, _DEFAULT_CONFIG)
        assert _get_edge_attr(result, "mdv") == False
        assert _get_edge_attr(result, "hgv") == False


# ---------------------------------------------------------------------------
# maxspeed standardization
# ---------------------------------------------------------------------------

class TestProcessTagsMaxspeed:
    def test_maxspeed_kph_converted_to_mph_string(self):
        g = _single_edge(maxspeed="50")
        result = process_tags(g, _DEFAULT_CONFIG)
        val = _get_edge_attr(result, "maxspeed")
        assert isinstance(val, str)
        assert "mph" in val

    def test_maxspeed_explicit_mph_preserved(self):
        g = _single_edge(maxspeed="25 mph")
        result = process_tags(g, _DEFAULT_CONFIG)
        assert _get_edge_attr(result, "maxspeed") == "25 mph"


# ---------------------------------------------------------------------------
# Missing required columns raise ValueError
# ---------------------------------------------------------------------------


class TestProcessTagsMissingColumns:
    """
    oneway is required — missing it raises ValueError.
    motor_vehicle / maxspeed / access are optional — missing them logs a
    warning and fills with None so standardize_* returns safe defaults.
    """

    def test_missing_oneway_raises_valueerror(self):
        """oneway is fundamental to routing; its absence is a hard error."""
        g = _single_edge()
        nodes, edges = ox.graph_to_gdfs(g)
        edges = edges.drop(columns=["oneway"], errors="ignore")
        g_stripped = ox.graph_from_gdfs(nodes, edges)

        with pytest.raises(ValueError, match="oneway"):
            process_tags(g_stripped, _DEFAULT_CONFIG)

    @pytest.mark.parametrize("missing_col", ["motor_vehicle", "maxspeed", "access"])
    def test_missing_optional_column_warns_and_continues(self, missing_col, caplog):
        """Optional columns missing → warning logged, function completes successfully."""
        import logging

        g = _single_edge()
        nodes, edges = ox.graph_to_gdfs(g)
        edges = edges.drop(columns=[missing_col], errors="ignore")
        g_stripped = ox.graph_from_gdfs(nodes, edges)

        with caplog.at_level(logging.WARNING, logger="osm_chordify.osm.graph"):
            result = process_tags(g_stripped, _DEFAULT_CONFIG)

        assert result is not None, "process_tags must return a graph, not crash"
        assert any(missing_col in msg for msg in caplog.messages), (
            f"Expected a warning mentioning '{missing_col}'"
        )

    def test_missing_optional_columns_produce_correct_defaults(self):
        """
        motor_vehicle=None → "yes", access=None → "yes", maxspeed=None → None.
        These are the safe defaults from each standardize_* function.
        """
        g = _single_edge()
        nodes, edges = ox.graph_to_gdfs(g)
        edges = edges.drop(columns=["motor_vehicle", "access", "maxspeed"], errors="ignore")
        g_stripped = ox.graph_from_gdfs(nodes, edges)

        result = process_tags(g_stripped, _DEFAULT_CONFIG)
        assert _get_edge_attr(result, "motor_vehicle") == "yes"
        assert _get_edge_attr(result, "access") == "yes"