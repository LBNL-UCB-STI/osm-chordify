"""
Tests for graph.py: topology validation, unique edge IDs, and the OSM download
pipeline.

Unit tests use synthetic osmnx-compatible graphs (no network access).
Integration tests download a small real area and are marked with
``@pytest.mark.integration`` so they can be excluded from fast CI runs:

    pytest -m "not integration"   # fast, no network
    pytest -m integration         # slow, needs network
"""

import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path

import geopandas as gpd
import networkx as nx
import osmnx as ox
import pytest
from shapely.geometry import LineString, Point

from osm_chordify.osm.graph import (
    create_unique_edge_id,
    _download_with_overpass_fallback,
    download_and_prepare_osm_network,
    validate_graph_topology,
)
from osm_chordify.osm.xml import save_graph_xml

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_XML_PRECISION = 7  # must match precision= in xml.py _save_graph_xml


@pytest.fixture(autouse=True)
def _force_all_oneway(monkeypatch):
    # XML export tests should mirror the supported production setting.
    monkeypatch.setattr(ox.settings, "all_oneway", True)


def _build_graph(node_coords: dict, edges: list) -> nx.MultiDiGraph:
    """
    Build a minimal osmnx-compatible WGS-84 MultiDiGraph.

    Parameters
    ----------
    node_coords : dict
        ``{node_id: (lon, lat)}``
    edges : list of (u, v, edge_id, osmid)
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
    for u, v, edge_id, osmid in edges:
        lon1, lat1 = node_coords[u]
        lon2, lat2 = node_coords[v]
        edge_rows.append(
            {
                "u": u,
                "v": v,
                "key": 0,
                "osmid": osmid,
                "edge_id": edge_id,
                "length": 100.0,
                "highway": "residential",
                "geometry": LineString([(lon1, lat1), (lon2, lat2)]),
            }
        )
    edges_gdf = (
        gpd.GeoDataFrame(edge_rows, geometry="geometry", crs="EPSG:4326")
        .set_index(["u", "v", "key"])
    )

    return ox.graph_from_gdfs(nodes_gdf, edges_gdf)


def _duplicate_coords_at_precision(g: nx.MultiDiGraph, precision: int = _XML_PRECISION) -> list:
    """
    Return list of (node_id_a, node_id_b) pairs that share the same
    rounded (lon, lat) — i.e. would become the same vertex in R5.
    """
    nodes_gdf, _ = ox.graph_to_gdfs(g)
    coord_map: dict = {}
    dupes = []
    for nid, row in nodes_gdf.iterrows():
        key = (round(row["x"], precision), round(row["y"], precision))
        if key in coord_map:
            dupes.append((coord_map[key], nid))
        else:
            coord_map[key] = nid
    return dupes


def _assert_unique_node_coordinates(
    g: nx.MultiDiGraph,
    *,
    precision: int = _XML_PRECISION,
    context: str = "graph",
) -> None:
    """Assert that no node IDs collapse to the same rounded coordinate pair."""
    dupes = _duplicate_coords_at_precision(g, precision=precision)
    assert dupes == [], (
        f"{context} contains {len(dupes)} node pair(s) with duplicate XML "
        f"coordinates. Examples: {dupes[:5]}"
    )


def _nodes_within_distance(g: nx.MultiDiGraph, threshold_m: float) -> list:
    """Return (node_id_a, node_id_b) pairs that are < threshold_m apart."""
    from shapely.strtree import STRtree

    g_proj = ox.project_graph(g)
    nodes_proj, _ = ox.graph_to_gdfs(g_proj)
    node_ids = list(nodes_proj.index)
    points = [Point(row.geometry.x, row.geometry.y) for _, row in nodes_proj.iterrows()]
    tree = STRtree(points)

    close = []
    for i, pt in enumerate(points):
        for j in tree.query(pt.buffer(threshold_m)):
            if j > i:
                close.append((node_ids[i], node_ids[j]))
    return close


# ---------------------------------------------------------------------------
# Unit tests — validate_graph_topology
# ---------------------------------------------------------------------------


class TestValidateGraphTopology:
    def test_removes_networkx_self_loops(self):
        """Edges where u == v must be removed."""
        g = _build_graph(
            node_coords={1: (-122.0, 37.0), 2: (-122.001, 37.001)},
            edges=[
                (1, 2, "e1", 10),
                (1, 1, "e_loop", 99),  # self-loop
            ],
        )
        assert any(u == v for u, v, _ in g.edges(keys=True))

        g_fixed = validate_graph_topology(g)

        assert not any(u == v for u, v, _ in g_fixed.edges(keys=True))

    def test_retains_protected_major_self_loops(self):
        """Long major-road self-loops are preserved for manual review."""
        g = _build_graph(
            node_coords={1: (-122.0, 37.0), 2: (-122.001, 37.001)},
            edges=[
                (1, 2, "e1", 10),
                (1, 1, "major_loop", 99),
            ],
        )
        nodes, edges = ox.graph_to_gdfs(g)
        edges.loc[(1, 1, 0), "length"] = 150.0
        edges.loc[(1, 1, 0), "highway"] = "motorway"
        g = ox.graph_from_gdfs(nodes, edges)

        g_fixed = validate_graph_topology(g)

        assert any(u == v for u, v, _ in g_fixed.edges(keys=True))

    def test_removes_isolated_nodes_and_raises_on_empty_result(self):
        """
        When cleanup leaves zero nodes the function must raise ValueError with a
        message explaining how many edges/nodes were removed vs the original counts.
        """
        g = _build_graph(
            node_coords={1: (-122.0, 37.0), 2: (-122.001, 37.001)},
            edges=[(1, 1, "loop_only", 99)],  # node 2 has no edges at all
        )
        with pytest.raises(ValueError, match="graph is empty after cleanup"):
            validate_graph_topology(g)

    def test_removes_isolated_nodes_from_non_empty_graph(self):
        """
        Isolated nodes are dropped while connected nodes are kept.
        The graph must not be empty for validate_graph_topology to proceed.
        """
        g = _build_graph(
            node_coords={
                1: (-122.0, 37.0),      # isolated — no edges
                2: (-122.001, 37.001),  # connected to 3 → survives
                3: (-122.002, 37.002),  # connected to 2 → survives
            },
            edges=[(2, 3, "e1", 1)],
        )
        g_fixed = validate_graph_topology(g)
        assert 1 not in g_fixed.nodes(), "Isolated node 1 should be removed"
        assert g_fixed.number_of_nodes() == 2

    def test_no_changes_to_clean_graph(self):
        """A clean graph (no self-loops, no isolates) is returned unchanged."""
        g = _build_graph(
            node_coords={
                1: (-122.0, 37.0),
                2: (-122.001, 37.0),
                3: (-122.001, 37.001),
            },
            edges=[(1, 2, "e1", 1), (2, 3, "e2", 2), (3, 1, "e3", 3)],
        )
        g_fixed = validate_graph_topology(g)
        assert g_fixed.number_of_nodes() == 3
        assert g_fixed.number_of_edges() == 3

    def test_warns_duplicate_edge_ids(self, caplog):
        """Duplicate edge_id values must trigger a warning."""
        import logging

        g = _build_graph(
            node_coords={1: (-122.0, 37.0), 2: (-122.001, 37.0), 3: (-122.002, 37.0)},
            edges=[
                (1, 2, "SAME_ID", 10),
                (2, 3, "SAME_ID", 11),  # duplicate edge_id
            ],
        )
        with caplog.at_level(logging.WARNING, logger="osm_chordify.osm.graph"):
            validate_graph_topology(g)

        assert any("duplicate" in msg.lower() for msg in caplog.messages)

    def test_warns_on_close_nodes(self, caplog):
        """Two nodes within 1 m of each other (but distinct IDs) must trigger a
        proximity warning — these are potential R5 self-loop candidates."""
        import logging

        # ~0.9 m separation at 37 °N: Δlat ≈ 8e-6° ≈ 0.89 m
        g = _build_graph(
            node_coords={
                1: (-122.0, 37.0),
                2: (-122.0000090, 37.0000080),  # ≈ 0.9 m away
                3: (-122.005, 37.005),           # far away, no warning
            },
            edges=[(1, 2, "e1", 10), (2, 3, "e2", 11)],
        )
        with caplog.at_level(logging.WARNING, logger="osm_chordify.osm.graph"):
            validate_graph_topology(g)

        assert any("close" in msg.lower() or "within" in msg.lower() for msg in caplog.messages)

    def test_no_proximity_warning_for_well_separated_nodes(self, caplog):
        """Nodes clearly more than 1 m apart must not trigger a proximity warning."""
        import logging

        g = _build_graph(
            node_coords={
                1: (-122.0, 37.0),
                2: (-122.01, 37.01),  # ~1.4 km away
                3: (-122.02, 37.02),
            },
            edges=[(1, 2, "e1", 10), (2, 3, "e2", 11)],
        )
        with caplog.at_level(logging.WARNING, logger="osm_chordify.osm.graph"):
            validate_graph_topology(g)

        assert not any(
            "close" in msg.lower() or "within" in msg.lower() for msg in caplog.messages
        )


class TestOverpassFallback:
    def test_tries_next_endpoint_after_failure(self, monkeypatch):
        config = {
            "osmnx_settings": {
                "overpass_url": "https://primary.example/api",
                "overpass_urls": [
                    "https://primary.example/api",
                    "https://secondary.example/api",
                ],
            }
        }
        attempts = []

        def download_fn():
            attempts.append(ox.settings.overpass_url)
            if len(attempts) == 1:
                raise RuntimeError("primary failed")
            return "ok"

        result = _download_with_overpass_fallback(download_fn, config, "test-layer")

        assert result == "ok"
        assert attempts == [
            "https://primary.example/api",
            "https://secondary.example/api",
        ]

    def test_raises_after_all_endpoints_fail(self):
        config = {
            "osmnx_settings": {
                "overpass_url": "https://primary.example/api",
                "overpass_urls": [
                    "https://primary.example/api",
                    "https://secondary.example/api",
                ],
            }
        }

        def download_fn():
            raise RuntimeError("boom")

        with pytest.raises(RuntimeError, match="All Overpass endpoints failed"):
            _download_with_overpass_fallback(download_fn, config, "test-layer")


class TestNetworkCoordinateIntegrity:
    """Coordinate-integrity checks on the output of validate_graph_topology()."""

    def test_validated_graph_has_unique_coordinates_at_xml_precision(self):
        """
        A validated graph must not contain two node IDs that collapse to the same
        rounded XML coordinate pair, since that would become a self-loop in R5.
        """
        g = _build_graph(
            node_coords={
                1: (-122.0000000, 37.0000000),
                2: (-122.0005000, 37.0005000),
                3: (-122.0010000, 37.0010000),
            },
            edges=[(1, 2, "e1", 10), (2, 3, "e2", 11)],
        )

        g_validated = validate_graph_topology(g)
        _assert_unique_node_coordinates(g_validated, context="Validated graph")


# ---------------------------------------------------------------------------
# Unit tests — R5 self-loop regression (coordinate duplicates)
# ---------------------------------------------------------------------------


class TestCoordinateDuplicateRegression:
    """Regression tests for the R5 self-loop bug.

    ``consolidate_intersections(reconnect_edges=True)`` can place connector nodes
    at the cluster centroid, creating two distinct node IDs with sub-millimetre
    separation.  After reprojection + rounding in the XML export, both share the
    same coordinate, so R5 maps them to the same vertex and any edge between them
    becomes a self-loop.  ``validate_graph_topology`` removes unprotected
    self-loops to prevent this.
    """

    def test_detects_nodes_with_identical_xml_coordinates(self):
        """
        Two nodes whose coordinates round identically at the XML precision must be
        flagged as duplicates by _duplicate_coords_at_precision().

        At precision=7 the rounding threshold is ~5e-8 ° (≈5 mm).  The nodes
        below differ by 2e-8 °, so they collapse to the same value when rounded.
        """
        g = _build_graph(
            node_coords={
                1: (-122.12345600, 37.12345600),
                2: (-122.12345602, 37.12345602),  # Δ=2e-8°, rounds to same at p=7
                3: (-122.20000000, 37.20000000),
            },
            edges=[(1, 2, "e1", 10), (2, 3, "e2", 11)],
        )
        dupes = _duplicate_coords_at_precision(g)
        assert len(dupes) > 0, "Expected to detect coordinate-identical nodes"

    def test_xml_export_has_no_duplicate_node_coordinates(self):
        """
        A clean graph (no near-duplicate nodes) must produce an XML file where
        every <node> has a unique (lat, lon) pair — the invariant R5 requires.
        """
        g = _build_graph(
            node_coords={
                1: (-122.0, 37.0),
                2: (-122.001, 37.0),
                3: (-122.001, 37.001),
            },
            edges=[(1, 2, "e1", 1), (2, 3, "e2", 2)],
        )
        with tempfile.NamedTemporaryFile(suffix=".osm", delete=False) as f:
            tmp = Path(f.name)

        try:
            save_graph_xml(g, filepath=tmp)
            tree = ET.parse(tmp)
            root = tree.getroot()

            coords = [
                (node.attrib["lat"], node.attrib["lon"])
                for node in root.findall("node")
            ]
            assert len(coords) == len(set(coords)), (
                "Duplicate (lat, lon) pairs found in XML — R5 would produce self-loops"
            )
        finally:
            tmp.unlink(missing_ok=True)

    def test_xml_way_nd_refs_all_reference_existing_nodes(self):
        """
        Every <nd ref=...> inside a <way> must point to a <node id=...> that
        actually exists in the XML — R5 would crash on dangling refs.
        """
        g = _build_graph(
            node_coords={
                1: (-122.0, 37.0),
                2: (-122.001, 37.0),
                3: (-122.001, 37.001),
            },
            edges=[(1, 2, "e1", 1), (2, 3, "e2", 2)],
        )
        with tempfile.NamedTemporaryFile(suffix=".osm", delete=False) as f:
            tmp = Path(f.name)

        try:
            save_graph_xml(g, filepath=tmp)
            tree = ET.parse(tmp)
            root = tree.getroot()

            node_ids = {node.attrib["id"] for node in root.findall("node")}
            dangling = [
                ref
                for way in root.findall("way")
                for nd in way.findall("nd")
                if (ref := nd.attrib["ref"]) not in node_ids
            ]
            assert not dangling, f"Dangling nd refs in XML: {dangling[:10]}"
        finally:
            tmp.unlink(missing_ok=True)


# ---------------------------------------------------------------------------
# Unit tests — XML export correctness
# ---------------------------------------------------------------------------


class TestXmlExport:
    """
    Regression tests for xml.py fixes that aren't covered by the coordinate
    duplicate tests above.
    """

    def test_oneway_bool_serialized_as_yes_no(self):
        """
        Regression: oneway was written as "True"/"False" strings in XML because
        a Series.replace call wasn't applied after astype(str).  R5/BEAM expect
        the OSM-standard "yes"/"no".
        """
        g = _build_graph(
            node_coords={
                1: (-122.0, 37.0),
                2: (-122.001, 37.0),
                3: (-122.002, 37.0),
            },
            edges=[(1, 2, "e1", 1), (2, 3, "e2", 2)],
        )
        # Force oneway column to Python booleans (the problematic case)
        nodes, edges = ox.graph_to_gdfs(g)
        edges["oneway"] = [True, False]
        g_with_bool_oneway = ox.graph_from_gdfs(nodes, edges)

        with tempfile.NamedTemporaryFile(suffix=".osm", delete=False) as f:
            tmp = Path(f.name)
        try:
            save_graph_xml(g_with_bool_oneway, filepath=tmp)
            root = ET.parse(tmp).getroot()
            oneway_values = {
                tag.attrib["v"]
                for way in root.findall("way")
                for tag in way.findall("tag")
                if tag.attrib["k"] == "oneway"
            }
            assert "True" not in oneway_values, "Raw 'True' written to XML — bool replace broken"
            assert "False" not in oneway_values, "Raw 'False' written to XML — bool replace broken"
            assert oneway_values <= {"yes", "no", "-1"}, (
                f"Unexpected oneway values in XML: {oneway_values}"
            )
        finally:
            tmp.unlink(missing_ok=True)

    def test_list_valued_osmid_does_not_appear_in_xml(self):
        """
        Regression (D2): simplify_graph(track_merged=True) produces list-valued
        osmid columns.  export_network must reduce them to scalars before writing
        XML; otherwise the attribute value is serialized as "[123, 456]".
        """
        g = _build_graph(
            node_coords={
                1: (-122.0, 37.0),
                2: (-122.001, 37.0),
                3: (-122.002, 37.0),
            },
            edges=[(1, 2, "e1", 1), (2, 3, "e2", 2)],
        )
        # Simulate the merged-edge list that simplify_graph produces
        nodes, edges = ox.graph_to_gdfs(g)
        edges["osmid"] = [[100, 200], [300, 400]]
        g_with_list_osmid = ox.graph_from_gdfs(nodes, edges)

        # Replicate the normalization from export_network (Option A fix)
        _, edges_norm = ox.graph_to_gdfs(g_with_list_osmid)
        for col in edges_norm.columns:
            if edges_norm[col].apply(lambda x: isinstance(x, list)).any():
                edges_norm[col] = edges_norm[col].apply(
                    lambda x: min(x) if isinstance(x, list) else x
                )
        g_normalized = ox.graph_from_gdfs(nodes, edges_norm)

        with tempfile.NamedTemporaryFile(suffix=".osm", delete=False) as f:
            tmp = Path(f.name)
        try:
            save_graph_xml(g_normalized, filepath=tmp)
            content = tmp.read_text()
            assert "[" not in content, (
                "List bracket '[' found in XML — osmid list was not normalized before export"
            )
        finally:
            tmp.unlink(missing_ok=True)

    def test_list_min_is_deterministic(self):
        """min(list) on the same list always yields the same scalar."""
        values = [300, 100, 200]
        assert min(values) == 100
        assert min(sorted(values)) == 100


# ---------------------------------------------------------------------------
# Unit tests — edge_id assignment after simplify_graph
# ---------------------------------------------------------------------------


class TestEdgeIdAssignment:
    """
    Regression: graph.py used `row['u_original']` and `row['v_original']` when
    assigning edge_ids after simplify_graph.  After graph_to_gdfs(), u/v live in
    the MultiIndex (row.name), NOT as regular columns — so those accesses raised
    KeyError and the pipeline crashed before producing any output.
    """

    def test_create_unique_edge_id_uses_index_not_column(self):
        """
        Simulate what download_and_prepare_osm_network does after simplify_graph:
        iterate edges from graph_to_gdfs, get u/v from the MultiIndex, and call
        create_unique_edge_id.  Must not raise KeyError.
        """
        g = _build_graph(
            node_coords={
                1: (-122.0, 37.0),
                2: (-122.001, 37.0),
                3: (-122.002, 37.0),
            },
            edges=[(1, 2, "e1", 10), (2, 3, "e2", 11)],
        )
        _, edges = ox.graph_to_gdfs(g)

        # This is the fixed pattern: u/v from index, not from columns
        try:
            edge_ids = [
                create_unique_edge_id(u, v, row["osmid"], k)
                for (u, v, k), row in edges.iterrows()
            ]
        except KeyError as exc:
            pytest.fail(
                f"create_unique_edge_id raised KeyError {exc!r} — "
                "u/v must be read from the MultiIndex, not from row columns"
            )

        assert len(edge_ids) == 2
        assert all(isinstance(eid, str) and len(eid) == 12 for eid in edge_ids)

    def test_u_and_v_are_not_regular_columns_after_graph_to_gdfs(self):
        """
        Documents the root cause: after graph_to_gdfs, 'u' and 'v' are part of
        the MultiIndex, not regular DataFrame columns.  row['u'] would KeyError.
        """
        g = _build_graph(
            node_coords={1: (-122.0, 37.0), 2: (-122.001, 37.0)},
            edges=[(1, 2, "e1", 10)],
        )
        _, edges = ox.graph_to_gdfs(g)
        assert "u" not in edges.columns, (
            "'u' appeared as a regular column — check if osmnx changed its graph_to_gdfs API"
        )
        assert "v" not in edges.columns, (
            "'v' appeared as a regular column — check if osmnx changed its graph_to_gdfs API"
        )
        # But they ARE accessible from the index
        u, v, k = next(iter(edges.index))
        assert u == 1
        assert v == 2


# ---------------------------------------------------------------------------
# Unit tests — create_unique_edge_id
# ---------------------------------------------------------------------------


class TestCreateUniqueEdgeId:
    def test_deterministic(self):
        assert create_unique_edge_id(1, 2, 999) == create_unique_edge_id(1, 2, 999)

    def test_differs_on_u(self):
        assert create_unique_edge_id(1, 2, 999) != create_unique_edge_id(3, 2, 999)

    def test_differs_on_v(self):
        assert create_unique_edge_id(1, 2, 999) != create_unique_edge_id(1, 3, 999)

    def test_differs_on_osmid(self):
        assert create_unique_edge_id(1, 2, 999) != create_unique_edge_id(1, 2, 888)

    def test_handles_list_osmid(self):
        result = create_unique_edge_id(1, 2, [100, 200])
        assert isinstance(result, str) and len(result) == 12

    def test_includes_key_when_provided(self):
        assert create_unique_edge_id(1, 2, 999, k=0) != create_unique_edge_id(1, 2, 999, k=1)

    def test_output_is_12_char_hex(self):
        result = create_unique_edge_id(1, 2, 999)
        assert len(result) == 12
        assert all(c in "0123456789abcdef" for c in result)


# ---------------------------------------------------------------------------
# Integration tests — require network access
# ---------------------------------------------------------------------------

_INTEGRATION_AREA = (37.7749, -122.4194)   # SF downtown
_INTEGRATION_DIST = 400                     # metres radius


@pytest.mark.integration
class TestIntegrationDownloadPipeline:
    """Download a small real-world area and verify the full pipeline produces a
    network that is safe for R5/BEAM."""

    @pytest.fixture(scope="class")
    def processed_graph(self):
        """Run the real download + consolidation + validation pipeline on a small
        area.  Shared across all tests in this class to avoid multiple downloads."""
        import osmnx as ox

        G = ox.graph_from_point(
            _INTEGRATION_AREA,
            dist=_INTEGRATION_DIST,
            network_type="drive",
            simplify=False,
            retain_all=True,
        )
        utm_epsg = ox.projection.project_gdf(ox.graph_to_gdfs(G)[0]).crs.to_epsg()
        G = ox.project_graph(G, to_crs=f"EPSG:{utm_epsg}")
        G = ox.add_edge_speeds(G)

        g_consolidated = ox.consolidate_intersections(
            G, tolerance=10, rebuild_graph=True, dead_ends=True, reconnect_edges=True
        )
        g_wgs84 = ox.project_graph(g_consolidated, to_latlong=True)

        nodes, edges = ox.graph_to_gdfs(g_wgs84)
        edges["edge_id"] = [f"e_{u}_{v}_{k}" for u, v, k in edges.index]
        return validate_graph_topology(ox.graph_from_gdfs(nodes, edges))

    # --- topology -----------------------------------------------------------

    def test_no_self_loops(self, processed_graph):
        loops = list(nx.selfloop_edges(processed_graph))
        assert loops == [], f"Self-loops found: {loops[:5]}"

    def test_no_isolated_nodes(self, processed_graph):
        isolates = list(nx.isolates(processed_graph))
        assert isolates == [], f"Isolated nodes found: {isolates[:5]}"

    def test_graph_is_connected(self, processed_graph):
        assert nx.is_weakly_connected(processed_graph), "Graph is not weakly connected"

    # --- R5 self-loop regression --------------------------------------------

    def test_no_duplicate_node_coordinates_at_xml_precision(self, processed_graph):
        """No node pairs should round to the same (lon, lat) in the exported OSM
        XML — that is the condition that triggers R5 self-loop errors.
        """
        dupes = _duplicate_coords_at_precision(processed_graph)
        _assert_unique_node_coordinates(
            processed_graph,
            context="Processed graph",
        )

    def test_no_suspiciously_close_node_pairs(self, processed_graph):
        """No node pairs should remain within 0.1 m of each other after validation.

        Sub-mm artefacts created by ``reconnect_edges=True`` in consolidation
        become self-loops during simplification.  ``validate_graph_topology``
        removes those self-loops.  This test verifies no residual sub-0.1 m
        duplicates survive into the final graph.
        """
        close = _nodes_within_distance(processed_graph, threshold_m=0.1)
        assert close == [], (
            f"Found {len(close)} node pair(s) within 0.1 m — "
            f"possible residual consolidation artifact. Examples: {close[:5]}"
        )

    # --- edge quality -------------------------------------------------------

    def test_all_edges_reference_valid_nodes(self, processed_graph):
        node_set = set(processed_graph.nodes())
        bad = [
            (u, v)
            for u, v, _ in processed_graph.edges(keys=True)
            if u not in node_set or v not in node_set
        ]
        assert bad == [], f"Edges reference missing nodes: {bad[:5]}"

    def test_nodes_have_required_attributes(self, processed_graph):
        for nid, data in processed_graph.nodes(data=True):
            assert "x" in data, f"Node {nid} missing 'x'"
            assert "y" in data, f"Node {nid} missing 'y'"
            assert -180 <= data["x"] <= 180, f"Node {nid} lon out of range: {data['x']}"
            assert -90 <= data["y"] <= 90, f"Node {nid} lat out of range: {data['y']}"

    def test_edges_have_positive_length(self, processed_graph):
        _, edges = ox.graph_to_gdfs(processed_graph)
        assert (edges["length"] > 0).all(), "Some edges have non-positive length"

    # --- XML export ---------------------------------------------------------

    def test_xml_export_has_unique_node_coordinates(self, processed_graph):
        """
        The exported OSM XML must never have two <node> elements at the same
        (lat, lon) — that is the exact condition that makes R5 create self-loops.
        """
        with tempfile.NamedTemporaryFile(suffix=".osm", delete=False) as f:
            tmp = Path(f.name)
        try:
            save_graph_xml(processed_graph, filepath=tmp)
            root = ET.parse(tmp).getroot()
            coords = [
                (node.attrib["lat"], node.attrib["lon"])
                for node in root.findall("node")
            ]
            assert len(coords) == len(set(coords)), (
                "Duplicate coordinates in exported XML — R5 would produce self-loops"
            )
        finally:
            tmp.unlink(missing_ok=True)

    def test_xml_export_way_nd_refs_valid(self, processed_graph):
        """All <nd ref> values must point to a <node id> present in the XML."""
        with tempfile.NamedTemporaryFile(suffix=".osm", delete=False) as f:
            tmp = Path(f.name)
        try:
            save_graph_xml(processed_graph, filepath=tmp)
            root = ET.parse(tmp).getroot()
            node_ids = {n.attrib["id"] for n in root.findall("node")}
            dangling = [
                ref
                for way in root.findall("way")
                for nd in way.findall("nd")
                if (ref := nd.attrib["ref"]) not in node_ids
            ]
            assert not dangling, f"Dangling nd refs: {dangling[:10]}"
        finally:
            tmp.unlink(missing_ok=True)

