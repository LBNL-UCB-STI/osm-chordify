import importlib
import importlib.util
from pathlib import Path
from statistics import mean

import geopandas as gpd
import networkx as nx
import osmnx as ox
import pytest
from shapely.geometry import LineString, Point
from shapely.strtree import STRtree

from osm_chordify.osm.export import export_network
from osm_chordify.osm.graph import (
    _download_with_overpass_fallback,
    create_unique_edge_id,
    format_validation_summary,
    process_ferry_edges,
    process_tags,
    summarize_graph_validation,
    validate_graph_topology,
)
from osm_chordify.osm.simplify import (
    bool_all,
    first_valid_value,
    mean_maxspeed,
    median_lanes,
    min_numeric_or_string,
    yes_no_all,
)
from osm_chordify.utils.geo import project_graph


ROOT = Path(__file__).resolve().parents[1]
_TEST_OVERPASS_URLS = [
    "https://overpass-api.de/api",
    "https://overpass.kumi.systems/api",
    "https://overpass.private.coffee/api",
]


def _build_graph(node_coords: dict, edges: list):
    node_rows = [
        {"osmid": node_id, "x": lon, "y": lat, "geometry": Point(lon, lat)}
        for node_id, (lon, lat) in node_coords.items()
    ]
    nodes = gpd.GeoDataFrame(node_rows, geometry="geometry", crs="EPSG:4326").set_index("osmid")

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
                "oneway": "yes",
                "geometry": LineString([(lon1, lat1), (lon2, lat2)]),
            }
        )
    edges = gpd.GeoDataFrame(edge_rows, geometry="geometry", crs="EPSG:4326").set_index(
        ["u", "v", "key"]
    )
    return ox.graph_from_gdfs(nodes, edges)


def _load_module(name: str, relative_path: str):
    path = ROOT / relative_path
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec is not None and spec.loader is not None
    spec.loader.exec_module(module)
    return module


@pytest.fixture(autouse=True)
def _force_all_oneway(monkeypatch):
    # XML/PBF export in this project is intentionally exercised under the same
    # supported OSMnx mode used by the production pipeline.
    monkeypatch.setattr(ox.settings, "all_oneway", True)


def _configure_fast_overpass_test_settings():
    """Use shorter per-endpoint timeouts and quiet logging for integration tests."""
    ox.settings.use_cache = True
    ox.settings.log_console = False
    ox.settings.all_oneway = True
    ox.settings.requests_timeout = 20
    ox.settings.overpass_max_attempts = 1
    ox.settings.overpass_rate_limit = False


class _PbfStatsHandler:
    def __init__(self):
        self.node_ids = set()
        self.way_count = 0
        self.dangling_refs = []

    def node(self, node):
        self.node_ids.add(node.id)

    def way(self, way):
        self.way_count += 1
        for node_ref in way.nodes:
            if node_ref.ref not in self.node_ids:
                self.dangling_refs.append(node_ref.ref)


def _make_pbf_stats_handler(osmium):
    class StatsHandler(osmium.SimpleHandler, _PbfStatsHandler):
        def __init__(self):
            osmium.SimpleHandler.__init__(self)
            _PbfStatsHandler.__init__(self)

    return StatsHandler


def _finalize_graph_for_export(graph, *, utm_epsg: int, network_config: dict):
    graph_projected = project_graph(graph, to_crs=utm_epsg)
    graph_with_speeds = ox.add_edge_speeds(graph_projected)
    graph_processed = process_tags(graph_with_speeds, network_config)

    graph_consolidated = ox.consolidate_intersections(
        graph_processed,
        tolerance=network_config["tolerance"],
        rebuild_graph=True,
        dead_ends=True,
        reconnect_edges=True,
    )
    graph_consolidated.graph["consolidated"] = False
    for _, node_data in graph_consolidated.nodes(data=True):
        node_data.pop("cluster", None)
        node_data.pop("osmid_original", None)
    graph_consolidated = ox.consolidate_intersections(
        graph_consolidated,
        tolerance=1.5,
        rebuild_graph=True,
        dead_ends=False,
        reconnect_edges=True,
    )

    graph_simplified = ox.simplification.simplify_graph(
        graph_consolidated,
        edge_attrs_differ=["highway", "lanes", "maxspeed"],
        remove_rings=False,
        track_merged=True,
        edge_attr_aggs={
            "length": sum,
            "travel_time": sum,
            "hgv": bool_all,
            "mdv": bool_all,
            "lanes": median_lanes,
            "speed_kph": mean,
            "maxspeed": mean_maxspeed,
            "oneway": yes_no_all,
            "access": yes_no_all,
            "reversed": bool_all,
            "maxweight": min_numeric_or_string,
            "bridge": first_valid_value,
            "tunnel": first_valid_value,
            "foot": yes_no_all,
            "bicycle": yes_no_all,
            "sidewalk": first_valid_value,
            "cycleway": first_valid_value,
            "maxheight": min_numeric_or_string,
            "maxwidth": min_numeric_or_string,
            "motor_vehicle": yes_no_all,
        },
    )

    nodes, edges = ox.graph_to_gdfs(graph_simplified)
    edges["edge_id"] = [
        create_unique_edge_id(u, v, row["osmid"], k)
        for (u, v, k), row in edges.iterrows()
    ]
    graph_hashed = ox.graph_from_gdfs(nodes, edges)
    graph_wgs84 = project_graph(graph_hashed, to_latlong=True)
    graph_validated = validate_graph_topology(graph_wgs84)
    largest_component = max(nx.weakly_connected_components(graph_validated), key=len)
    return nx.MultiDiGraph(graph_validated.subgraph(largest_component).copy())


def _prepare_minimal_graph_for_export(graph):
    graph_projected = project_graph(graph, to_crs=32610)
    graph_with_speeds = ox.add_edge_speeds(graph_projected)
    graph_processed = process_tags(
        graph_with_speeds,
        config={
            "weight_limits": {
                "unit": "lbs",
                "mdv_max": 26000,
                "hdv_max": 80000,
            }
        },
    )
    nodes, edges = ox.graph_to_gdfs(graph_processed)
    edges["edge_id"] = [f"e_{u}_{v}_{k}" for u, v, k in edges.index]
    graph_hashed = ox.graph_from_gdfs(nodes, edges)
    graph_wgs84 = project_graph(graph_hashed, to_latlong=True)
    graph_validated = validate_graph_topology(graph_wgs84)
    largest_component = max(nx.weakly_connected_components(graph_validated), key=len)
    return nx.MultiDiGraph(graph_validated.subgraph(largest_component).copy())


def _assert_graph_integrity(graph, *, graph_name: str):
    summary = summarize_graph_validation(graph, close_threshold_m=0.5)
    print(format_validation_summary(summary, title=f"{graph_name} graph validation checks"))

    assert summary["self_loops"] == 0
    assert summary["isolated_nodes"] == 0
    assert summary["weakly_connected"], "Final graph is not weakly connected"
    assert summary["duplicate_xml_coordinate_pairs"] == 0, (
        "Found node pairs with duplicate XML-rounded coordinates: "
        f"{summary['duplicate_examples'][:5]}"
    )
    assert summary["close_node_pairs_lt_0_5m"] == 0, (
        "Found node pairs within 0.5m in final graph: "
        f"{summary['close_examples'][:5]}"
    )
    assert summary["missing_edge_id"] == 0
    assert summary["missing_length"] == 0
    assert summary["nonpositive_length"] == 0
    assert summary["missing_speed_kph"] == 0
    assert summary["invalid_oneway_values"] == 0
    assert summary["missing_geometry"] == 0


def _assert_pbf_integrity(pbf_path: Path, handler_cls):
    assert pbf_path.exists()
    assert pbf_path.stat().st_size > 0

    handler = handler_cls()
    handler.apply_file(str(pbf_path))

    assert len(handler.node_ids) > 0
    assert handler.way_count > 0
    assert handler.dangling_refs == []


def test_export_network_writes_geojson_without_external_cli(tmp_path):
    graph = _build_graph(
        node_coords={
            1: (-122.0, 37.0),
            2: (-122.001, 37.0),
            3: (-122.001, 37.001),
        },
        edges=[(1, 2, "e1", 10), (2, 3, "e2", 11)],
    )

    exported = export_network(
        graph,
        output_dir=tmp_path,
        name="sample",
        formats={"geojson"},
    )

    assert Path(exported["geojson"]).exists()
    edges = gpd.read_file(exported["geojson"])
    assert len(edges) == 2
    assert "edge_id" in edges.columns
    assert "pbf" not in exported


def test_export_network_pbf_requires_python_osmium(tmp_path, monkeypatch):
    graph = _build_graph(
        node_coords={1: (-122.0, 37.0), 2: (-122.001, 37.0)},
        edges=[(1, 2, "e1", 10)],
    )

    real_import_module = importlib.import_module

    def fake_import_module(name, package=None):
        if name == "osmium":
            raise ModuleNotFoundError("No module named 'osmium'")
        return real_import_module(name, package)

    monkeypatch.setattr(importlib, "import_module", fake_import_module)

    try:
        export_network(graph, output_dir=tmp_path, name="sample", formats={"pbf"})
    except ModuleNotFoundError as exc:
        assert "Install it with `pip install osmium`" in str(exc)
    else:
        raise AssertionError("Expected PBF export to fail cleanly when pyosmium is absent")


def test_export_network_pbf_uses_pyosmium_backend(tmp_path, monkeypatch):
    graph = _build_graph(
        node_coords={1: (-122.0, 37.0), 2: (-122.001, 37.0)},
        edges=[(1, 2, "e1", 10)],
    )

    class FakeWriter:
        def __init__(self, path):
            self.path = Path(path)
            self.items = []

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            self.path.write_text("\n".join(self.items))

        def add(self, obj):
            self.items.append(str(obj))

    class FakeOsmium:
        @staticmethod
        def SimpleWriter(path):
            return FakeWriter(path)

        @staticmethod
        def FileProcessor(path):
            return iter(["node:1", "way:1"])

    real_import_module = importlib.import_module

    def fake_import_module(name, package=None):
        if name == "osmium":
            return FakeOsmium
        return real_import_module(name, package)

    monkeypatch.setattr(importlib, "import_module", fake_import_module)

    exported = export_network(
        graph,
        output_dir=tmp_path,
        name="sample",
        formats={"pbf"},
    )

    pbf_path = Path(exported["pbf"])
    assert pbf_path.exists()
    assert "node:1" in pbf_path.read_text()


@pytest.mark.integration
def test_export_network_pbf_smoke_real_download(tmp_path):
    print("\n=== Generic Real-Download Smoke Test ===")
    osmium = importlib.import_module("osmium")
    StatsHandler = _make_pbf_stats_handler(osmium)

    _configure_fast_overpass_test_settings()
    network_config = {
        "osmnx_settings": {
            "overpass_url": _TEST_OVERPASS_URLS[0],
            "overpass_urls": _TEST_OVERPASS_URLS,
        }
    }

    graph = _download_with_overpass_fallback(
        lambda: ox.graph_from_point(
            (47.6062, -122.3321),
            dist=200,
            network_type="drive",
            simplify=False,
            retain_all=True,
        ),
        network_config,
        "real-download-test",
    )
    graph = _prepare_minimal_graph_for_export(graph)
    print("Generic real-download smoke checks:")
    print("- downloaded live OSM sample successfully")
    print("- processed graph for export successfully")
    print(f"- final nodes: {graph.number_of_nodes()}")
    print(f"- final edges: {graph.number_of_edges()}")

    exported = export_network(
        graph,
        output_dir=tmp_path,
        name="seattle-mini",
        formats={"pbf"},
    )

    _assert_pbf_integrity(Path(exported["pbf"]), StatsHandler)
    print("- exported and parsed .osm.pbf successfully")
    print("=== End Generic Real-Download Smoke Test ===")


@pytest.mark.integration
def test_export_network_pbf_validation_seattle_example_config(tmp_path):
    print("\n=== Seattle Example Validation Test ===")
    build_seattle = _load_module("build_seattle", "examples/build_seattle.py")

    osmium = importlib.import_module("osmium")
    StatsHandler = _make_pbf_stats_handler(osmium)

    for setting, value in build_seattle.osm_config["osmnx_settings"].items():
        setattr(ox.settings, setting, value)
    ox.settings.log_console = False
    ox.settings.requests_timeout = 20
    ox.settings.overpass_max_attempts = 1

    point = (47.6055, -122.3370)  # Seattle waterfront; captures ferry + road context
    graphs = []

    for layer_name, layer_cfg in build_seattle.osm_config["graph_layers"].items():
        network_type = "all" if layer_name == "ferry" else "drive"
        graph = ox.graph_from_point(
            point,
            dist=300,
            network_type=network_type,
            simplify=False,
            retain_all=True,
            truncate_by_edge=(layer_name != "ferry"),
            custom_filter=layer_cfg["custom_filter"],
        )
        if layer_name == "ferry":
            graph = process_ferry_edges(graph)
            if graph.number_of_edges() == 0:
                continue
        graphs.append(graph)

    assert graphs, "Seattle example sample download returned no graph layers"

    graph_combined = nx.compose_all(graphs)
    graph_final = _finalize_graph_for_export(
        graph_combined,
        utm_epsg=build_seattle.geo_config["utm_epsg"],
        network_config=build_seattle.osm_config,
    )
    _assert_graph_integrity(graph_final, graph_name="Seattle example sample")

    exported = export_network(
        graph_final,
        output_dir=tmp_path,
        name="seattle-example-mini",
        formats={"pbf"},
    )
    _assert_pbf_integrity(Path(exported["pbf"]), StatsHandler)
    print("=== End Seattle Example Validation Test ===")
