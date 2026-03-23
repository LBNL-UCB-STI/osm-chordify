"""Compare validation and diagnostics metrics across two built OSM PBF artifacts."""

import argparse
import sys
from pathlib import Path

EXAMPLES_DIR = Path(__file__).resolve().parent
if str(EXAMPLES_DIR) not in sys.path:
    sys.path.insert(0, str(EXAMPLES_DIR))

from _bootstrap import bootstrap_example_paths

bootstrap_example_paths(__file__)

from common import (
    build_validation_summary,
    format_validation_summary,
    infer_sidecar_paths,
    load_built_graph,
)
from osm_chordify.main import collect_osm_diagnostics, format_osm_diagnostics_summary


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("pbf_a", help="Path to the first .osm.pbf file.")
    parser.add_argument("pbf_b", help="Path to the second .osm.pbf file.")
    parser.add_argument("--epsg-utm-a", type=int, required=True, help="UTM EPSG for the first PBF.")
    parser.add_argument("--epsg-utm-b", type=int, default=None, help="UTM EPSG for the second PBF.")
    parser.add_argument("--graph-a", default=None, help="Optional .pkl or .graphml for the first artifact.")
    parser.add_argument("--graph-b", default=None, help="Optional .pkl or .graphml for the second artifact.")
    parser.add_argument("--osm-xml-a", default=None, help="Optional .osm XML for the first artifact.")
    parser.add_argument("--osm-xml-b", default=None, help="Optional .osm XML for the second artifact.")
    return parser.parse_args()


def _resolve_sidecars(pbf_path: Path, graph_path: str | None, osm_xml_path: str | None):
    graph = Path(graph_path).expanduser().resolve() if graph_path else None
    osm_xml = Path(osm_xml_path).expanduser().resolve() if osm_xml_path else None
    if graph is None or osm_xml is None:
        default_graph, default_osm = infer_sidecar_paths(pbf_path)
        graph = graph or default_graph
        osm_xml = osm_xml or default_osm
    return graph, osm_xml


def _print_metric_deltas(label: str, summary_a: dict, summary_b: dict, keys: list[str]):
    print(f"{label}:")
    differences = 0
    for key in keys:
        a = summary_a.get(key)
        b = summary_b.get(key)
        if a != b:
            print(f"- {key}: A={a} | B={b}")
            differences += 1
    if differences == 0:
        print("- no differences detected")


def main():
    args = parse_args()
    pbf_a = Path(args.pbf_a).expanduser().resolve()
    pbf_b = Path(args.pbf_b).expanduser().resolve()
    epsg_b = args.epsg_utm_b or args.epsg_utm_a

    graph_a, osm_xml_a = _resolve_sidecars(pbf_a, args.graph_a, args.osm_xml_a)
    graph_b, osm_xml_b = _resolve_sidecars(pbf_b, args.graph_b, args.osm_xml_b)

    build_summary_a = None
    build_summary_b = None
    if graph_a and graph_a.exists():
        print("=== Artifact A Built Graph Validation ===")
        build_summary_a = build_validation_summary(
            load_built_graph(graph_a),
            osm_path=str(osm_xml_a) if osm_xml_a and osm_xml_a.exists() else None,
        )
        print(format_validation_summary(build_summary_a, title="Artifact A validation checks"))
    else:
        print("Artifact A: no sibling built graph (.pkl or .graphml) found; skipping build-graph validation.")

    if graph_b and graph_b.exists():
        print("=== Artifact B Built Graph Validation ===")
        build_summary_b = build_validation_summary(
            load_built_graph(graph_b),
            osm_path=str(osm_xml_b) if osm_xml_b and osm_xml_b.exists() else None,
        )
        print(format_validation_summary(build_summary_b, title="Artifact B validation checks"))
    else:
        print("Artifact B: no sibling built graph (.pkl or .graphml) found; skipping build-graph validation.")

    print("=== Artifact A OSM PBF Diagnostics ===")
    pbf_summary_a = collect_osm_diagnostics(str(pbf_a), args.epsg_utm_a)
    print(format_osm_diagnostics_summary(pbf_summary_a, title="Artifact A OSM PBF diagnostics"))

    print("=== Artifact B OSM PBF Diagnostics ===")
    pbf_summary_b = collect_osm_diagnostics(str(pbf_b), epsg_b)
    print(format_osm_diagnostics_summary(pbf_summary_b, title="Artifact B OSM PBF diagnostics"))

    print("=== Metric Differences ===")
    if build_summary_a and build_summary_b:
        _print_metric_deltas(
            "Built validation deltas",
            build_summary_a,
            build_summary_b,
            [
                "nodes",
                "edges",
                "duplicate_xml_coordinate_pairs",
                "close_node_pairs_lt_0_5m",
                "short_links_lt_10m",
                "very_short_links_lt_5m",
                "long_links_gt_10km",
                "speed_kph_min",
                "speed_kph_max",
                "xml_duplicate_coordinate_pairs",
                "xml_dangling_nd_refs",
            ],
        )
        if build_summary_a.get("highway_type_counts") != build_summary_b.get("highway_type_counts"):
            print("- highway_type_counts differ")

    _print_metric_deltas(
        "PBF diagnostic deltas",
        pbf_summary_a,
        pbf_summary_b,
        [
            "retrieved_nodes",
            "retrieved_edges",
            "connected_components",
            "graph_nodes",
            "graph_edges",
            "largest_component_nodes",
            "largest_component_edges",
            "short_links_lt_15m",
            "long_links_gt_10km",
            "short_links_min_m",
            "short_links_max_m",
            "short_links_mean_m",
            "short_links_median_m",
            "links_le_500m_count",
            "links_le_500m_min_m",
            "links_le_500m_max_m",
            "links_le_500m_mean_m",
            "links_le_500m_median_m",
        ],
    )
    if pbf_summary_a.get("top_component_sizes") != pbf_summary_b.get("top_component_sizes"):
        print(f"- top_component_sizes: A={pbf_summary_a.get('top_component_sizes')} | B={pbf_summary_b.get('top_component_sizes')}")


if __name__ == "__main__":
    main()
