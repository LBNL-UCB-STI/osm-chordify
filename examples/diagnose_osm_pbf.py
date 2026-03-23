"""Run validation and diagnostics for a built OSM PBF artifact."""

import argparse
import os
import sys
from pathlib import Path

EXAMPLES_DIR = Path(__file__).resolve().parent
if str(EXAMPLES_DIR) not in sys.path:
    sys.path.insert(0, str(EXAMPLES_DIR))

from _bootstrap import bootstrap_example_paths

bootstrap_example_paths(__file__)

from common import infer_sidecar_paths, load_built_graph, validate_built_network
from osm_chordify import diagnose_osm


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("pbf_path", help="Path to the .osm.pbf file to diagnose.")
    parser.add_argument(
        "--epsg-utm",
        type=int,
        required=True,
        help="UTM EPSG code used for projected length diagnostics.",
    )
    parser.add_argument(
        "--graph-path",
        default=None,
        help="Optional path to a sibling .pkl or .graphml built graph for validation metrics.",
    )
    parser.add_argument(
        "--osm-xml",
        default=None,
        help="Optional path to a sibling .osm XML file for XML-level validation metrics.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    pbf_path = Path(args.pbf_path).expanduser().resolve()
    graph_path = Path(args.graph_path).expanduser().resolve() if args.graph_path else None
    osm_xml_path = Path(args.osm_xml).expanduser().resolve() if args.osm_xml else None

    if graph_path is None or osm_xml_path is None:
        default_graph_path, default_osm_path = infer_sidecar_paths(pbf_path)
        graph_path = graph_path or default_graph_path
        osm_xml_path = osm_xml_path or default_osm_path

    if graph_path and graph_path.exists():
        print("=== Built Graph Validation ===")
        graph = load_built_graph(graph_path)
        validate_built_network(
            graph,
            osm_path=str(osm_xml_path) if osm_xml_path and osm_xml_path.exists() else None,
        )
        print("=== End Built Graph Validation ===")
    else:
        print("No sibling built graph (.pkl or .graphml) found; skipping build-graph validation.")

    if os.environ.get("OSM_CHORDIFY_SKIP_DIAGNOSE") != "1":
        print("=== OSM PBF Diagnostics ===")
        diagnose_osm(str(pbf_path), args.epsg_utm)
        print("=== End OSM PBF Diagnostics ===")


if __name__ == "__main__":
    main()
