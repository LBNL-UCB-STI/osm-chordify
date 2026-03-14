"""Shared configuration helpers for runnable example scripts."""

import argparse
import json
import os
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

EXAMPLES_DIR = Path(__file__).resolve().parent
if str(EXAMPLES_DIR) not in sys.path:
    sys.path.insert(0, str(EXAMPLES_DIR))

from _bootstrap import bootstrap_example_paths

bootstrap_example_paths(__file__)

import osmnx as ox

from osm_chordify import build_osm_by_pop_density, create_osm_highway_filter
from osm_chordify.osm.graph import (
    format_validation_summary,
    summarize_graph_validation,
)


OSM_HIGHWAYS = [
    "motorway", "motorway_link", "trunk", "trunk_link",
    "primary", "primary_link", "secondary", "secondary_link",
    "tertiary", "tertiary_link", "unclassified",
]
OSM_RESIDENTIAL = ["residential"]
DEFAULT_OVERPASS_URLS = [
    "https://overpass-api.de/api",
    "https://overpass.kumi.systems/api",
    "https://overpass.private.coffee/api",
]

def base_osmnx_settings() -> dict:
    """Return the shared OSMnx settings used by the example builds."""
    return {
        "log_console": True,
        "use_cache": True,
        "cache_only_mode": False,
        "all_oneway": True,
        "requests_timeout": 180,
        "overpass_memory": None,
        "max_query_area_size": 50 * 1000 * 50 * 1000,
        "overpass_rate_limit": False,
        "overpass_max_attempts": 3,
        "useful_tags_way": list(ox.settings.useful_tags_way) + [
            "maxweight", "hgv", "maxweight:hgv", "maxlength",
            "motorcar", "motor_vehicle", "goods", "truck",
        ],
        "overpass_url": DEFAULT_OVERPASS_URLS[0],
        "overpass_urls": DEFAULT_OVERPASS_URLS,
    }


def highway_filter(include_residential: bool = False) -> str:
    """Return the standard highway custom filter for example downloads."""
    highway_types = set(OSM_HIGHWAYS)
    if include_residential:
        highway_types.update(OSM_RESIDENTIAL)
    return create_osm_highway_filter(sorted(highway_types))


def base_osm_config(*, tolerance: int, strongly_connected_components: bool) -> dict:
    """Return the shared top-level OSM config fields for example builds."""
    return {
        "osmnx_settings": base_osmnx_settings(),
        "weight_limits": {"unit": "lbs", "mdv_max": 26000, "hdv_max": 80000},
        "download_enabled": True,
        "tolerance": tolerance,
        "strongly_connected_components": strongly_connected_components,
    }


def parse_build_args(default_output_dir: str):
    """Parse common CLI arguments for example build scripts."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default=default_output_dir)
    parser.add_argument("--census-api-key", default=os.environ.get("CENSUS_API_KEY", ""))
    parser.add_argument(
        "--skip-validation",
        action="store_true",
        help="Skip post-build graph/XML integrity checks.",
    )
    return parser.parse_args()


def run_example_build(*, default_output_dir: str, osm_config: dict, area_config: dict, geo_config: dict):
    """Run a build example with shared CLI parsing and optional validation."""
    args = parse_build_args(default_output_dir)
    if args.census_api_key:
        os.environ["CENSUS_API_KEY"] = args.census_api_key

    result = build_osm_by_pop_density(
        work_dir=args.output_dir,
        osm_config=osm_config,
        area_config=area_config,
        geo_config=geo_config,
    )

    if not args.skip_validation:
        validate_built_network(
            result["graph"],
            osm_path=result["exported"].get("osm"),
        )

    return result


def validate_exported_osm_xml(osm_path: str) -> dict:
    """Validate that exported OSM XML has unique coordinates and valid node refs."""
    root = ET.parse(osm_path).getroot()
    coords = [(node.attrib["lat"], node.attrib["lon"]) for node in root.findall("node")]
    node_ids = {node.attrib["id"] for node in root.findall("node")}
    dangling = [
        ref
        for way in root.findall("way")
        for nd in way.findall("nd")
        if (ref := nd.attrib["ref"]) not in node_ids
    ]
    return {
        "xml_duplicate_coordinate_pairs": len(coords) - len(set(coords)),
        "xml_dangling_nd_refs": len(dangling),
        "xml_dangling_examples": dangling[:5],
    }

def validate_built_network(graph, osm_path: str | None = None, close_threshold_m: float = 0.5) -> dict:
    """Validate graph- and XML-level integrity for a built network."""
    summary = summarize_graph_validation(graph, close_threshold_m=close_threshold_m)

    if osm_path:
        summary.update(validate_exported_osm_xml(osm_path))

    print(format_validation_summary(summary, title="Validation checks"))

    problems = []
    if summary["self_loops"]:
        problems.append(f"self_loops={summary['self_loops']}")
    if summary["isolated_nodes"]:
        problems.append(f"isolated_nodes={summary['isolated_nodes']}")
    if not summary["weakly_connected"]:
        problems.append("graph is not weakly connected")
    if summary["duplicate_xml_coordinate_pairs"]:
        problems.append(
            f"duplicate_xml_coordinate_pairs={summary['duplicate_xml_coordinate_pairs']}"
        )
    if summary["close_node_pairs_lt_0_5m"]:
        problems.append(
            f"close_node_pairs_lt_0_5m={summary['close_node_pairs_lt_0_5m']}"
        )
    if summary["missing_edge_id"]:
        problems.append(f"missing_edge_id={summary['missing_edge_id']}")
    if summary["missing_length"]:
        problems.append(f"missing_length={summary['missing_length']}")
    if summary["nonpositive_length"]:
        problems.append(f"nonpositive_length={summary['nonpositive_length']}")
    if summary["missing_speed_kph"]:
        problems.append(f"missing_speed_kph={summary['missing_speed_kph']}")
    if summary["invalid_oneway_values"]:
        problems.append(f"invalid_oneway_values={summary['invalid_oneway_values']}")
    if summary["missing_geometry"]:
        problems.append(f"missing_geometry={summary['missing_geometry']}")
    if summary.get("xml_duplicate_coordinate_pairs", 0):
        problems.append(
            f"xml_duplicate_coordinate_pairs={summary['xml_duplicate_coordinate_pairs']}"
        )
    if summary.get("xml_dangling_nd_refs", 0):
        problems.append(f"xml_dangling_nd_refs={summary['xml_dangling_nd_refs']}")

    if problems:
        raise RuntimeError(
            "Built network validation failed: " + ", ".join(problems) +
            "\n" + json.dumps(summary, indent=2)
        )

    return summary
