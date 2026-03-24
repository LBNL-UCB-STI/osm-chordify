#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Orchestration pipeline for OSM network download, export, and polygon
intersection.

Combines network download/build with polygon-grid intersection and
tabular-network mapping into a single module.

@author: haitamlaarabi, cristian.poliziani, zaneedell
"""

import argparse
import json
import logging
import os
import pickle
import xml.etree.ElementTree as ET

import geopandas as gpd
import pandas as pd
from osm_chordify.utils.io import save_dataframe
from osm_chordify.utils.network import map_network_to_intersection

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# OSM download & build
# ---------------------------------------------------------------------------

def build_osm_by_pop_density(
    work_dir, osm_config, area_config, geo_config,
):
    """Download, process, validate, and export an OSM network.

    Parameters
    ----------
    work_dir : str
        Root working directory.
    osm_config : dict
        OSM configuration (osmnx_settings, graph_layers, etc.).
    area_config : dict
        Area definition (name, state_fips, county_fips, census_year, etc.).
    geo_config : dict
        Geo settings (utm_epsg, etc.).
    """
    from osm_chordify.osm.diagnostics import check_invalid_coordinates
    from osm_chordify.osm.export import export_network
    from osm_chordify.osm.graph import download_and_prepare_osm_network
    from osm_chordify.utils.geo import name_osm_network

    area_name = area_config["name"]
    osm_name = name_osm_network(
        area_name,
        osm_config["graph_layers"],
        osm_config.get("strongly_connected_components", False),
    )
    os.makedirs(work_dir, exist_ok=True)
    os.makedirs(os.path.join(work_dir, "geo"), exist_ok=True)
    os.makedirs(os.path.join(work_dir, "network"), exist_ok=True)
    osm_dir = f'{work_dir}/network/{osm_name}'

    logger.info("Downloading and preparing OSM-based %s network...", osm_name)
    g_osm = download_and_prepare_osm_network(
        osm_config, area_config, geo_config, work_dir,
    )

    has_invalid, invalid_nodes = check_invalid_coordinates(g_osm)
    if has_invalid:
        logger.warning("Found %d nodes with invalid coordinates.", len(invalid_nodes))
    else:
        logger.info("All node coordinates are valid.")

    exported = export_network(
        g_osm,
        output_dir=osm_dir,
        name=osm_name,
        edge_tag_aggs=[('length', 'sum')],
    )
    return {
        "graph": g_osm,
        "name": osm_name,
        "output_dir": osm_dir,
        "exported": exported,
    }


# ---------------------------------------------------------------------------
# Network mapping
# ---------------------------------------------------------------------------

def map_osm_with_beam_network(
    osm_path, network_path, network_osm_id_col="attributeOrigId",
    output_path=None,
):
    """Map OSM edges to a BEAM network on a shared OSM ID.

    Joins the BEAM network table to OSM edge geometries.  Accepts
    ``.gpkg``, ``.geojson``, or ``.pbf`` as the OSM source.  For
    ``.gpkg`` the ``edges`` layer is read directly; for ``.pbf`` /
    ``.geojson`` the ``other_tags`` field is parsed to extract
    ``edge_id`` and ``edge_length``.

    All columns from both inputs are included in the output.

    Parameters
    ----------
    osm_path : str
        Path to the OSM network file (``.gpkg``, ``.geojson``, or
        ``.osm.pbf``).
    network_path : str
        Path to the BEAM network CSV / CSV.GZ file.
    network_osm_id_col : str
        Column in the BEAM network that holds the OSM ID for joining.
        Default ``"attributeOrigId"``.
    output_path : str, optional
        If provided, save the result to this file (GeoJSON, GPKG, or CSV).

    Returns
    -------
    pd.DataFrame
        Merged result with all columns from both inputs.
    """
    logger.info("Loading network from %s", network_path)
    network_df = pd.read_csv(network_path)

    logger.info("Loading OSM data from %s", osm_path)
    if osm_path.endswith(".gpkg"):
        from osm_chordify.osm.intersect import load_osm_edges
        osm_gdf = load_osm_edges(osm_path)
    elif osm_path.lower().endswith(".pbf"):
        osm_gdf = _read_osm_pbf_lines(osm_path)
        osm_gdf = _parse_osm_tags(osm_gdf)
    else:
        osm_gdf = gpd.read_file(osm_path)
        if "other_tags" in osm_gdf.columns:
            osm_gdf = _parse_osm_tags(osm_gdf)

    osm_id_col = "edge_osm_id" if "edge_osm_id" in osm_gdf.columns else "osm_id"

    merged_df = map_network_to_intersection(
        network_df=network_df,
        intersection_gdf=osm_gdf,
        network_osm_id_col=network_osm_id_col,
        intersection_osm_id_col=osm_id_col,
    )

    # Verification summary
    n_network_ids = network_df[network_osm_id_col].dropna().nunique()
    n_osm_ids = osm_gdf[osm_id_col].nunique()
    n_matched = merged_df[network_osm_id_col].nunique() if len(merged_df) else 0
    match_rate = n_matched / n_network_ids if n_network_ids else 0
    logger.info(
        "Match summary — network IDs: %d, OSM IDs: %d, matched: %d (%.1f%%)",
        n_network_ids, n_osm_ids, n_matched, match_rate * 100,
    )

    if output_path:
        save_dataframe(merged_df, output_path)
    logger.info("Network mapping complete")
    return merged_df


def map_network_csv_to_osm_pbf(
    network_csv_path,
    osm_pbf_path,
    output_path,
    network_osm_id_col="attributeOrigId",
):
    """Map a network CSV/CSV.GZ to an OSM PBF and save a spatial join."""
    if not str(osm_pbf_path).lower().endswith(".pbf"):
        raise ValueError(
            f"osm_pbf_path must point to a .pbf/.osm.pbf file, got: {osm_pbf_path}"
        )
    if os.path.splitext(str(output_path))[1].lower() not in {".parquet", ".gpkg", ".geojson"}:
        raise ValueError("output_path must end with .parquet, .gpkg, or .geojson")

    return map_osm_with_beam_network(
        osm_path=osm_pbf_path,
        network_path=network_csv_path,
        network_osm_id_col=network_osm_id_col,
        output_path=output_path,
    )


def _read_osm_pbf_lines(osm_path):
    """Read the lines layer from an OSM PBF with robust GDAL indexing settings."""
    previous = os.environ.get("OSM_USE_CUSTOM_INDEXING")
    os.environ["OSM_USE_CUSTOM_INDEXING"] = "NO"
    try:
        return gpd.read_file(osm_path, layer="lines")
    finally:
        if previous is None:
            os.environ.pop("OSM_USE_CUSTOM_INDEXING", None)
        else:
            os.environ["OSM_USE_CUSTOM_INDEXING"] = previous


def _parse_osm_tags(osm_gdf):
    """Extract ``edge_id`` and ``edge_length`` from ``other_tags``."""
    from osm_chordify.osm.tags import parse_other_tags, extract_tag_as_float

    osm_gdf = osm_gdf.copy()
    osm_gdf["osm_id"] = osm_gdf["osm_id"].astype(int)
    parsed = osm_gdf["other_tags"].apply(parse_other_tags)
    osm_gdf["edge_id"] = parsed.apply(lambda t: t.get("edge_id"))
    osm_gdf["edge_length"] = parsed.apply(
        lambda t: extract_tag_as_float(t, "length")
    )
    return osm_gdf


def match_road_network_geometries(
    network_a,
    network_a_epsg,
    network_b,
    network_b_epsg,
    matching="flexible",
    output_path=None,
    output_epsg=None,
):
    """Match link geometries between two road networks.

    Spatially matches edges from *network_a* to *network_b* based on
    geometric proximity and overlap.  All attributes from both networks
    are included in the result.

    Both *network_a* and *network_b* accept either a file path (GPKG,
    GeoJSON, Shapefile, …) or a :class:`~geopandas.GeoDataFrame`.
    When a ``.gpkg`` path is given, the ``edges`` layer is read
    automatically.

    Parameters
    ----------
    network_a : str, os.PathLike, or gpd.GeoDataFrame
        First road network (edges with line geometries).
    network_a_epsg : int
        EPSG code for *network_a*'s CRS (e.g. ``4326`` for WGS 84).
    network_b : str, os.PathLike, or gpd.GeoDataFrame
        Second road network to match against.
    network_b_epsg : int
        EPSG code for *network_b*'s CRS.
    output_epsg : int, optional
        EPSG code for the output CRS.  Defaults to *network_a_epsg*
        if not provided.
    matching : str, optional
        Matching strategy — ``"strict"`` requires high geometric overlap
        between edges, ``"flexible"`` (default) allows partial and
        nearby matches.
    output_path : str, optional
        If provided, save the result to this file.

    Returns
    -------
    gpd.GeoDataFrame
        Matched edge pairs with all attributes from both networks
        (prefixed ``a_`` and ``b_``) and match-quality metrics.
    """
    raise NotImplementedError("match_road_network_geometries is not yet implemented")


# ---------------------------------------------------------------------------
# OSM diagnostics
# ---------------------------------------------------------------------------

def collect_osm_diagnostics(pbf_path, epsg_utm):
    """Collect connectivity and length diagnostics for an OSM ``.pbf`` file."""
    import warnings

    import networkx as nx
    import pandas as pd

    try:
        from pyrosm import OSM
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            "PBF diagnostics require the optional 'pyrosm' dependency. "
            "Install it with `pip install \"osm-chordify[diagnostics]\"` "
            "or `pip install pyrosm`."
        ) from exc

    warnings.filterwarnings('ignore', category=FutureWarning)
    chained_assignment_error = getattr(pd.errors, "ChainedAssignmentError", None)
    if chained_assignment_error is not None:
        warnings.filterwarnings('ignore', category=chained_assignment_error)
    setting_with_copy_warning = getattr(pd.errors, "SettingWithCopyWarning", None)
    if setting_with_copy_warning is not None:
        warnings.filterwarnings('ignore', category=setting_with_copy_warning)

    pbf_path = os.path.expanduser(pbf_path)
    osm = OSM(pbf_path)

    try:
        nodes, edges = osm.get_network(network_type="driving", nodes=True)
        retrieved_nodes = len(nodes)
        retrieved_edges = len(edges)
    except Exception:
        edges = osm.get_network(network_type="driving")
        retrieved_nodes = None
        retrieved_edges = len(edges)

    edges_projected = edges.to_crs(f'EPSG:{epsg_utm}')
    edges_projected['length_m'] = edges_projected.geometry.length

    G = nx.Graph()
    node_coords_to_id = {}
    next_node_id = 0

    for _, row in edges_projected.iterrows():
        coords = list(row.geometry.coords)
        start_coord = coords[0]
        end_coord = coords[-1]

        if start_coord not in node_coords_to_id:
            node_coords_to_id[start_coord] = next_node_id
            next_node_id += 1

        if end_coord not in node_coords_to_id:
            node_coords_to_id[end_coord] = next_node_id
            next_node_id += 1

    components = list(nx.connected_components(G))
    components_sorted = sorted(components, key=len, reverse=True)
    largest = components_sorted[0]
    largest_subgraph = G.subgraph(largest)

    short_links = edges_projected[edges_projected['length_m'] < 15].copy()
    short_links = short_links.sort_values('length_m')

    long_links = edges_projected[edges_projected['length_m'] > 10000].copy()
    long_links = long_links.sort_values('length_m', ascending=False)
    edges_no_outliers = edges_projected[edges_projected['length_m'] <= 500].copy()

    return {
        "pbf_path": pbf_path,
        "retrieved_nodes": retrieved_nodes,
        "retrieved_edges": retrieved_edges,
        "connected_components": len(components),
        "graph_nodes": G.number_of_nodes(),
        "graph_edges": G.number_of_edges(),
        "largest_component_nodes": len(largest),
        "largest_component_edges": largest_subgraph.number_of_edges(),
        "top_component_sizes": [len(c) for c in components_sorted[: min(10, len(components_sorted))]],
        "short_links_lt_15m": len(short_links),
        "long_links_gt_10km": len(long_links),
        "short_links_min_m": float(short_links['length_m'].min()) if len(short_links) else None,
        "short_links_max_m": float(short_links['length_m'].max()) if len(short_links) else None,
        "short_links_mean_m": float(short_links['length_m'].mean()) if len(short_links) else None,
        "short_links_median_m": float(short_links['length_m'].median()) if len(short_links) else None,
        "links_le_500m_count": len(edges_no_outliers),
        "links_le_500m_min_m": float(edges_no_outliers['length_m'].min()) if len(edges_no_outliers) else None,
        "links_le_500m_max_m": float(edges_no_outliers['length_m'].max()) if len(edges_no_outliers) else None,
        "links_le_500m_mean_m": float(edges_no_outliers['length_m'].mean()) if len(edges_no_outliers) else None,
        "links_le_500m_median_m": float(edges_no_outliers['length_m'].median()) if len(edges_no_outliers) else None,
        "edges_no_outliers": edges_no_outliers,
    }


def format_osm_diagnostics_summary(summary: dict, title: str = "OSM PBF diagnostics") -> str:
    """Render a readable diagnostics report for an OSM PBF artifact."""
    return "\n".join(
        [
            f"{title}:",
            f"- retrieved nodes: {summary['retrieved_nodes']}",
            f"- retrieved edges: {summary['retrieved_edges']}",
            f"- connected components: {summary['connected_components']}",
            f"- graph nodes: {summary['graph_nodes']}",
            f"- graph edges: {summary['graph_edges']}",
            f"- largest component nodes: {summary['largest_component_nodes']}",
            f"- largest component edges: {summary['largest_component_edges']}",
            f"- top component sizes: {summary['top_component_sizes']}",
            f"- short links <15m: {summary['short_links_lt_15m']}",
            f"- long links >10km: {summary['long_links_gt_10km']}",
            f"- short link min/max/mean/median (m): {summary['short_links_min_m']}/{summary['short_links_max_m']}/{summary['short_links_mean_m']}/{summary['short_links_median_m']}",
            f"- links <=500m count: {summary['links_le_500m_count']}",
            f"- links <=500m min/max/mean/median (m): {summary['links_le_500m_min_m']}/{summary['links_le_500m_max_m']}/{summary['links_le_500m_mean_m']}/{summary['links_le_500m_median_m']}",
        ]
    )


def infer_sidecar_paths(pbf_path):
    """Infer sibling graph and OSM XML artifacts for a built PBF."""
    pbf_path = os.path.abspath(os.path.expanduser(pbf_path))
    base_name = os.path.basename(pbf_path)
    parent = os.path.dirname(pbf_path)
    if base_name.endswith(".osm.pbf"):
        stem = base_name[:-8]
        osm_path = os.path.join(parent, f"{stem}.osm")
        pkl_path = os.path.join(parent, f"{stem}.pkl")
        graphml_path = os.path.join(parent, f"{stem}.graphml")
    else:
        stem, _ = os.path.splitext(base_name)
        osm_path = os.path.join(parent, f"{stem}.osm")
        pkl_path = os.path.join(parent, f"{stem}.pkl")
        graphml_path = os.path.join(parent, f"{stem}.graphml")

    graph_path = pkl_path if os.path.exists(pkl_path) else graphml_path if os.path.exists(graphml_path) else None
    return graph_path, osm_path if os.path.exists(osm_path) else None


def load_built_graph(graph_path):
    """Load a built graph from a pickle or GraphML artifact."""
    import osmnx as ox

    graph_path = os.path.abspath(os.path.expanduser(graph_path))
    if graph_path.endswith(".pkl"):
        with open(graph_path, "rb") as f:
            return pickle.load(f)
    if graph_path.endswith(".graphml"):
        return ox.load_graphml(graph_path)
    raise ValueError(f"Unsupported graph sidecar format: {graph_path}")


def build_validation_summary(graph, osm_path=None, close_threshold_m: float = 0.5):
    """Build graph/XML validation metrics without raising."""
    from osm_chordify.osm.graph import summarize_graph_validation

    summary = summarize_graph_validation(graph, close_threshold_m=close_threshold_m)
    if osm_path:
        root = ET.parse(osm_path).getroot()
        coords = [(node.attrib["lat"], node.attrib["lon"]) for node in root.findall("node")]
        node_ids = {node.attrib["id"] for node in root.findall("node")}
        dangling = [
            ref
            for way in root.findall("way")
            for nd in way.findall("nd")
            if (ref := nd.attrib["ref"]) not in node_ids
        ]
        summary.update(
            {
                "xml_duplicate_coordinate_pairs": len(coords) - len(set(coords)),
                "xml_dangling_nd_refs": len(dangling),
                "xml_dangling_examples": dangling[:5],
            }
        )
    return summary


def validate_built_network(graph, osm_path=None, close_threshold_m: float = 0.5):
    """Validate graph- and XML-level integrity for a built network."""
    from osm_chordify.osm.graph import format_validation_summary

    summary = build_validation_summary(graph, osm_path=osm_path, close_threshold_m=close_threshold_m)
    print(format_validation_summary(summary, title="Validation checks"))

    problems = []
    if summary["unprotected_self_loops"]:
        problems.append(f"unprotected_self_loops={summary['unprotected_self_loops']}")
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


def diagnose_built_osm_pbf(
    pbf_path,
    epsg_utm,
    graph_path=None,
    osm_xml_path=None,
    skip_pbf_diagnostics=False,
):
    """Run built-graph validation plus PBF diagnostics for one artifact."""
    pbf_path = os.path.abspath(os.path.expanduser(pbf_path))
    graph_path = os.path.abspath(os.path.expanduser(graph_path)) if graph_path else None
    osm_xml_path = os.path.abspath(os.path.expanduser(osm_xml_path)) if osm_xml_path else None

    if graph_path is None or osm_xml_path is None:
        default_graph_path, default_osm_path = infer_sidecar_paths(pbf_path)
        graph_path = graph_path or default_graph_path
        osm_xml_path = osm_xml_path or default_osm_path

    if graph_path and os.path.exists(graph_path):
        print("=== Built Graph Validation ===")
        graph = load_built_graph(graph_path)
        validate_built_network(
            graph,
            osm_path=osm_xml_path if osm_xml_path and os.path.exists(osm_xml_path) else None,
        )
        print("=== End Built Graph Validation ===")
    else:
        print("No sibling built graph (.pkl or .graphml) found; skipping build-graph validation.")

    if not skip_pbf_diagnostics:
        print("=== OSM PBF Diagnostics ===")
        diagnose_osm(pbf_path, epsg_utm)
        print("=== End OSM PBF Diagnostics ===")


def compare_osm_pbf_artifacts(
    pbf_a,
    pbf_b,
    epsg_utm_a,
    epsg_utm_b=None,
    graph_a=None,
    graph_b=None,
    osm_xml_a=None,
    osm_xml_b=None,
):
    """Compare validation and diagnostics metrics across two built OSM PBF artifacts."""
    from osm_chordify.osm.graph import format_validation_summary

    def _resolve_sidecars(pbf_path, graph_path, osm_xml_path):
        graph = os.path.abspath(os.path.expanduser(graph_path)) if graph_path else None
        osm_xml = os.path.abspath(os.path.expanduser(osm_xml_path)) if osm_xml_path else None
        if graph is None or osm_xml is None:
            default_graph, default_osm = infer_sidecar_paths(pbf_path)
            graph = graph or default_graph
            osm_xml = osm_xml or default_osm
        return graph, osm_xml

    def _print_metric_deltas(label, summary_a, summary_b, keys):
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

    pbf_a = os.path.abspath(os.path.expanduser(pbf_a))
    pbf_b = os.path.abspath(os.path.expanduser(pbf_b))
    epsg_utm_b = epsg_utm_b or epsg_utm_a

    graph_a, osm_xml_a = _resolve_sidecars(pbf_a, graph_a, osm_xml_a)
    graph_b, osm_xml_b = _resolve_sidecars(pbf_b, graph_b, osm_xml_b)

    build_summary_a = None
    build_summary_b = None
    if graph_a and os.path.exists(graph_a):
        print("=== Artifact A Built Graph Validation ===")
        build_summary_a = build_validation_summary(
            load_built_graph(graph_a),
            osm_path=osm_xml_a if osm_xml_a and os.path.exists(osm_xml_a) else None,
        )
        print(format_validation_summary(build_summary_a, title="Artifact A validation checks"))
    else:
        print("Artifact A: no sibling built graph (.pkl or .graphml) found; skipping build-graph validation.")

    if graph_b and os.path.exists(graph_b):
        print("=== Artifact B Built Graph Validation ===")
        build_summary_b = build_validation_summary(
            load_built_graph(graph_b),
            osm_path=osm_xml_b if osm_xml_b and os.path.exists(osm_xml_b) else None,
        )
        print(format_validation_summary(build_summary_b, title="Artifact B validation checks"))
    else:
        print("Artifact B: no sibling built graph (.pkl or .graphml) found; skipping build-graph validation.")

    print("=== Artifact A OSM PBF Diagnostics ===")
    pbf_summary_a = collect_osm_diagnostics(pbf_a, epsg_utm_a)
    print(format_osm_diagnostics_summary(pbf_summary_a, title="Artifact A OSM PBF diagnostics"))

    print("=== Artifact B OSM PBF Diagnostics ===")
    pbf_summary_b = collect_osm_diagnostics(pbf_b, epsg_utm_b)
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


def diagnose_osm(pbf_path, epsg_utm):
    """Load an OSM ``.pbf`` file and run connectivity / length diagnostics.

    Checks connected components, reports short (< 15 m) and long (> 10 km)
    link anomalies, prints summary statistics, and saves a link-length
    histogram to ``all_links_histogram.png``.

    Parameters
    ----------
    pbf_path : str
        Path to the ``.osm.pbf`` file to analyse.
    epsg_utm : int
        EPSG code for the UTM projection used for length calculation.
    """

    try:
        import matplotlib.pyplot as plt
    except ModuleNotFoundError:
        plt = None

    summary = collect_osm_diagnostics(pbf_path, epsg_utm)
    logger.info(format_osm_diagnostics_summary(summary, title="OSM PBF diagnostics"))

    if plt is None:
        logger.warning(
            "matplotlib is not installed; skipping link-length histogram generation."
        )
    elif len(summary["edges_no_outliers"]) > 0:
        plt.figure(figsize=(12, 6))
        plt.hist(summary["edges_no_outliers"]['length_m'], bins=50, edgecolor='black', alpha=0.7)
        plt.xlabel('Length (meters)')
        plt.ylabel('Number of Links')
        plt.title('Distribution of All Link Lengths in Network')
        plt.grid(True, alpha=0.3)

        mean_len = summary["edges_no_outliers"]['length_m'].mean()
        median_len = summary["edges_no_outliers"]['length_m'].median()

        plt.axvline(mean_len, color='red', linestyle='--',
                    linewidth=2, label=f'Mean: {mean_len:.2f}m')
        plt.axvline(median_len, color='green', linestyle='--',
                    linewidth=2, label=f'Median: {median_len:.2f}m')
        plt.legend()

        plt.tight_layout()
        plt.savefig('all_links_histogram.png', dpi=300)
        logger.info("Histogram saved to 'all_links_histogram.png'")
        plt.show()
    else:
        logger.info("Skipping histogram: no links after outlier filtering.")


def _load_json_config(config_path):
    """Load a JSON config file for CLI use."""
    with open(config_path, "r", encoding="utf-8") as f:
        return json.load(f)


def _build_arg_parser():
    """Build the argparse CLI for ``python -m osm_chordify.main``."""
    parser = argparse.ArgumentParser(
        description="CLI entrypoints for osm-chordify build, mapping, intersection, and diagnostics workflows."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    build_parser = subparsers.add_parser(
        "build",
        help="Build and export an OSM network from JSON config files.",
    )
    build_parser.add_argument("--work-dir", required=True)
    build_parser.add_argument("--osm-config", required=True)
    build_parser.add_argument("--area-config", required=True)
    build_parser.add_argument("--geo-config", required=True)

    intersect_parser = subparsers.add_parser(
        "intersect",
        help="Intersect a road network with polygon zones.",
    )
    intersect_parser.add_argument("--road-network", required=True)
    intersect_parser.add_argument("--road-network-epsg", required=True, type=int)
    intersect_parser.add_argument("--zones", required=True)
    intersect_parser.add_argument("--output-path", required=True)
    intersect_parser.add_argument("--output-epsg", type=int)

    map_parser = subparsers.add_parser(
        "map",
        help="Map a BEAM/network CSV to OSM geometries and optionally save the result.",
    )
    map_parser.add_argument("--osm-path", required=True)
    map_parser.add_argument("--network-path", required=True)
    map_parser.add_argument("--output-path")
    map_parser.add_argument("--network-osm-id-col", default="attributeOrigId")

    diagnose_parser = subparsers.add_parser(
        "diagnose",
        help="Run diagnostics on a built OSM PBF file.",
    )
    diagnose_parser.add_argument("--pbf-path", required=True)
    diagnose_parser.add_argument("--epsg-utm", required=True, type=int)

    diagnose_built_parser = subparsers.add_parser(
        "diagnose-built",
        help="Run built-graph validation plus PBF diagnostics for one built artifact.",
    )
    diagnose_built_parser.add_argument("--pbf-path", required=True)
    diagnose_built_parser.add_argument("--epsg-utm", required=True, type=int)
    diagnose_built_parser.add_argument("--graph-path")
    diagnose_built_parser.add_argument("--osm-xml")
    diagnose_built_parser.add_argument("--skip-pbf-diagnostics", action="store_true")

    compare_parser = subparsers.add_parser(
        "compare-pbf",
        help="Compare validation and diagnostics metrics across two built OSM PBF artifacts.",
    )
    compare_parser.add_argument("--pbf-a", required=True)
    compare_parser.add_argument("--pbf-b", required=True)
    compare_parser.add_argument("--epsg-utm-a", required=True, type=int)
    compare_parser.add_argument("--epsg-utm-b", type=int)
    compare_parser.add_argument("--graph-a")
    compare_parser.add_argument("--graph-b")
    compare_parser.add_argument("--osm-xml-a")
    compare_parser.add_argument("--osm-xml-b")

    map_pbf_parser = subparsers.add_parser(
        "map-pbf",
        help="Map a network CSV/CSV.GZ to an OSM PBF and save a spatial join.",
    )
    map_pbf_parser.add_argument("--network-csv-path", required=True)
    map_pbf_parser.add_argument("--osm-pbf-path", required=True)
    map_pbf_parser.add_argument("--output-path", required=True)
    map_pbf_parser.add_argument("--network-osm-id-col", default="attributeOrigId")

    return parser


def _run_cli(args):
    """Dispatch parsed CLI args to the existing workflow functions."""
    if args.command == "build":
        result = build_osm_by_pop_density(
            work_dir=args.work_dir,
            osm_config=_load_json_config(args.osm_config),
            area_config=_load_json_config(args.area_config),
            geo_config=_load_json_config(args.geo_config),
        )
        logger.info("Build complete: %s", result["output_dir"])
        return result

    if args.command == "intersect":
        from osm_chordify.osm.intersect import intersect_road_network_with_zones

        return intersect_road_network_with_zones(
            road_network=args.road_network,
            road_network_epsg=args.road_network_epsg,
            zones=args.zones,
            output_path=args.output_path,
            output_epsg=args.output_epsg,
        )

    if args.command == "map":
        return map_osm_with_beam_network(
            osm_path=args.osm_path,
            network_path=args.network_path,
            network_osm_id_col=args.network_osm_id_col,
            output_path=args.output_path,
        )

    if args.command == "diagnose":
        return diagnose_osm(args.pbf_path, args.epsg_utm)

    if args.command == "diagnose-built":
        return diagnose_built_osm_pbf(
            pbf_path=args.pbf_path,
            epsg_utm=args.epsg_utm,
            graph_path=args.graph_path,
            osm_xml_path=args.osm_xml,
            skip_pbf_diagnostics=args.skip_pbf_diagnostics,
        )

    if args.command == "compare-pbf":
        return compare_osm_pbf_artifacts(
            pbf_a=args.pbf_a,
            pbf_b=args.pbf_b,
            epsg_utm_a=args.epsg_utm_a,
            epsg_utm_b=args.epsg_utm_b,
            graph_a=args.graph_a,
            graph_b=args.graph_b,
            osm_xml_a=args.osm_xml_a,
            osm_xml_b=args.osm_xml_b,
        )

    if args.command == "map-pbf":
        return map_network_csv_to_osm_pbf(
            network_csv_path=args.network_csv_path,
            osm_pbf_path=args.osm_pbf_path,
            output_path=args.output_path,
            network_osm_id_col=args.network_osm_id_col,
        )

    raise ValueError(f"Unsupported command: {args.command}")


def main(argv=None):
    """CLI entrypoint for ``python -m osm_chordify.main``."""
    parser = _build_arg_parser()
    args = parser.parse_args(argv)
    return _run_cli(args)


if __name__ == "__main__":
    main()


# ---------------------------------------------------------------------------
