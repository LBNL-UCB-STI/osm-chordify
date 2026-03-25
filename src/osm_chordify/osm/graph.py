"""Graph processing: ferry edges, topology, validation, and download pipeline."""

import hashlib
import json
import logging
import os
import pickle
from statistics import mean

import networkx as nx
import osmnx as ox
import pandas as pd
import shapely.ops
from networkx.algorithms import strongly_connected_components, weakly_connected_components

from osm_chordify.osm.simplify import (
    bool_all,
    bool_any,
    first_valid_value,
    mean_maxspeed,
    median_lanes,
    min_numeric_or_string,
    yes_no_all,
)
from osm_chordify.osm.tags import (
    standardize_access,
    standardize_hgv,
    standardize_maxspeed,
    standardize_motor_vehicle,
    standardize_oneway,
    standardize_weight,
)
from osm_chordify.utils.data_collection import (
    collect_census_data,
    collect_geographic_boundaries,
    filter_boundaries_by_density,
)
from osm_chordify.utils.geo import build_area_mask_geometry, project_graph, to_convex_hull

logger = logging.getLogger(__name__)

PROTECTED_HIGHWAY_TYPES = {
    "motorway",
    "motorway_link",
    "trunk",
    "trunk_link",
    "primary",
    "primary_link",
}


def _get_overpass_urls(network_config: dict) -> list[str]:
    """Return the ordered list of Overpass endpoints to try."""
    urls = network_config["osmnx_settings"].get("overpass_urls")
    if urls:
        return list(dict.fromkeys(urls))

    url = network_config["osmnx_settings"].get("overpass_url")
    return [url] if url else []


def _download_with_overpass_fallback(download_fn, network_config: dict, context: str):
    """Try multiple Overpass endpoints until one succeeds."""
    overpass_urls = _get_overpass_urls(network_config)
    if not overpass_urls:
        return download_fn()

    original_url = ox.settings.overpass_url
    errors = []
    try:
        for overpass_url in overpass_urls:
            ox.settings.overpass_url = overpass_url
            try:
                logger.info("Downloading %s via Overpass endpoint %s", context, overpass_url)
                return download_fn()
            except Exception as exc:  # pragma: no cover - exercised via monkeypatch in tests
                logger.warning(
                    "Overpass download failed for %s via %s: %s",
                    context,
                    overpass_url,
                    exc,
                )
                errors.append(f"{overpass_url}: {exc}")
    finally:
        ox.settings.overpass_url = original_url

    raise RuntimeError(
        f"All Overpass endpoints failed for {context}. "
        f"Tried: {overpass_urls}. Errors: {' | '.join(errors)}"
    )


def _raw_graph_cache_path(network_config: dict, area_config: dict, work_dir: str) -> str:
    """Build a cache path that changes when the area/layer download config changes."""
    cache_descriptor = {
        "study_area": area_config["name"],
        "census_year": area_config["census_year"],
        "state_fips": area_config["state_fips"],
        "county_fips": area_config["county_fips"],
        "graph_layers": network_config["graph_layers"],
    }
    cache_hash = hashlib.md5(
        json.dumps(cache_descriptor, sort_keys=True).encode("utf-8")
    ).hexdigest()[:10]
    return os.path.join(
        work_dir,
        "network",
        f"raw_osm_graph_{area_config['name']}_{cache_hash}.pkl",
    )


def process_ferry_edges(ferry_graph) -> nx.MultiDiGraph:
    """Process ferry edges to make them compatible with car network"""
    if ferry_graph.number_of_edges() == 0:
        logger.info("No ferry edges found in the graph.")
        return nx.MultiDiGraph()

    # Extract nodes and edges
    ferry_nodes, ferry_edges = ox.graph_to_gdfs(ferry_graph)
    logger.info("Total ferry edges: %d", len(ferry_edges))

    # Print available columns to debug
    logger.debug("Available columns: %s", ferry_edges.columns.tolist())

    # Create default masks - assume access is allowed unless explicitly denied
    # This is more lenient and works better with OSM data which often lacks explicit tags
    car_mask = pd.Series(True, index=ferry_edges.index)

    # Check for explicit denials first
    if 'motorcar' in ferry_edges.columns:
        car_mask &= ~(ferry_edges['motorcar'] == 'no')
        logger.info("After motorcar check: %d car-accessible edges", car_mask.sum())

    if 'motor_vehicle' in ferry_edges.columns:
        motor_vehicle_denied = ferry_edges['motor_vehicle'] == 'no'
        car_mask &= ~motor_vehicle_denied
        logger.info("After motor_vehicle check: %d car-accessible edges", car_mask.sum())

    # Select ferry edges that allow passenger cars
    selected_edges = ferry_edges[car_mask].copy()

    if selected_edges.empty:
        logger.info("No ferry routes found that allow passenger cars")
        return nx.MultiDiGraph()

    logger.info("Found %d suitable ferry edges", len(selected_edges))

    # Set ferry attributes
    selected_edges['reversed'] = False
    selected_edges['maxspeed'] = "10 mph"
    selected_edges['highway'] = "unclassified"
    selected_edges['oneway'] = "no"
    selected_edges['lanes'] = "2"
    selected_edges["hgv"] = False  # Mark as not accessible to heavy-duty
    selected_edges["mdv"] = True  # Mark as accessible to medium-duty

    # Keep only nodes that are used by the filtered edges
    used_nodes = set(selected_edges.index.get_level_values(0)).union(
        set(selected_edges.index.get_level_values(1))
    )
    selected_nodes = ferry_nodes.loc[list(used_nodes)]

    # Reconstruct graph and project
    g_ferry_reconstructed = ox.graph_from_gdfs(selected_nodes, selected_edges)

    return g_ferry_reconstructed


def _normalize_highway_values(value) -> set[str]:
    """Normalize an edge highway value into a comparable set of strings."""
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return set()
    if isinstance(value, list):
        values = value
    else:
        values = [value]
    normalized = set()
    for item in values:
        if item is None or pd.isna(item):
            continue
        normalized.add(str(item))
    return normalized


def _is_truthy_osm_tag(value) -> bool:
    """Interpret common OSM-style truthy values."""
    if isinstance(value, bool):
        return value
    if value is None or pd.isna(value):
        return False
    return str(value).strip().lower() in {"yes", "true", "1"}


def _edge_is_protected(edge_data: dict) -> bool:
    """Return True when an edge should be preserved during conservative cleanup."""
    if edge_data.get("protected_backbone"):
        return True
    if edge_data.get("layer_role") == "backbone":
        return True
    highway_values = _normalize_highway_values(edge_data.get("highway"))
    if highway_values & PROTECTED_HIGHWAY_TYPES:
        return True
    if _is_truthy_osm_tag(edge_data.get("bridge")) or _is_truthy_osm_tag(edge_data.get("tunnel")):
        return True
    return False


def _split_self_loops_for_cleanup(G: nx.MultiDiGraph):
    """Partition self-loop edges into removable artifacts and retained protected loops."""
    removable = []
    retained = []
    for u, v, key in nx.selfloop_edges(G, keys=True):
        edge_data = G[u][v][key]
        should_remove = not _edge_is_protected(edge_data)
        if should_remove:
            removable.append((u, v, key))
        else:
            retained.append((u, v, key))
    return removable, retained


def process_tags(_g: nx.MultiDiGraph, config: dict) -> nx.MultiDiGraph:
    """Process vehicle classifications based on FHWA weight classes."""
    logger.info("Processing vehicle classifications...")

    # Get weight limits and unit from config
    weight_config = config["weight_limits"]
    target_unit = weight_config["unit"]
    mdv_max = weight_config["mdv_max"]
    hdv_max = weight_config["hdv_max"]

    # Get graph data while preserving MultiIndex
    nodes, edges = ox.graph_to_gdfs(_g)

    # Standardize tags.
    # Type design note:
    #   oneway / motor_vehicle / access  → string ("yes" / "-1" / "no") — OSM tag
    #                                       values kept as strings for XML export.
    #   hgv / mdv                        → bool — used only for internal boolean
    #                                       logic (masks, bool_all aggregation).
    #                                       Converted explicitly with astype(bool)
    #                                       at the end of this function.
    # oneway is fundamental to routing direction — osmnx always provides it.
    # A missing oneway column means the graph was not built correctly.
    if 'oneway' not in edges.columns:
        raise ValueError(
            "process_tags: graph edges are missing the required 'oneway' column. "
            "Ensure the graph was downloaded with the standard osmnx tag set."
        )

    # motor_vehicle, maxspeed, and access are optional OSM tags — many edges
    # legitimately lack them.  Fill with None so standardize_* functions return
    # their safe defaults ("yes", None, "yes" respectively).
    for _opt_col in ('motor_vehicle', 'maxspeed', 'access'):
        if _opt_col not in edges.columns:
            logger.warning(
                "process_tags: optional column '%s' not found in graph edges — "
                "defaulting all edges to None for this tag.",
                _opt_col,
            )
            edges[_opt_col] = None

    edges['oneway'] = edges['oneway'].apply(standardize_oneway)
    edges['motor_vehicle'] = edges['motor_vehicle'].apply(standardize_motor_vehicle)
    edges['maxspeed'] = edges['maxspeed'].apply(standardize_maxspeed)
    edges['access'] = edges['access'].apply(standardize_access)
    # Initialize hgv and mdv as True by default if they don't exist
    edges["mdv"] = True
    if "hgv" not in edges.columns:
        edges["hgv"] = True
    edges["hgv"] = edges["hgv"].apply(standardize_hgv)

    # Copy HGV weight restrictions if present
    if "maxweight:hgv" in edges.columns:
        hgv_mask = ~edges["maxweight:hgv"].isna()
        if hgv_mask.any():
            edges.loc[hgv_mask, "maxweight"] = edges.loc[hgv_mask, "maxweight:hgv"].copy()

    if "maxweight" in edges.columns:
        logger.info("Processing weight restrictions...")
        # Convert weights to standard unit specified in config
        edges["maxweight"] = edges["maxweight"].apply(lambda x: standardize_weight(x, target_unit))

        # Update hgv and mdv based on weight restrictions
        # Medium-duty vehicles are restricted when weight is below MDV limit
        mdv_restricted_mask = edges["maxweight"].notna() & (edges["maxweight"] <= mdv_max)
        edges.loc[mdv_restricted_mask, "mdv"] = False

        # Heavy-duty vehicles are restricted when weight is below HDV limit
        # Create a mask for MDVs being restricted
        mdv_is_restricted = edges["mdv"] == False
        # Combine masks properly
        hdv_restricted_mask = mdv_is_restricted | (edges["maxweight"].notna() & (edges["maxweight"] <= hdv_max))
        edges.loc[hdv_restricted_mask, "hgv"] = False

    # Process other restrictions like maxlength
    if "maxlength" in edges.columns:
        # If maxlength is set, assume heavy vehicles are restricted
        length_restricted_mask = ~edges["maxlength"].isna()
        edges.loc[length_restricted_mask, "hgv"] = False

    # Ensure hgv, mdv and oneway are strictly boolean
    edges["hgv"] = edges["hgv"].astype(bool)
    edges["mdv"] = edges["mdv"].astype(bool)

    # Convert back to MultiDiGraph
    g_updated = ox.graph_from_gdfs(nodes, edges)

    return g_updated


def create_unique_edge_id(u, v, osmid, k=None):
    """
    Create a unique edge ID by combining start node, end node, and osmid.

    Parameters
    ----------
    u : node ID of the edge's source
    v : node ID of the edge's target
    osmid : original OSM way ID
    k : optional key for MultiDiGraphs (default: None)

    Returns
    -------
    str : A unique edge identifier
    """
    # Handle the case where osmid might be a list.
    # Sort before joining so the hash is deterministic regardless of the
    # order osmnx uses when merging edges during simplify_graph.
    if isinstance(osmid, list):
        osmid_str = '_'.join(map(str, sorted(osmid)))
    else:
        osmid_str = str(osmid)

    # Include the key if provided (for MultiDiGraphs)
    if k is not None:
        unique_id = f"{u}_{v}_{k}_{osmid_str}"
    else:
        unique_id = f"{u}_{v}_{osmid_str}"

    # Optionally hash it if you want a shorter fixed-length ID
    hash_object = hashlib.md5(unique_id.encode())
    return hash_object.hexdigest()[:12]  # 12 characters should be sufficient


def validate_graph_topology(G):
    """
    Validate graph topology and fix common issues.

    Fixes:
    - Self-loop edges (u == v)
    - Isolated nodes (nodes with no edges)
    - Checks for duplicate edge IDs

    Parameters
    ----------
    G : networkx.MultiDiGraph
        Input graph

    Returns
    -------
    G : networkx.MultiDiGraph
        Validated and fixed graph
    """
    original_nodes = G.number_of_nodes()
    original_edges = G.number_of_edges()
    self_loops_removed = 0
    protected_self_loops_retained = 0
    isolated_nodes_removed = 0

    # 1. Remove self-loops
    removable_self_loops, retained_self_loops = _split_self_loops_for_cleanup(G)
    if removable_self_loops:
        logger.info(
            "Found %d removable self-loop edges - removing non-protected loops",
            len(removable_self_loops),
        )
        G.remove_edges_from(removable_self_loops)
        self_loops_removed = len(removable_self_loops)
    if retained_self_loops:
        protected_self_loops_retained = len(retained_self_loops)
        logger.warning(
            "Retaining %d protected self-loop edges on major/backbone infrastructure",
            protected_self_loops_retained,
        )

    # 2. Remove isolated nodes
    isolated = list(nx.isolates(G))
    if isolated:
        logger.info("Found %d isolated nodes - removing", len(isolated))
        G.remove_nodes_from(isolated)
        isolated_nodes_removed = len(isolated)

    # Early-exit if the graph is now empty — nothing left to validate.
    if G.number_of_nodes() == 0:
        raise ValueError(
            "validate_graph_topology: graph is empty after cleanup.\n"
            f"  Original : {original_nodes} nodes, {original_edges} edges\n"
            f"  Removed  : {self_loops_removed} self-loop edge(s), "
            f"{isolated_nodes_removed} isolated node(s)\n"
            f"  Remaining: 0 nodes, 0 edges\n"
            "This usually means the input graph consisted entirely of self-loops "
            "or isolated nodes and is not a valid road network."
        )

    # 3. Check for and resolve duplicate edge IDs.
    #    Duplicates occur when the same OSM way appears in multiple layers
    #    (e.g. a primary road covered by both main and residential downloads).
    #    Joining on a non-unique edge_id silently duplicates rows, so we make
    #    them unique by appending a counter suffix to all but the first occurrence.
    nodes, edges = ox.graph_to_gdfs(G)
    id_counts = edges['edge_id'].value_counts()
    duplicates = id_counts[id_counts > 1]

    if len(duplicates) > 0:
        logger.warning("Found %d duplicate edge IDs - appending suffix to deduplicate", len(duplicates))
        seen: dict = {}
        new_ids = []
        for eid in edges['edge_id']:
            count = seen.get(eid, 0)
            new_ids.append(eid if count == 0 else f"{eid}_{count}")
            seen[eid] = count + 1
        edges['edge_id'] = new_ids
        G = ox.graph_from_gdfs(nodes, edges)

    # 4. Detect node pairs that are suspiciously close in physical space.
    #    After consolidation, any remaining near-duplicate nodes (different IDs,
    #    nearly identical coordinates) can still cause R5 to produce self-loops
    #    if its internal vertex snapping collapses them to the same location.
    #    This check projects to UTM for metric distances and uses a spatial index
    #    to avoid O(n²) comparisons.
    CLOSE_NODE_THRESHOLD_M = 1.0  # warn if two nodes are within this many meters
    try:
        from shapely.geometry import Point
        from shapely.strtree import STRtree

        g_proj = ox.project_graph(G)
        nodes_proj, _ = ox.graph_to_gdfs(g_proj)
        node_ids = list(nodes_proj.index)
        points = [Point(geom.x, geom.y) for geom in nodes_proj.geometry]
        tree = STRtree(points)

        close_pairs = []
        for i, pt in enumerate(points):
            for j in tree.query(pt.buffer(CLOSE_NODE_THRESHOLD_M)):
                if j > i:
                    close_pairs.append((node_ids[i], node_ids[j]))

        if close_pairs:
            logger.warning(
                "Found %d node pairs within %.1fm of each other - possible residual "
                "consolidation artifact. Examples: %s",
                len(close_pairs), CLOSE_NODE_THRESHOLD_M, close_pairs[:5]
            )
        else:
            logger.info("No suspiciously close node pairs found (threshold=%.1fm)", CLOSE_NODE_THRESHOLD_M)
    except Exception as e:
        logger.debug("Proximity check skipped: %s", e)

    # Summary if changes were made
    if self_loops_removed > 0 or isolated_nodes_removed > 0 or protected_self_loops_retained > 0:
        final_nodes = G.number_of_nodes()
        final_edges = G.number_of_edges()
        logger.info(
            "Network validation: %d->%d nodes, %d->%d edges (removed_self_loops=%d, retained_protected_self_loops=%d)",
            original_nodes,
            final_nodes,
            original_edges,
            final_edges,
            self_loops_removed,
            protected_self_loops_retained,
        )

    return G


def summarize_edge_quality(G):
    """Summarize edge-level quality checks and anomaly counts."""
    _, edges = ox.graph_to_gdfs(G)

    valid_oneway_values = {"yes", "no", "-1"}
    oneway_series = edges["oneway"].astype(str) if "oneway" in edges.columns else None
    speed_series = edges["speed_kph"] if "speed_kph" in edges.columns else None
    length_series = edges["length"] if "length" in edges.columns else None

    return {
        "missing_edge_id": int(edges["edge_id"].isna().sum()) if "edge_id" in edges.columns else len(edges),
        "missing_length": int(length_series.isna().sum()) if length_series is not None else len(edges),
        "nonpositive_length": int((length_series <= 0).sum()) if length_series is not None else len(edges),
        "short_links_lt_10m": int((length_series < 10).sum()) if length_series is not None else 0,
        "very_short_links_lt_5m": int((length_series < 5).sum()) if length_series is not None else 0,
        "long_links_gt_10km": int((length_series > 10_000).sum()) if length_series is not None else 0,
        "missing_speed_kph": int(speed_series.isna().sum()) if speed_series is not None else len(edges),
        "speed_kph_min": float(speed_series.min()) if speed_series is not None and len(speed_series) else None,
        "speed_kph_max": float(speed_series.max()) if speed_series is not None and len(speed_series) else None,
        "invalid_oneway_values": int((~oneway_series.isin(valid_oneway_values)).sum())
        if oneway_series is not None else len(edges),
        "missing_geometry": int(edges.geometry.isna().sum()) if "geometry" in edges.columns else len(edges),
    }


def duplicate_coords_at_precision(G, precision: int = 7) -> list[tuple]:
    """Return node-ID pairs that collapse to the same rounded coordinate pair."""
    nodes_gdf, _ = ox.graph_to_gdfs(G)
    coord_map = {}
    dupes = []
    for node_id, row in nodes_gdf.iterrows():
        key = (round(row["x"], precision), round(row["y"], precision))
        if key in coord_map:
            dupes.append((coord_map[key], node_id))
        else:
            coord_map[key] = node_id
    return dupes


def nodes_within_distance(G, threshold_m: float) -> list[tuple]:
    """Return node-ID pairs that are closer than the given threshold in meters."""
    graph_projected = ox.project_graph(G)
    nodes_proj, _ = ox.graph_to_gdfs(graph_projected)
    node_ids = list(nodes_proj.index)
    points = [row.geometry for _, row in nodes_proj.iterrows()]
    tree = shapely.STRtree(points)

    close = []
    for i, pt in enumerate(points):
        for j in tree.query(pt.buffer(threshold_m)):
            if j > i:
                close.append((node_ids[i], node_ids[j]))
    return close


def summarize_graph_validation(G, close_threshold_m: float = 0.5) -> dict:
    """Summarize graph-level validation checks and anomaly counts."""
    _, edges = ox.graph_to_gdfs(G)
    highway_counts = (
        edges["highway"].explode().astype(str).value_counts().head(8).to_dict()
        if "highway" in edges.columns
        else {}
    )
    duplicate_pairs = duplicate_coords_at_precision(G)
    close_pairs = nodes_within_distance(G, threshold_m=close_threshold_m)
    self_loops = list(nx.selfloop_edges(G))
    isolates = list(nx.isolates(G))

    summary = {
        "nodes": G.number_of_nodes(),
        "edges": G.number_of_edges(),
        "self_loops": len(self_loops),
        "protected_self_loops": sum(1 for u, v, k in self_loops if _edge_is_protected(G[u][v][k])),
        "unprotected_self_loops": sum(1 for u, v, k in self_loops if not _edge_is_protected(G[u][v][k])),
        "isolated_nodes": len(isolates),
        "weakly_connected": nx.is_weakly_connected(G),
        "duplicate_xml_coordinate_pairs": len(duplicate_pairs),
        "close_node_pairs_lt_0_5m": len(close_pairs),
        "duplicate_examples": duplicate_pairs[:5],
        "close_examples": close_pairs[:5],
        "highway_type_counts": highway_counts,
    }
    summary.update(summarize_edge_quality(G))
    return summary


def format_validation_summary(summary: dict, title: str = "Validation checks") -> str:
    """Render a readable validation report for users and tests."""
    lines = [
        f"{title}:",
        f"- nodes: {summary['nodes']}",
        f"- edges: {summary['edges']}",
        f"- self-loops: {summary['self_loops']}",
        f"- protected self-loops: {summary['protected_self_loops']}",
        f"- unprotected self-loops: {summary['unprotected_self_loops']}",
        f"- isolated nodes: {summary['isolated_nodes']}",
        f"- weakly connected: {summary['weakly_connected']}",
        f"- duplicate XML-rounded coordinate pairs: {summary['duplicate_xml_coordinate_pairs']}",
        f"- node pairs within 0.5m: {summary['close_node_pairs_lt_0_5m']}",
        f"- missing edge_id values: {summary['missing_edge_id']}",
        f"- missing length values: {summary['missing_length']}",
        f"- nonpositive lengths: {summary['nonpositive_length']}",
        f"- short links <10m: {summary['short_links_lt_10m']}",
        f"- very short links <5m: {summary['very_short_links_lt_5m']}",
        f"- long links >10km: {summary['long_links_gt_10km']}",
        f"- missing speed_kph values: {summary['missing_speed_kph']}",
        f"- invalid oneway values: {summary['invalid_oneway_values']}",
        f"- missing edge geometry: {summary['missing_geometry']}",
    ]
    if summary.get("speed_kph_min") is not None:
        lines.append(
            f"- speed_kph min/max: {summary['speed_kph_min']:.2f}/{summary['speed_kph_max']:.2f}"
        )
    if summary.get("highway_type_counts") is not None:
        lines.append(f"- highway type counts: {summary['highway_type_counts']}")
    if "xml_duplicate_coordinate_pairs" in summary:
        lines.append(
            f"- XML duplicate coordinate pairs: {summary['xml_duplicate_coordinate_pairs']}"
        )
    if "xml_dangling_nd_refs" in summary:
        lines.append(f"- XML dangling nd refs: {summary['xml_dangling_nd_refs']}")
    if summary.get("duplicate_examples"):
        lines.append(f"- duplicate coordinate examples: {summary['duplicate_examples']}")
    if summary.get("close_examples"):
        lines.append(f"- close node examples: {summary['close_examples']}")
    if summary.get("xml_dangling_examples"):
        lines.append(f"- XML dangling nd examples: {summary['xml_dangling_examples']}")
    return "\n".join(lines)


def _log_bridge_tunnel_stats(G: nx.MultiDiGraph, stage: str) -> None:
    """Log the number of bridge/tunnel edges at a given pipeline stage.

    Used to track when bridges are severed during consolidation, simplification,
    or topology cleanup.  Comparing counts across stages pinpoints the culprit.
    """
    bridge_count = 0
    tunnel_count = 0
    for u, v, k, data in G.edges(data=True, keys=True):
        if _is_truthy_osm_tag(data.get("bridge")):
            bridge_count += 1
        if _is_truthy_osm_tag(data.get("tunnel")):
            tunnel_count += 1
    logger.info(
        "  [%s] bridge edges: %d, tunnel edges: %d",
        stage, bridge_count, tunnel_count,
    )


def _warn_if_major_infrastructure_removed(G: nx.MultiDiGraph, removed_node_ids: set) -> None:
    """Log per-fragment detail for every disconnected component that was removed.

    Called after the largest-connected-component extraction.  For each fragment
    we report node/edge counts, highway type breakdown, bridge/tunnel edge count,
    total edge length, and geographic bbox — enough to judge whether a fragment
    represents a meaningful road segment or a trivial isolated stub.
    """
    major_types = {"motorway", "motorway_link", "trunk", "trunk_link", "primary", "primary_link"}

    # Partition removed nodes into their own weakly-connected sub-components.
    removed_subgraph = G.subgraph(removed_node_ids)
    fragments = sorted(
        nx.weakly_connected_components(removed_subgraph), key=len, reverse=True
    )

    total_major_edges = 0
    total_bt_edges = 0

    for i, comp in enumerate(fragments):
        sub = G.subgraph(comp)
        hw_counts: dict[str, int] = {}
        bt_edges = 0
        total_length = 0.0
        lons, lats = [], []

        for n in comp:
            nd = G.nodes[n]
            x, y = nd.get("x"), nd.get("y")
            if x is not None and y is not None:
                lons.append(x)
                lats.append(y)

        for u, v, k, data in sub.edges(data=True, keys=True):
            hw = data.get("highway", "unknown")
            if isinstance(hw, list):
                hw = hw[0] if hw else "unknown"
            hw = str(hw).strip()
            hw_counts[hw] = hw_counts.get(hw, 0) + 1
            if _is_truthy_osm_tag(data.get("bridge")) or _is_truthy_osm_tag(data.get("tunnel")):
                bt_edges += 1
            length = data.get("length")
            if isinstance(length, (int, float)):
                total_length += length

        major_found = {k: v for k, v in hw_counts.items() if k in major_types}
        total_major_edges += sum(major_found.values())
        total_bt_edges += bt_edges

        if not major_found and bt_edges == 0:
            continue  # minor local-road stub — skip verbose logging

        bbox_str = ""
        if lons:
            bbox_str = (
                f" bbox lon [{min(lons):.4f},{max(lons):.4f}]"
                f" lat [{min(lats):.4f},{max(lats):.4f}]"
            )

        logger.warning(
            "  Fragment %d: %d nodes, %d edges, length=%.0fm"
            " | major hw: %s | bridge/tunnel: %d%s",
            i + 1,
            len(comp),
            sub.number_of_edges(),
            total_length,
            major_found,
            bt_edges,
            bbox_str,
        )

    if total_major_edges > 0 or total_bt_edges > 0:
        logger.warning(
            "  Total across all fragments: major-road edges=%d, bridge/tunnel edges=%d",
            total_major_edges, total_bt_edges,
        )
    else:
        logger.debug("  All removed fragments contain only local roads — no major infrastructure affected")


def adjust_and_add_graph(graphs, current_graph):
    # Get nodes and edges of current graph
    current_nodes, current_edges = ox.graph_to_gdfs(current_graph)

    # Collect all unique columns from existing graphs
    existing_columns = set()
    for existing_graph in graphs:
        _, existing_edges = ox.graph_to_gdfs(existing_graph)
        existing_columns.update(existing_edges.columns)

    # Add missing columns to current graph's edges
    # Use None (not "") so pandas .notna() / .isna() checks work correctly
    # on all downstream code that tests for missing values.
    for col in existing_columns:
        if col not in current_edges.columns:
            current_edges[col] = None

    # Also ensure existing graphs have columns from current graph
    current_columns = set(current_edges.columns)
    for i, existing_graph in enumerate(graphs):
        existing_nodes, existing_edges = ox.graph_to_gdfs(existing_graph)

        columns_added = False
        for col in current_columns:
            if col not in existing_edges.columns:
                existing_edges[col] = None
                columns_added = True

        # Only rebuild the graph if columns were added
        if columns_added:
            graphs[i] = ox.graph_from_gdfs(existing_nodes, existing_edges)

    # Add the graph to the list if it has edges
    graphs.append(ox.graph_from_gdfs(current_nodes, current_edges))


def download_and_prepare_osm_network(_network_config: dict, _area_config: dict, _geo_config: dict,
                                     work_dir) -> nx.MultiDiGraph:
    """Download and prepare OSM network based on study area configuration."""
    logger.info("=== Preparing OSM Network: %s ===", _area_config['name'])

    # Apply OSMNX settings
    for setting, value in _network_config["osmnx_settings"].items():
        if setting == "overpass_urls":
            continue
        setattr(ox.settings, setting, value)

    # Extract configuration
    study_area = _area_config['name']
    base_name = f"{work_dir}/geo/{study_area}"
    census_year = _area_config["census_year"]
    state_fips_code = _area_config["state_fips"]
    county_fips_codes = _area_config["county_fips"]
    tolerance = _network_config["tolerance"]
    utm_epsg = _geo_config["utm_epsg"]
    should_strongly_connect = _network_config.get("strongly_connected_components", False)

    # Cache handling
    raw_graph_cache_path = _raw_graph_cache_path(_network_config, _area_config, work_dir)
    g_combined = None

    # Try loading from cache
    if os.path.exists(raw_graph_cache_path):
        logger.info("Loading cached graph from: %s", raw_graph_cache_path)
        try:
            with open(raw_graph_cache_path, 'rb') as f:
                g_combined = pickle.load(f)
            logger.info("  Loaded: %d nodes, %d edges", g_combined.number_of_nodes(), g_combined.number_of_edges())
        except Exception as e:
            logger.warning("  Cache load failed: %s. Downloading fresh.", e)
            g_combined = None

    # Download if cache miss
    if g_combined is None:
        logger.info("Downloading OSM network layers...")
        graphs = []

        # Process each configured layer
        for layer_name, layer_config in _network_config["graph_layers"].items():
            geo_level = layer_config["geo_level"]
            min_density = layer_config.get("min_density_per_km2", 0)
            custom_filter = layer_config["custom_filter"]
            buffer_in_meters = layer_config["buffer_zone_in_meters"]
            layer_role = layer_config.get("layer_role", layer_name)

            logger.info("  Processing layer: %s", layer_name)

            # Collect geographic boundaries
            region_boundary_wgs84 = collect_geographic_boundaries(
                state_fips_code=state_fips_code,
                county_fips_codes=county_fips_codes,
                year=census_year,
                area_name=study_area,
                geo_level=geo_level,
                work_dir=work_dir
            )

            # Configure layer-specific parameters
            if layer_role in {"main", "backbone", "connector"}:
                graph_layer = build_area_mask_geometry(
                    region_boundary_wgs84,
                    include_water=True,
                    buffer_m=buffer_in_meters,
                    buffer_epsg=utm_epsg,
                )
                network_type = layer_config.get("network_type", "drive")
                simplify = layer_config.get("simplify", False)
                retain_all = layer_config.get("retain_all", True)
                truncate_by_edge = layer_config.get("truncate_by_edge", True)

            elif layer_role == "residential":
                if min_density > 0:
                    logger.info("    Filtering by density: %s pop/km²", min_density)

                pop_data = collect_census_data(
                    state_fips_code, county_fips_codes, census_year,
                    census_data_file=f"{base_name}_acs_census_{geo_level}_{census_year}.csv",
                    geo_level=geo_level
                )

                filtered_boundaries = filter_boundaries_by_density(
                    region_boundary_wgs84, pop_data, utm_epsg, geo_level, min_density,
                    density_geo_file=f"{base_name}_{geo_level}_{census_year}_{min_density}ppsk_wgs84.geojson",
                )

                graph_layer = shapely.ops.unary_union([
                    to_convex_hull(geom, utm_epsg, buffer_in_meters)
                    for geom in filtered_boundaries.geometry
                ])
                network_type, simplify, retain_all, truncate_by_edge = "drive", False, True, True

            elif layer_role == "ferry":
                graph_layer = build_area_mask_geometry(
                    region_boundary_wgs84,
                    include_water=True,
                    buffer_m=buffer_in_meters,
                    buffer_epsg=utm_epsg,
                )
                network_type, simplify, retain_all, truncate_by_edge = "all", True, True, False

            else:
                raise ValueError(f"Invalid layer role: {layer_role} for layer {layer_name}")

            # Download OSM network for this layer
            g = _download_with_overpass_fallback(
                lambda: ox.graph_from_polygon(
                    graph_layer,
                    network_type=network_type,
                    simplify=simplify,
                    retain_all=retain_all,
                    truncate_by_edge=truncate_by_edge,
                    custom_filter=custom_filter,
                ),
                _network_config,
                context=f"{study_area}:{layer_name}",
            )
            logger.info("    Downloaded: %d nodes, %d edges", g.number_of_nodes(), g.number_of_edges())

            # Special processing for ferry layer
            if layer_role == "ferry":
                g = process_ferry_edges(g)
                if g.number_of_edges() == 0:
                    logger.info("    No suitable ferry connections found - skipping layer")
                    continue

            layer_nodes, layer_edges = ox.graph_to_gdfs(g)
            layer_edges["source_layer"] = layer_name
            layer_edges["layer_role"] = layer_role
            layer_edges["protected_backbone"] = layer_config.get(
                "protected_backbone",
                layer_role == "backbone",
            )
            g = ox.graph_from_gdfs(layer_nodes, layer_edges)

            # Add to list
            adjust_and_add_graph(graphs, g)

        # Combine all layers
        g_combined = nx.compose_all(graphs)
        logger.info("Combined layers: %d nodes, %d edges", g_combined.number_of_nodes(), g_combined.number_of_edges())

        # Save to cache
        try:
            os.makedirs(os.path.dirname(raw_graph_cache_path), exist_ok=True)
            with open(raw_graph_cache_path, 'wb') as f:
                pickle.dump(g_combined, f)
            logger.info("Cached to: %s", raw_graph_cache_path)
        except Exception as e:
            logger.warning("Cache save failed: %s", e)

    # Process the combined graph
    logger.info("Processing network...")

    # Project to UTM
    g_projected = project_graph(g_combined, to_crs=utm_epsg)
    _log_bridge_tunnel_stats(g_projected, "raw (after projection)")

    # Add speeds and process tags
    g_with_speeds = ox.add_edge_speeds(g_projected)
    g_processed_tags = process_tags(g_with_speeds, _network_config)

    # Consolidate intersections
    g_consolidated = ox.consolidate_intersections(
        g_processed_tags,
        tolerance=tolerance,
        rebuild_graph=True,
        dead_ends=True,
        reconnect_edges=True
    )
    _log_bridge_tunnel_stats(g_consolidated, f"after consolidation (tol={tolerance}m)")
    logger.info(
        "After consolidation: %d nodes, %d edges",
        g_consolidated.number_of_nodes(), g_consolidated.number_of_edges(),
    )

    # Simplify graph
    g_simplified = ox.simplification.simplify_graph(
        g_consolidated,
        # Preserve critical corridor structure at bridge/tunnel and layer
        # transitions instead of merging straight through them.
        edge_attrs_differ=[
            "highway",
            "lanes",
            "maxspeed",
            "bridge",
            "tunnel",
            "layer_role",
            "protected_backbone",
        ],
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
            'bridge': first_valid_value,
            'tunnel': first_valid_value,
            'foot': yes_no_all,
            'bicycle': yes_no_all,
            'sidewalk': first_valid_value,
            'cycleway': first_valid_value,
            'maxheight': min_numeric_or_string,
            'maxwidth': min_numeric_or_string,
            'motor_vehicle': yes_no_all,
            'source_layer': first_valid_value,
            'layer_role': first_valid_value,
            'protected_backbone': bool_any,
        }
    )

    _log_bridge_tunnel_stats(g_simplified, "after simplification")
    logger.info(
        "After simplification: %d nodes, %d edges",
        g_simplified.number_of_nodes(), g_simplified.number_of_edges(),
    )

    nodes, edges = ox.graph_to_gdfs(g_simplified)
    edges['edge_id'] = [
        create_unique_edge_id(u, v, row['osmid'], k)
        for (u, v, k), row in edges.iterrows()
    ]
    g_hashed = ox.graph_from_gdfs(nodes, edges)

    # Project back to WGS84
    g_wgs84 = project_graph(g_hashed, to_latlong=True)

    # Validate topology and fix issues
    g_validated = validate_graph_topology(g_wgs84)
    logger.info(
        "After topology validation: %d nodes, %d edges",
        g_validated.number_of_nodes(), g_validated.number_of_edges(),
    )

    # Extract largest connected component
    logger.info("Extracting largest connected component...")
    if should_strongly_connect:
        largest_component = max(strongly_connected_components(g_validated), key=len)
    else:
        largest_component = max(weakly_connected_components(g_validated), key=len)

    # Report if nodes were removed, and warn if major infrastructure is affected.
    removed_node_ids = set(g_validated.nodes()) - set(largest_component)
    removed_nodes = len(removed_node_ids)
    if removed_nodes > 0:
        logger.info("  Removed %d nodes in disconnected components", removed_nodes)
        _warn_if_major_infrastructure_removed(g_validated, removed_node_ids)

    # Create final graph
    g_final = nx.MultiDiGraph(g_validated.subgraph(largest_component).copy())

    logger.info("=== Final network: %d nodes, %d edges ===", g_final.number_of_nodes(), g_final.number_of_edges())

    return g_final
