"""Graph processing: ferry edges, topology, unique IDs, download pipeline."""

import hashlib
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
from osm_chordify.utils.geo import project_graph, to_convex_hull

logger = logging.getLogger(__name__)


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

    # Standardize tags
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

    if 'oneway' in edges.columns:
        edges['oneway'] = edges['oneway'].astype(str).str.lower()
        # Map non-standard values to standard BEAM-compatible values
        edges['oneway'] = edges['oneway'].replace({
            'reverse': '-1',
            'true': 'yes',
            '-1.0': '-1',
            '1.0': 'yes'
        })

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
    # Handle the case where osmid might be a list
    if isinstance(osmid, list):
        osmid_str = '_'.join(map(str, osmid))
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
    isolated_nodes_removed = 0

    # 1. Remove self-loops
    self_loops = list(nx.selfloop_edges(G))
    if self_loops:
        logger.info("Found %d self-loop edges - removing", len(self_loops))
        # Show example with proper MultiDiGraph key handling
        if len(self_loops) > 0:
            # Handle both tuple formats from selfloop_edges
            first_loop = self_loops[0]
            if len(first_loop) == 3:
                u, v, key = first_loop
            else:
                u, v = first_loop
                key = 0

            # Access edge data properly for MultiDiGraph
            if G.is_multigraph():
                edge_data = G[u][v][key]
            else:
                edge_data = G[u][v]

            logger.debug("  Example: Node %s->%s, length=%sm", u, u, edge_data.get('length', 'N/A'))

        # Remove all self-loops
        G.remove_edges_from(self_loops)
        self_loops_removed = len(self_loops)

    # 2. Remove isolated nodes
    isolated = list(nx.isolates(G))
    if isolated:
        logger.info("Found %d isolated nodes - removing", len(isolated))
        G.remove_nodes_from(isolated)
        isolated_nodes_removed = len(isolated)

    # 3. Check for duplicate edge IDs
    nodes, edges = ox.graph_to_gdfs(G)
    id_counts = edges['edge_id'].value_counts()
    duplicates = id_counts[id_counts > 1]

    if len(duplicates) > 0:
        logger.warning("Found %d duplicate edge IDs", len(duplicates))

    # Summary if changes were made
    if self_loops_removed > 0 or isolated_nodes_removed > 0:
        final_nodes = G.number_of_nodes()
        final_edges = G.number_of_edges()
        logger.info("Network validation: %d->%d nodes, %d->%d edges",
                    original_nodes, final_nodes, original_edges, final_edges)

    return G


def adjust_and_add_graph(graphs, current_graph):
    # Get nodes and edges of current graph
    current_nodes, current_edges = ox.graph_to_gdfs(current_graph)

    # Collect all unique columns from existing graphs
    existing_columns = set()
    for existing_graph in graphs:
        _, existing_edges = ox.graph_to_gdfs(existing_graph)
        existing_columns.update(existing_edges.columns)

    # Add missing columns to current graph's edges
    for col in existing_columns:
        if col not in current_edges.columns:
            current_edges[col] = ""

    # Also ensure existing graphs have columns from current graph
    current_columns = set(current_edges.columns)
    for i, existing_graph in enumerate(graphs):
        existing_nodes, existing_edges = ox.graph_to_gdfs(existing_graph)

        columns_added = False
        for col in current_columns:
            if col not in existing_edges.columns:
                existing_edges[col] = ""
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
    raw_graph_cache_path = os.path.join(work_dir, 'network', f'raw_osm_graph_{study_area}.pkl')
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
            if layer_name == "main":
                graph_layer = to_convex_hull(region_boundary_wgs84, utm_epsg, buffer_in_meters)
                network_type, simplify, retain_all, truncate_by_edge = "drive", False, True, True

            elif layer_name == "residential":
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

            elif layer_name == "ferry":
                graph_layer = to_convex_hull(region_boundary_wgs84, utm_epsg, buffer_in_meters)
                network_type, simplify, retain_all, truncate_by_edge = "all", True, True, False

            else:
                raise ValueError(f"Invalid layer name: {layer_name}")

            # Download OSM network for this layer
            g = ox.graph_from_polygon(
                graph_layer,
                network_type=network_type,
                simplify=simplify,
                retain_all=retain_all,
                truncate_by_edge=truncate_by_edge,
                custom_filter=custom_filter
            )
            logger.info("    Downloaded: %d nodes, %d edges", g.number_of_nodes(), g.number_of_edges())

            # Special processing for ferry layer
            if layer_name == "ferry":
                g = process_ferry_edges(g)
                if g.number_of_edges() == 0:
                    logger.info("    No suitable ferry connections found - skipping layer")
                    continue

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

    # Simplify graph
    g_simplified = ox.simplification.simplify_graph(
        g_consolidated,
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
            'bridge': first_valid_value,
            'tunnel': first_valid_value,
            'foot': yes_no_all,
            'bicycle': yes_no_all,
            'sidewalk': first_valid_value,
            'cycleway': first_valid_value,
            'maxheight': min_numeric_or_string,
            'maxwidth': min_numeric_or_string,
            'motor_vehicle': yes_no_all,
        }
    )

    # Create unique edge IDs
    nodes, edges = ox.graph_to_gdfs(g_simplified)
    edges['edge_id'] = edges.apply(
        lambda row: create_unique_edge_id(
            row['u_original'], row['v_original'], row['osmid'], row.get('key', None)
        ),
        axis=1
    )
    g_hashed = ox.graph_from_gdfs(nodes, edges)

    # Project back to WGS84
    g_wgs84 = project_graph(g_hashed, to_latlong=True)

    # Validate topology and fix issues
    g_validated = validate_graph_topology(g_wgs84)

    # Extract largest connected component
    logger.info("Extracting largest connected component...")
    if should_strongly_connect:
        largest_component = max(strongly_connected_components(g_validated), key=len)
    else:
        largest_component = max(weakly_connected_components(g_validated), key=len)

    # Report if nodes were removed
    removed_nodes = g_validated.number_of_nodes() - len(largest_component)
    if removed_nodes > 0:
        logger.info("  Removed %d nodes in disconnected components", removed_nodes)

    # Create final graph
    g_final = nx.MultiDiGraph(g_validated.subgraph(largest_component).copy())

    logger.info("=== Final network: %d nodes, %d edges ===", g_final.number_of_nodes(), g_final.number_of_edges())

    return g_final
