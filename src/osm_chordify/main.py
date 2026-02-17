#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Orchestration pipeline for OSM network download, export, and polygon
intersection.

Combines network download/build with polygon-grid intersection and
tabular-network mapping into a single module.

@author: haitamlaarabi, cristian.poliziani, zaneedell
"""

import logging
import os

import geopandas as gpd
import pandas as pd

from osm_chordify.osm.diagnostics import check_invalid_coordinates
from osm_chordify.osm.export import export_network
from osm_chordify.osm.graph import download_and_prepare_osm_network
from osm_chordify.utils.geo import name_osm_network
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
    area_name = area_config["name"]
    osm_name = name_osm_network(
        area_name,
        osm_config["graph_layers"],
        osm_config.get("strongly_connected_components", False),
    )
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

    export_network(
        g_osm,
        output_dir=osm_dir,
        name=osm_name,
        edge_tag_aggs=[('length', 'sum')],
    )


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
        osm_gdf = gpd.read_file(osm_path, layer="lines")
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
    import warnings

    import matplotlib.pyplot as plt
    import networkx as nx
    from pyrosm import OSM

    warnings.filterwarnings('ignore', category=FutureWarning)

    pbf_path = os.path.expanduser(pbf_path)
    osm = OSM(pbf_path)

    # Get driving network edges (and nodes implicitly)
    try:
        nodes, edges = osm.get_network(network_type="driving", nodes=True)
        logger.info("Retrieved %d nodes and %d edges from OSM", len(nodes), len(edges))
    except Exception:
        edges = osm.get_network(network_type="driving")
        logger.info("Retrieved %d edges from OSM (nodes not separately available)", len(edges))

    # Project to a suitable CRS for length calculation
    edges_projected = edges.to_crs(f'EPSG:{epsg_utm}')
    edges_projected['length_m'] = edges_projected.geometry.length

    # -----------------------------------------------------------------
    # Connected components
    # -----------------------------------------------------------------
    logger.info("=" * 80)
    logger.info("BUILDING CONNECTIVITY GRAPH")
    logger.info("=" * 80)

    G = nx.Graph()

    logger.info("Extracting node coordinates from edge geometries...")
    node_coords_to_id = {}
    next_node_id = 0

    edge_count = 0
    for idx, row in edges_projected.iterrows():
        geom = row.geometry
        coords = list(geom.coords)
        start_coord = coords[0]
        end_coord = coords[-1]

        if start_coord not in node_coords_to_id:
            node_coords_to_id[start_coord] = next_node_id
            next_node_id += 1

        if end_coord not in node_coords_to_id:
            node_coords_to_id[end_coord] = next_node_id
            next_node_id += 1

        start_id = node_coords_to_id[start_coord]
        end_id = node_coords_to_id[end_coord]

        G.add_edge(start_id, end_id)
        edge_count += 1

        if edge_count % 10000 == 0:
            logger.info("  Processed %d edges...", edge_count)

    logger.info("Finished processing %d edges", edge_count)
    logger.info("Created %d unique nodes from edge endpoints", len(node_coords_to_id))

    components = list(nx.connected_components(G))
    num_components = len(components)

    logger.info("=" * 80)
    logger.info("CONNECTED COMPONENT ANALYSIS")
    logger.info("=" * 80)
    logger.info("Number of connected components (undirected): %d", num_components)

    components_sorted = sorted(components, key=len, reverse=True)

    largest = components_sorted[0]
    largest_size = len(largest)
    total_nodes = G.number_of_nodes()
    total_edges = G.number_of_edges()

    largest_subgraph = G.subgraph(largest)
    largest_edges = largest_subgraph.number_of_edges()

    logger.info("Total nodes: %d", total_nodes)
    logger.info("Total edges: %d", total_edges)
    logger.info("Largest component nodes: %d (%.2f%% of all nodes)",
                largest_size, largest_size / total_nodes * 100)
    logger.info("Largest component edges: %d (%.2f%% of all edges)",
                largest_edges, largest_edges / total_edges * 100)

    top_k = min(10, num_components)
    top_sizes = [len(c) for c in components_sorted[:top_k]]
    logger.info("Sizes of top %d components (in nodes): %s", top_k, top_sizes)

    if num_components > 1:
        logger.warning("Network is not fully connected.")
        logger.info("  - There are %d disconnected components", num_components)
        logger.info("  - Consider cleaning/removing small components or verifying "
                     "that activity locations lie on the largest component.")

        small_components = components_sorted[1:]
        num_small = len(small_components)
        total_small_nodes = sum(len(c) for c in small_components)
        logger.info("Small component summary:")
        logger.info("  - %d small components contain %d nodes total",
                     num_small, total_small_nodes)
        logger.info("  - Smallest component has %d nodes", top_sizes[-1])
        if num_components <= 20:
            logger.info("  - All component sizes: %s", top_sizes)
    else:
        logger.info("Network appears fully connected (single component).")

    # -----------------------------------------------------------------
    # Link-length checks
    # -----------------------------------------------------------------

    short_links = edges_projected[edges_projected['length_m'] < 15].copy()
    short_links = short_links.sort_values('length_m')

    long_links = edges_projected[edges_projected['length_m'] > 10000].copy()
    long_links = long_links.sort_values('length_m', ascending=False)

    logger.info("=" * 80)
    logger.info("Found %d links under 15 meters:", len(short_links))

    for idx, row in short_links.iterrows():
        highway = str(row.get('highway', 'N/A'))
        osm_id = str(row.get('id', 'N/A'))
        name = str(row.get('name', 'Unnamed'))

        logger.info(
            f"Length: {row['length_m']:.2f}m | "
            f"Highway: {highway:15s} | "
            f"OSM ID: {osm_id:12s} | "
            f"Name: {name}"
        )

    logger.info("=" * 80)
    logger.info("Found %d links longer than 10 km:", len(long_links))

    for idx, row in long_links.iterrows():
        highway = str(row.get('highway', 'N/A'))
        osm_id = str(row.get('id', 'N/A'))
        name = str(row.get('name', 'Unnamed'))

        logger.info(
            f"Length: {row['length_m'] / 1000:.2f}km | "
            f"Highway: {highway:15s} | "
            f"OSM ID: {osm_id:12s} | "
            f"Name: {name}"
        )

    logger.info("=" * 80)
    logger.info("Statistics for links under 15 meters:")
    if len(short_links) > 0:
        logger.info("Minimum length: %.2f meters", short_links['length_m'].min())
        logger.info("Maximum length: %.2f meters", short_links['length_m'].max())
        logger.info("Average length: %.2f meters", short_links['length_m'].mean())
        logger.info("Median length: %.2f meters", short_links['length_m'].median())
    else:
        logger.info("No links shorter than 15 meters.")

    logger.info("Statistics for entire network (links <= 500m):")
    edges_no_outliers = edges_projected[edges_projected['length_m'] <= 500].copy()
    logger.info("Total links: %d", len(edges_no_outliers))
    if len(edges_no_outliers) > 0:
        logger.info("Minimum length: %.2f meters", edges_no_outliers['length_m'].min())
        logger.info("Maximum length: %.2f km", edges_no_outliers['length_m'].max() / 1000)
        logger.info("Average length: %.2f meters", edges_no_outliers['length_m'].mean())
        logger.info("Median length: %.2f meters", edges_no_outliers['length_m'].median())
    else:
        logger.info("No links with length <= 500m to report.")

    # Histogram
    if len(edges_no_outliers) > 0:
        plt.figure(figsize=(12, 6))
        plt.hist(edges_no_outliers['length_m'], bins=50, edgecolor='black', alpha=0.7)
        plt.xlabel('Length (meters)')
        plt.ylabel('Number of Links')
        plt.title('Distribution of All Link Lengths in Network')
        plt.grid(True, alpha=0.3)

        mean_len = edges_no_outliers['length_m'].mean()
        median_len = edges_no_outliers['length_m'].median()

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


# ---------------------------------------------------------------------------