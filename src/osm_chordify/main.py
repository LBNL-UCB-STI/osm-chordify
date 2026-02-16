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
import warnings

import geopandas as gpd
import matplotlib.pyplot as plt
import networkx as nx
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

    print(f'Downloading and preparing OSM-based {osm_name} network...')
    g_osm = download_and_prepare_osm_network(
        osm_config, area_config, geo_config, work_dir,
    )

    has_invalid, invalid_nodes = check_invalid_coordinates(g_osm)
    if has_invalid:
        print(f"WARNING: Found {len(invalid_nodes)} nodes with invalid coordinates.")
    else:
        print("All node coordinates are valid.")

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
    network_path, intersection_path, network_osm_id_col="attributeOrigId",
    output_path=None,
):
    """
    Map OSM edges to a BEAM network via a zone-network intersection.

    Joins the BEAM network table to the intersection result on the shared
    OSM ID.  All columns from both inputs are included in the output.

    Args:
        network_path (str): Path to BEAM network CSV/CSV.GZ file.
        intersection_path (str): Path to zone-network intersection file.
        network_osm_id_col (str): Column in the BEAM network that holds
            the OSM ID for joining.  Default ``"attributeOrigId"``.
        output_path (str, optional): Path to save the output file.

    Returns:
        pd.DataFrame: Merged result with all columns from both inputs.
    """
    logger.info("Loading network from %s", network_path)
    network_df = pd.read_csv(network_path)

    logger.info("Loading intersection from %s", intersection_path)
    intersection_gdf = gpd.read_file(intersection_path)

    merged_df = map_network_to_intersection(
        network_df=network_df,
        intersection_gdf=intersection_gdf,
        network_osm_id_col=network_osm_id_col,
    )

    if output_path:
        save_dataframe(merged_df, output_path)
    logger.info("Network mapping complete")
    return merged_df


def match_road_network_geometries(
    network_a,
    network_b,
    epsg_utm=None,
    matching="flexible",
    output_path=None,
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
    network_b : str, os.PathLike, or gpd.GeoDataFrame
        Second road network to match against.
    epsg_utm : int, optional
        EPSG code for the projected CRS used for distance calculations.
        Required when the input data is not already in a projected CRS.
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
    from pyrosm import OSM

    warnings.filterwarnings('ignore', category=FutureWarning)

    pbf_path = os.path.expanduser(pbf_path)
    osm = OSM(pbf_path)

    # Get driving network edges (and nodes implicitly)
    try:
        nodes, edges = osm.get_network(network_type="driving", nodes=True)
        print(f"Retrieved {len(nodes)} nodes and {len(edges)} edges from OSM")
    except Exception:
        edges = osm.get_network(network_type="driving")
        print(f"Retrieved {len(edges)} edges from OSM (nodes not separately available)")

    # Project to a suitable CRS for length calculation
    edges_projected = edges.to_crs(f'EPSG:{epsg_utm}')
    edges_projected['length_m'] = edges_projected.geometry.length

    # -----------------------------------------------------------------
    # Connected components
    # -----------------------------------------------------------------
    print("\n" + "=" * 80)
    print("BUILDING CONNECTIVITY GRAPH")
    print("=" * 80)

    G = nx.Graph()

    print("Extracting node coordinates from edge geometries...")
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
            print(f"  Processed {edge_count} edges...")

    print(f"Finished processing {edge_count} edges")
    print(f"Created {len(node_coords_to_id)} unique nodes from edge endpoints")

    components = list(nx.connected_components(G))
    num_components = len(components)

    print("\n" + "=" * 80)
    print("CONNECTED COMPONENT ANALYSIS")
    print("=" * 80)
    print(f"Number of connected components (undirected): {num_components}")

    components_sorted = sorted(components, key=len, reverse=True)

    largest = components_sorted[0]
    largest_size = len(largest)
    total_nodes = G.number_of_nodes()
    total_edges = G.number_of_edges()

    largest_subgraph = G.subgraph(largest)
    largest_edges = largest_subgraph.number_of_edges()

    print(f"Total nodes: {total_nodes}")
    print(f"Total edges: {total_edges}")
    print(f"Largest component nodes: {largest_size} "
          f"({largest_size / total_nodes * 100:.2f}% of all nodes)")
    print(f"Largest component edges: {largest_edges} "
          f"({largest_edges / total_edges * 100:.2f}% of all edges)")

    top_k = min(10, num_components)
    top_sizes = [len(c) for c in components_sorted[:top_k]]
    print(f"\nSizes of top {top_k} components (in nodes): {top_sizes}")

    if num_components > 1:
        print("\nWARNING: Network is not fully connected.")
        print(f"  - There are {num_components} disconnected components")
        print("  - Consider cleaning/removing small components or verifying "
              "that activity locations lie on the largest component.")

        small_components = components_sorted[1:]
        num_small = len(small_components)
        total_small_nodes = sum(len(c) for c in small_components)
        print(f"\nSmall component summary:")
        print(f"  - {num_small} small components contain {total_small_nodes} nodes total")
        print(f"  - Smallest component has {top_sizes[-1]} nodes")
        if num_components <= 20:
            print(f"  - All component sizes: {top_sizes}")
    else:
        print("\nNetwork appears fully connected (single component).")

    # -----------------------------------------------------------------
    # Link-length checks
    # -----------------------------------------------------------------

    short_links = edges_projected[edges_projected['length_m'] < 15].copy()
    short_links = short_links.sort_values('length_m')

    long_links = edges_projected[edges_projected['length_m'] > 10000].copy()
    long_links = long_links.sort_values('length_m', ascending=False)

    print(f"\n{'=' * 80}")
    print(f"Found {len(short_links)} links under 15 meters:\n")

    for idx, row in short_links.iterrows():
        highway = str(row.get('highway', 'N/A'))
        osm_id = str(row.get('id', 'N/A'))
        name = str(row.get('name', 'Unnamed'))

        print(
            f"Length: {row['length_m']:.2f}m | "
            f"Highway: {highway:15s} | "
            f"OSM ID: {osm_id:12s} | "
            f"Name: {name}"
        )

    print(f"\n{'=' * 80}")
    print(f"Found {len(long_links)} links longer than 10 km:\n")

    for idx, row in long_links.iterrows():
        highway = str(row.get('highway', 'N/A'))
        osm_id = str(row.get('id', 'N/A'))
        name = str(row.get('name', 'Unnamed'))

        print(
            f"Length: {row['length_m'] / 1000:.2f}km | "
            f"Highway: {highway:15s} | "
            f"OSM ID: {osm_id:12s} | "
            f"Name: {name}"
        )

    print(f"\n{'=' * 80}")
    print(f"Statistics for links under 15 meters:")
    if len(short_links) > 0:
        print(f"Minimum length: {short_links['length_m'].min():.2f} meters")
        print(f"Maximum length: {short_links['length_m'].max():.2f} meters")
        print(f"Average length: {short_links['length_m'].mean():.2f} meters")
        print(f"Median length: {short_links['length_m'].median():.2f} meters")
    else:
        print("No links shorter than 15 meters.")

    print(f"\nStatistics for entire network (links <= 500m):")
    edges_no_outliers = edges_projected[edges_projected['length_m'] <= 500].copy()
    print(f"Total links: {len(edges_no_outliers)}")
    if len(edges_no_outliers) > 0:
        print(f"Minimum length: {edges_no_outliers['length_m'].min():.2f} meters")
        print(f"Maximum length: {edges_no_outliers['length_m'].max() / 1000:.2f} km")
        print(f"Average length: {edges_no_outliers['length_m'].mean():.2f} meters")
        print(f"Median length: {edges_no_outliers['length_m'].median():.2f} meters")
    else:
        print("No links with length <= 500m to report.")

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
        print(f"\nHistogram saved to 'all_links_histogram.png'")
        plt.show()
    else:
        print("\nSkipping histogram: no links after outlier filtering.")


# ---------------------------------------------------------------------------