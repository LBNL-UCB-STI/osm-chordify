#!/usr/bin/env python3
"""
Download and prepare OSM network data.

@author: haitamlaarabi, cristian.poliziani, zaneedell
"""
import os

from osm_chordify.osm.graph import download_and_prepare_osm_network
from osm_chordify.osm.diagnostics import check_invalid_coordinates, scan_network_directories_for_ways
from osm_chordify.osm.export import export_network
from osm_chordify.study_area_config import get_area_config, generate_network_name


def download_and_build_network(study_area_config):
    """Download, process, validate, and export an OSM network.

    This is the high-level script entrypoint that orchestrates:
    1. Downloading and processing the network via ``download_and_prepare_osm_network``
    2. Validating coordinates
    3. Exporting to all formats via ``export_network``

    Parameters
    ----------
    study_area_config : dict
        Full study-area configuration (as returned by ``get_area_config``).
    """
    config_name = generate_network_name(study_area_config)
    work_dir = study_area_config["work_dir"]
    network_dir = f'{work_dir}/network/{config_name}'

    # Step 1: Download and process
    print(f'Downloading and preparing OSM-based {config_name} network...')
    g_network = download_and_prepare_osm_network(
        study_area_config["network"],
        study_area_config["area"],
        study_area_config["geo"],
        work_dir
    )

    # Step 2: Validate
    has_invalid, invalid_nodes = check_invalid_coordinates(g_network)
    if has_invalid:
        print(f"WARNING: Found {len(invalid_nodes)} nodes with invalid coordinates.")
    else:
        print("All node coordinates are valid.")

    # Step 3: Export all formats
    export_network(
        g_network,
        output_dir=network_dir,
        name=config_name,
        edge_tag_aggs=[('length', 'sum')],
    )


def main():
    area = "sfbay"  # Options: sfbay, seattle
    min_density_per_km2 = 5500  # 5500 for sfbay, 120 for seattle
    strongly_connected_components = False

    # Update study area configuration
    study_area_config = get_area_config(area)
    study_area_config["network"]["strongly_connected_components"] = strongly_connected_components
    study_area_config["network"]["graph_layers"]["residential"]["min_density_per_km2"] = min_density_per_km2
    download_and_build_network(study_area_config)

    # Scan network directories for ways
    study_area_config = get_area_config(area)
    scan_network_directories_for_ways(os.path.expanduser(f'{study_area_config["work_dir"]}/network'))


if __name__ == "__main__":
    main()
