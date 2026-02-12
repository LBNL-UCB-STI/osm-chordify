#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script to intersect polygon grids with OSM edge geometries and then map
BEAM network links to the result.

Polygon grids (ISRM, TAZ, census tracts, …) are declared in the study-area
config under ``polygon_grids``.  Each entry specifies its shapefile path and
the name of its unique-ID column.

This is a thin orchestration script.  The heavy lifting lives in:
- :mod:`osm_chordify.osm.intersect` (generic polygon × edge intersection)
- :mod:`osm_chordify.utils.beam`     (BEAM network mapping)
- :mod:`osm_chordify.utils.io`       (file saving)
"""

import logging
import os

import geopandas as gpd
import pandas as pd

from osm_chordify.osm.intersect import (
    load_osm_edges,
    intersect_edges_with_polygons,
    intersect_network_with_polygons,
)
from osm_chordify.study_area_config import get_area_config, generate_network_name
from osm_chordify.utils.beam import map_beam_network_to_polygon_intersection
from osm_chordify.utils.io import save_dataframe, save_geodataframe

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def process_polygon_osm_intersection(
    polygon_grid_path,
    polygon_id_col,
    osm_geojson_path,
    osm_gpkg_path,
    epsg_utm,
    output_path,
):
    """
    Intersect a polygon grid with OSM edge geometries.

    Loads the polygon grid and OSM edges, delegates to
    :func:`~osm_chordify.osm.intersect.intersect_edges_with_polygons`,
    and saves the result.

    Args:
        polygon_grid_path (str): Path to polygon grid shapefile.
        polygon_id_col (str): Name of the unique-ID column in the grid
            (e.g. ``"isrm"``, ``"taz_id"``).
        osm_geojson_path (str): Path to OSM GeoJSON file.
        osm_gpkg_path (str): Path to OSM GPKG network file.
        epsg_utm (int): EPSG code for UTM projection.
        output_path (str): Path to output file.

    Returns:
        gpd.GeoDataFrame: Intersection results.
    """
    logger.info("Loading polygon grid from %s (id_col=%s)", polygon_grid_path, polygon_id_col)
    polygons_gdf = gpd.read_file(polygon_grid_path)
    if polygon_id_col not in polygons_gdf.columns:
        raise ValueError(
            f"Polygon grid file is missing '{polygon_id_col}' column"
        )

    logger.info("Loading OSM edges")
    edges_gdf = load_osm_edges(osm_geojson_path, osm_gpkg_path)

    result_gdf = intersect_edges_with_polygons(
        polygons_gdf=polygons_gdf,
        edges_gdf=edges_gdf,
        polygon_id_col=polygon_id_col,
        epsg_utm=epsg_utm,
        copy_edge_attrs=True,
        copy_polygon_attrs=True,
    )

    if result_gdf.empty:
        logger.warning("No intersections found")
        return result_gdf

    # Save
    save_geodataframe(result_gdf, output_path)
    logger.info("Polygon-OSM intersection complete")
    return result_gdf


def map_beam_to_polygon_osm(network_path, polygon_osm_path, polygon_id_col, output_path):
    """
    Map BEAM network data to a polygon-OSM intersection result.

    Loads the BEAM CSV and the polygon-OSM file, delegates to
    :func:`~osm_chordify.utils.beam.map_beam_network_to_polygon_intersection`,
    and saves.

    Args:
        network_path (str): Path to BEAM ``network.csv.gz``.
        polygon_osm_path (str): Path to polygon-OSM intersection file.
        polygon_id_col (str): Name of the polygon-ID column in the
            intersection file (must match what was used during intersection).
        output_path (str): Path to save the output file.

    Returns:
        pd.DataFrame: Merged result.
    """
    logger.info("Loading BEAM network from %s", network_path)
    network_df = pd.read_csv(network_path)

    logger.info("Loading polygon-OSM intersection from %s", polygon_osm_path)
    polygon_osm_gdf = gpd.read_file(polygon_osm_path)

    merged_df = map_beam_network_to_polygon_intersection(
        network_df=network_df,
        intersection_gdf=polygon_osm_gdf,
        polygon_id_col=polygon_id_col,
    )

    save_dataframe(merged_df, output_path)
    logger.info("BEAM-polygon mapping complete")
    return merged_df


def main():
    """Main execution function with hardcoded paths."""
    area = "sfbay"
    study_area_config = get_area_config(area)
    network_config = study_area_config["network"]
    geo_config = study_area_config["geo"]
    network_config["graph_layers"]["residential"]["min_density_per_km2"] = 5500

    network_name = generate_network_name(study_area_config)
    work_dir = study_area_config["work_dir"]
    network_dir = f'{work_dir}/network/{network_name}'
    utm_epsg = geo_config["utm_epsg"]

    osm_geojson_path = os.path.expanduser(f"{network_dir}/{network_name}.osm.geojson")
    osm_gpkg_path = os.path.expanduser(f"{network_dir}/{network_name}.gpkg")

    polygon_grid = study_area_config.get("polygon_grid", {})
    if not polygon_grid:
        logger.info("No polygon_grid configured for %s — skipping", area)
        return

    polygon_id_col = polygon_grid["id_col"]
    grid_path = os.path.expanduser(f"{work_dir}/{polygon_grid['path']}")

    out_dir = os.path.expanduser(f"{work_dir}/polygon-{network_name}")
    os.makedirs(out_dir, exist_ok=True)
    intersection_path = os.path.expanduser(
        f"{out_dir}/polygon-{network_name}.geojson"
    )

    if not os.path.exists(intersection_path):
        process_polygon_osm_intersection(
            polygon_grid_path=grid_path,
            polygon_id_col=polygon_id_col,
            osm_geojson_path=osm_geojson_path,
            osm_gpkg_path=osm_gpkg_path,
            epsg_utm=utm_epsg,
            output_path=intersection_path,
        )

    # Network-to-polygon spatial intersection
    polygon_network = study_area_config.get("polygon_network", {})
    if polygon_network and polygon_network.get("path"):
        network_file = os.path.expanduser(
            f"{work_dir}/{polygon_network['path']}"
        )
        network_output_path = os.path.expanduser(
            f"{out_dir}/polygon-network-intersection.geojson"
        )

        if not os.path.exists(network_output_path):
            process_network_polygon_intersection(
                polygon_grid_path=grid_path,
                polygon_id_col=polygon_id_col,
                network_path=network_file,
                network_id_col=polygon_network.get("id_col", "link_id"),
                network_length_col=polygon_network.get("length_col", "length"),
                epsg_utm=utm_epsg,
                output_path=network_output_path,
            )


def process_network_polygon_intersection(
    polygon_grid_path,
    polygon_id_col,
    network_path,
    network_id_col,
    network_length_col,
    epsg_utm,
    output_path,
):
    """
    Spatially intersect network link geometries with a polygon grid.

    Loads the polygon grid and network geo file, delegates to
    :func:`~osm_chordify.osm.intersect.intersect_network_with_polygons`,
    and saves.

    Args:
        polygon_grid_path (str): Path to polygon grid file (shp/gpkg/geojson).
        polygon_id_col (str): Name of the unique-ID column in the grid.
        network_path (str): Path to network geo file with link geometries.
        network_id_col (str): Column name for the link ID in the network.
        network_length_col (str): Column name for the link length in the network.
        epsg_utm (int): EPSG code for UTM projection.
        output_path (str): Path to save the output file.

    Returns:
        gpd.GeoDataFrame: Intersection results.
    """
    logger.info("Loading polygon grid from %s", polygon_grid_path)
    polygons_gdf = gpd.read_file(polygon_grid_path)

    logger.info("Loading network from %s", network_path)
    network_gdf = gpd.read_file(network_path)

    result_gdf = intersect_network_with_polygons(
        polygons_gdf=polygons_gdf,
        network_gdf=network_gdf,
        polygon_id_col=polygon_id_col,
        network_id_col=network_id_col,
        network_length_col=network_length_col,
        epsg_utm=epsg_utm,
    )

    if result_gdf.empty:
        logger.warning("No intersections found")
        return result_gdf

    save_geodataframe(result_gdf, output_path)
    logger.info("Network-polygon intersection complete")
    return result_gdf


if __name__ == "__main__":
    main()
