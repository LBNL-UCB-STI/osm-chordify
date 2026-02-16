"""Generic polygon–OSM-edge intersection.

This module provides two functions:

* :func:`load_osm_edges` — loads an OSM GeoJSON + GPKG pair, parses tags,
  and returns a single GeoDataFrame with ``osm_id``, ``edge_id``,
  ``edge_length``, and ``geometry``.

* :func:`intersect_osm_with_zones` — intersects *any* polygon grid
  (ISRM, TAZ, census tracts, …) with OSM edges and returns proportional
  lengths per polygon.
"""

import logging

import geopandas as gpd
import pandas as pd
from tqdm import tqdm

from osm_chordify.osm.tags import extract_tag_as_float, parse_other_tags

logger = logging.getLogger(__name__)

WGS84_EPSG = 4326


def load_osm_edges(osm_geojson_path, osm_gpkg_path):
    """Load OSM edges from a GeoJSON + GPKG pair.

    The GeoJSON file supplies ``osm_id`` and ``other_tags`` (from which
    ``edge_id`` and ``edge_length`` are extracted).  The GPKG file supplies
    the ``edge_id`` ↔ geometry mapping.

    Parameters
    ----------
    osm_geojson_path : str
        Path to the ``*.osm.geojson`` file.
    osm_gpkg_path : str
        Path to the ``*.gpkg`` network file (must contain an ``edges`` layer).

    Returns
    -------
    gpd.GeoDataFrame
        Columns: ``osm_id``, ``edge_id``, ``edge_length``, ``geometry``.

    Raises
    ------
    ValueError
        If required columns are missing from the input files.
    """
    # Load GPKG edges
    gpkg_gdf = gpd.read_file(osm_gpkg_path, layer="edges")
    if "edge_id" not in gpkg_gdf.columns:
        raise ValueError("OSM GPKG file is missing 'edge_id' column")

    # Load GeoJSON
    osm_gdf = gpd.read_file(osm_geojson_path)
    for col in ("osm_id", "other_tags"):
        if col not in osm_gdf.columns:
            raise ValueError(f"OSM GeoJSON file is missing '{col}' column")

    osm_gdf["osm_id"] = osm_gdf["osm_id"].astype(int)
    osm_gdf["parsed_tags"] = osm_gdf["other_tags"].apply(parse_other_tags)
    osm_gdf["edge_id"] = osm_gdf["parsed_tags"].apply(
        lambda x: x.get("edge_id", None)
    )
    osm_gdf["edge_length"] = osm_gdf["parsed_tags"].apply(
        extract_tag_as_float, key="length"
    )

    valid = osm_gdf.dropna(subset=["edge_length"])
    logger.info("Found %d edges with valid length information", len(valid))

    edges_gdf = pd.merge(
        valid[["osm_id", "edge_id", "edge_length"]],
        gpkg_gdf[["edge_id", "geometry"]],
        on="edge_id",
        how="inner",
    )
    edges_gdf = gpd.GeoDataFrame(edges_gdf, geometry="geometry", crs=gpkg_gdf.crs)
    logger.info("Mapped %d edges to geometries", len(edges_gdf))
    return edges_gdf


def intersect_osm_with_zones(
    polygons_gdf,
    edges_gdf,
    polygon_id_col,
    epsg_utm,
    copy_edge_attrs=False,
    copy_polygon_attrs=False,
):
    """Intersect polygon geometries with OSM edge geometries.

    All geometric operations are performed in the supplied UTM projection;
    the returned GeoDataFrame is in WGS 84.

    Parameters
    ----------
    polygons_gdf : gpd.GeoDataFrame
        Polygon grid.  Must contain *polygon_id_col* and a geometry column.
    edges_gdf : gpd.GeoDataFrame
        OSM edges (e.g. from :func:`load_osm_edges`).  Must contain
        ``osm_id``, ``edge_id``, ``edge_length``, and ``geometry``.
    polygon_id_col : str
        Name of the unique-ID column in *polygons_gdf* (e.g. ``"isrm"``,
        ``"taz_id"``).
    epsg_utm : int
        EPSG code for the UTM projection used for length calculations.
    copy_edge_attrs : bool, optional
        If ``True``, carry over all extra edge attributes (prefixed with
        ``edge_``).
    copy_polygon_attrs : bool, optional
        If ``True``, carry over all extra polygon attributes (prefixed with
        ``polygon_``).

    Returns
    -------
    gpd.GeoDataFrame
        Columns: ``polygon_id``, ``osm_id``, ``edge_id``,
        ``original_edge_length``, ``proportion``, ``proportional_length``,
        ``geometry``.  CRS is WGS 84.

    Raises
    ------
    ValueError
        If *polygon_id_col* is not found in *polygons_gdf* or required
        edge columns are missing.
    """
    if polygon_id_col not in polygons_gdf.columns:
        raise ValueError(
            f"Polygon GeoDataFrame is missing '{polygon_id_col}' column"
        )
    for col in ("osm_id", "edge_id", "edge_length"):
        if col not in edges_gdf.columns:
            raise ValueError(f"Edges GeoDataFrame is missing '{col}' column")

    # Project to UTM
    logger.info("Projecting geometries to EPSG:%d", epsg_utm)
    edges_utm = edges_gdf.to_crs(epsg=epsg_utm)
    polys_utm = polygons_gdf.to_crs(epsg=epsg_utm)

    sindex = edges_utm.sindex
    skip_edge_keys = {"geometry", "osm_id", "edge_length", "edge_id"}
    skip_poly_keys = {"geometry", polygon_id_col}

    results = []
    logger.info("Intersecting %d polygons with %d edges", len(polys_utm), len(edges_utm))

    for _, poly_row in tqdm(
        polys_utm.iterrows(), total=len(polys_utm), desc="Processing polygons"
    ):
        poly_id = poly_row[polygon_id_col]
        poly_geom = poly_row.geometry

        candidates_idx = list(sindex.intersection(poly_geom.bounds))
        if not candidates_idx:
            continue

        candidates = edges_utm.iloc[candidates_idx]
        intersecting = candidates[candidates.geometry.intersects(poly_geom)]
        if intersecting.empty:
            continue

        for _, edge_row in intersecting.iterrows():
            edge_geom = edge_row.geometry
            intersection_geom = edge_geom.intersection(poly_geom)
            if intersection_geom.is_empty:
                continue

            geom_length = edge_geom.length
            proportion = (
                round(intersection_geom.length / geom_length, 2)
                if geom_length > 0
                else 0
            )
            original_length = edge_row["edge_length"]

            record = {
                "polygon_id": poly_id,
                "osm_id": edge_row["osm_id"],
                "edge_id": edge_row["edge_id"],
                "original_edge_length": original_length,
                "proportion": proportion,
                "proportional_length": original_length * proportion,
                "geometry": intersection_geom,
            }

            if copy_edge_attrs:
                for k, v in edge_row.items():
                    if k not in skip_edge_keys and k not in record:
                        record[f"edge_{k}"] = v

            if copy_polygon_attrs:
                for k, v in poly_row.items():
                    if k not in skip_poly_keys and k not in record:
                        record[f"polygon_{k}"] = v

            results.append(record)

    logger.info("Intersection produced %d results", len(results))

    if not results:
        logger.warning("No intersections found — returning empty GeoDataFrame")
        return gpd.GeoDataFrame(
            columns=[
                "polygon_id", "osm_id", "edge_id",
                "original_edge_length", "proportion",
                "proportional_length", "geometry",
            ],
            geometry="geometry",
            crs=WGS84_EPSG,
        )

    result_gdf = gpd.GeoDataFrame(results, geometry="geometry", crs=epsg_utm)
    result_gdf = result_gdf.to_crs(epsg=WGS84_EPSG)
    logger.info("Converted results back to WGS 84")
    return result_gdf
