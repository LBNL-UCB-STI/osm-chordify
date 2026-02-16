"""Road-network / zone-polygon intersection.

* :func:`load_osm_edges` — loads edges from an OSM GPKG file.
* :func:`intersect_road_network_with_zones` — intersects road-network edges
  with zone polygons (TAZ, census tracts, ISRM cells, …) and returns
  proportional lengths per zone.
"""

import logging
import os

import geopandas as gpd
import pandas as pd
from tqdm import tqdm

logger = logging.getLogger(__name__)

WGS84_EPSG = 4326


def load_osm_edges(osm_gpkg_path):
    """Load OSM edges from a GPKG file.

    Reads the ``edges`` layer and normalises column names to
    ``osm_id``, ``edge_id``, ``edge_length``, and ``geometry``.

    Parameters
    ----------
    osm_gpkg_path : str
        Path to the ``*.gpkg`` network file (must contain an ``edges`` layer).

    Returns
    -------
    gpd.GeoDataFrame
        Columns: ``osm_id``, ``edge_id``, ``edge_length``, ``geometry``.

    Raises
    ------
    ValueError
        If required columns are missing from the GPKG file.
    """
    edges_gdf = gpd.read_file(osm_gpkg_path, layer="edges")

    for col in ("osmid", "edge_id", "length"):
        if col not in edges_gdf.columns:
            raise ValueError(f"GPKG edges layer is missing '{col}' column")

    edges_gdf = edges_gdf.rename(columns={"osmid": "osm_id", "length": "edge_length"})
    edges_gdf["osm_id"] = edges_gdf["osm_id"].astype(int)

    edges_gdf = edges_gdf[["osm_id", "edge_id", "edge_length", "geometry"]]
    logger.info("Loaded %d edges from %s", len(edges_gdf), osm_gpkg_path)
    return edges_gdf


def _load_edges(road_network):
    """Resolve *road_network* to a GeoDataFrame of edges."""
    if isinstance(road_network, gpd.GeoDataFrame):
        return road_network
    if not isinstance(road_network, (str, os.PathLike)):
        raise TypeError(
            f"road_network must be a GeoDataFrame or a file path, got {type(road_network).__name__}"
        )
    path = str(road_network)
    if path.endswith(".gpkg"):
        return load_osm_edges(path)
    return gpd.read_file(path)


def _load_zones(zones):
    """Resolve *zones* to a GeoDataFrame of polygons."""
    if isinstance(zones, gpd.GeoDataFrame):
        return zones
    if not isinstance(zones, (str, os.PathLike)):
        raise TypeError(
            f"zones must be a GeoDataFrame or a file path, got {type(zones).__name__}"
        )
    return gpd.read_file(str(zones))


def intersect_road_network_with_zones(
    road_network,
    zones,
    epsg_utm,
    proportional_cols=None,
    output_path=None,
):
    """Intersect road-network edges with zone polygons.

    Both *road_network* and *zones* accept either a file path (GeoJSON,
    GPKG, Shapefile, …) or a :class:`~geopandas.GeoDataFrame`.  When a
    ``.gpkg`` path is given for *road_network*, the ``edges`` layer is read
    and columns are normalised automatically via :func:`load_osm_edges`.

    All attributes from both inputs are carried through to the result,
    prefixed with ``edge_`` and ``zone_`` respectively to avoid collisions.

    The ``proportion`` column indicates what fraction of each edge's
    geometry falls within the intersecting zone.  Any columns listed in
    *proportional_cols* are multiplied by this proportion and added as
    new ``proportional_<col>`` columns.

    Parameters
    ----------
    road_network : str, os.PathLike, or gpd.GeoDataFrame
        Road-network edges with line geometries.
    zones : str, os.PathLike, or gpd.GeoDataFrame
        Zone polygons.
    epsg_utm : int
        EPSG code for the UTM projection used for length calculations.
    proportional_cols : str or list of str, optional
        Edge attribute column(s) to scale by the intersection proportion.
        For each column ``col``, a ``proportional_<col>`` column is added
        to the result (e.g. ``proportional_cols=["edge_length", "vmt"]``
        produces ``proportional_edge_length`` and ``proportional_vmt``).
    output_path : str, optional
        If provided, save the result to this file (GeoJSON, GPKG, etc.).

    Returns
    -------
    gpd.GeoDataFrame
        One row per edge-zone intersection piece.  Columns include all
        original edge attributes (prefixed ``edge_``), all original zone
        attributes (prefixed ``zone_``), plus ``edge_length_m``,
        ``proportion``, ``proportional_length_m``, any requested
        ``proportional_*`` columns, and ``geometry``.  CRS is WGS 84.
    """
    edges_gdf = _load_edges(road_network)
    polygons_gdf = _load_zones(zones)

    if proportional_cols is None:
        proportional_cols = []
    elif isinstance(proportional_cols, str):
        proportional_cols = [proportional_cols]

    # Validate that requested columns exist
    for col in proportional_cols:
        if col not in edges_gdf.columns:
            raise ValueError(
                f"proportional_cols: '{col}' not found in road network columns "
                f"{list(edges_gdf.columns)}"
            )

    # Project to UTM
    logger.info("Projecting geometries to EPSG:%d", epsg_utm)
    edges_utm = edges_gdf.to_crs(epsg=epsg_utm)
    polys_utm = polygons_gdf.to_crs(epsg=epsg_utm)

    sindex = edges_utm.sindex
    edge_attr_cols = [c for c in edges_gdf.columns if c != "geometry"]
    zone_attr_cols = [c for c in polygons_gdf.columns if c != "geometry"]

    results = []
    logger.info("Intersecting %d zones with %d edges", len(polys_utm), len(edges_utm))

    for _, poly_row in tqdm(
        polys_utm.iterrows(), total=len(polys_utm), desc="Processing zones"
    ):
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

            edge_length = edge_geom.length
            intersection_length = intersection_geom.length
            proportion = (
                round(intersection_length / edge_length, 4)
                if edge_length > 0
                else 0
            )

            record = {
                "edge_length_m": round(edge_length, 2),
                "proportion": proportion,
                "proportional_length_m": round(intersection_length, 2),
                "geometry": intersection_geom,
            }

            for col in edge_attr_cols:
                record[f"edge_{col}"] = edge_row[col]

            for col in proportional_cols:
                record[f"proportional_{col}"] = edge_row[col] * proportion

            for col in zone_attr_cols:
                record[f"zone_{col}"] = poly_row[col]

            results.append(record)

    logger.info("Intersection produced %d results", len(results))

    if not results:
        logger.warning("No intersections found — returning empty GeoDataFrame")
        return gpd.GeoDataFrame(
            columns=[
                "edge_length_m", "proportion",
                "proportional_length_m", "geometry",
            ],
            geometry="geometry",
            crs=WGS84_EPSG,
        )

    result_gdf = gpd.GeoDataFrame(results, geometry="geometry", crs=epsg_utm)
    result_gdf = result_gdf.to_crs(epsg=WGS84_EPSG)
    logger.info("Converted results back to WGS 84")

    if output_path:
        from osm_chordify.utils.io import save_geodataframe
        save_geodataframe(result_gdf, output_path)
        logger.info("Saved intersection to %s", output_path)

    return result_gdf
