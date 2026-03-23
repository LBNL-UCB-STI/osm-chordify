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
    road_network_epsg,
    zones,
    proportional_cols=None,
    output_path=None,
    output_epsg=None,
):
    """Intersect road-network edges with zone polygons.

    Both *road_network* and *zones* accept either a file path (GeoJSON,
    GPKG, Shapefile, …) or a :class:`~geopandas.GeoDataFrame`.  When a
    ``.gpkg`` path is given for *road_network*, the ``edges`` layer is read
    and columns are normalised automatically via :func:`load_osm_edges`.

    All attributes from both inputs are carried through to the result,
    prefixed with ``edge_`` and ``zone_`` respectively to avoid collisions.

    The ``proportion`` column indicates what fraction of each edge's
    geometry falls within the intersecting zone. Any columns listed in
    *proportional_cols* are multiplied by this proportion and added as
    new ``proportional_<col>`` columns.

    Parameters
    ----------
    road_network : str, os.PathLike, or gpd.GeoDataFrame
        Road-network edges with line geometries.
    road_network_epsg : int
        EPSG code for the road network's CRS (e.g. ``26910`` for UTM 10N).
    zones : str, os.PathLike, or gpd.GeoDataFrame
        Zone polygons.
    proportional_cols : str or list of str, optional
        Edge attribute column(s) to scale by the intersection proportion.
        For each column ``col``, a ``proportional_<col>`` column is added
        to the result (e.g. ``proportional_cols=["edge_length", "vmt"]``
        produces ``proportional_edge_length`` and ``proportional_vmt``).
        Special values are also supported:
        - ``"edge_length_m"``: include full edge length in meters.
        - ``"proportional_length_m"``: include intersected segment length in meters.
    output_path : str, optional
        If provided, save the result to this file (GeoJSON, GPKG, etc.).
    output_epsg : int, optional
        EPSG code for the output CRS.  Defaults to *road_network_epsg*
        if not provided.

    Returns
    -------
    gpd.GeoDataFrame
        One row per edge-zone intersection piece.  Columns include all
        original edge attributes (prefixed ``edge_``), all original zone
        attributes (prefixed ``zone_``), plus ``proportion``,
        optionally ``edge_length_m`` / ``proportional_length_m`` when
        requested, any requested ``proportional_*`` columns, and ``geometry``.
    """
    edges_gdf = _load_edges(road_network)
    polygons_gdf = _load_zones(zones)

    if output_epsg is None:
        output_epsg = road_network_epsg

    if proportional_cols is None:
        proportional_cols = []
    elif isinstance(proportional_cols, str):
        proportional_cols = [proportional_cols]

    include_edge_length_m = "edge_length_m" in proportional_cols
    include_proportional_length_m = "proportional_length_m" in proportional_cols
    proportional_cols = [
        c for c in proportional_cols if c not in {"edge_length_m", "proportional_length_m"}
    ]

    # Validate that requested columns exist
    for col in proportional_cols:
        if col not in edges_gdf.columns:
            raise ValueError(
                f"proportional_cols: '{col}' not found in road network columns "
                f"{list(edges_gdf.columns)}"
            )

    # Project both inputs to road_network_epsg for intersection
    logger.info("Projecting geometries to EPSG:%d", road_network_epsg)
    edges_proj = edges_gdf.to_crs(epsg=road_network_epsg)
    polys_proj = polygons_gdf.to_crs(epsg=road_network_epsg)

    sindex = edges_proj.sindex
    edge_attr_cols = [c for c in edges_gdf.columns if c != "geometry"]
    zone_attr_cols = [c for c in polygons_gdf.columns if c != "geometry"]

    results = []
    processed_zones = 0
    hit_zones = 0
    intersected_edges = 0
    logger.info("Intersecting %d zones with %d edges", len(polys_proj), len(edges_proj))

    with tqdm(
        total=len(polys_proj),
        desc="Intersecting zones",
        unit="zone",
        dynamic_ncols=True,
    ) as pbar:
        for _, poly_row in polys_proj.iterrows():
            processed_zones += 1
            poly_geom = poly_row.geometry

            candidates_idx = list(sindex.intersection(poly_geom.bounds))
            if not candidates_idx:
                pbar.update(1)
                if processed_zones == 1 or processed_zones % 25 == 0:
                    pbar.set_postfix(
                        hit_zones=hit_zones,
                        pieces=len(results),
                        edge_hits=intersected_edges,
                    )
                continue

            candidates = edges_proj.iloc[candidates_idx]
            intersecting = candidates[candidates.geometry.intersects(poly_geom)]
            if intersecting.empty:
                pbar.update(1)
                if processed_zones == 1 or processed_zones % 25 == 0:
                    pbar.set_postfix(
                        hit_zones=hit_zones,
                        pieces=len(results),
                        edge_hits=intersected_edges,
                    )
                continue

            hit_zones += 1

            for _, edge_row in intersecting.iterrows():
                intersected_edges += 1
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
                    "proportion": proportion,
                    "geometry": intersection_geom,
                }
                if include_edge_length_m:
                    record["edge_length_m"] = round(edge_length, 2)
                if include_proportional_length_m:
                    record["proportional_length_m"] = round(intersection_length, 2)

                for col in edge_attr_cols:
                    record[f"edge_{col}"] = edge_row[col]

                for col in proportional_cols:
                    record[f"proportional_{col}"] = edge_row[col] * proportion

                for col in zone_attr_cols:
                    record[f"zone_{col}"] = poly_row[col]

                results.append(record)

            pbar.update(1)
            if processed_zones == 1 or processed_zones % 25 == 0:
                pbar.set_postfix(
                    hit_zones=hit_zones,
                    pieces=len(results),
                    edge_hits=intersected_edges,
                )

        pbar.set_postfix(
            hit_zones=hit_zones,
            pieces=len(results),
            edge_hits=intersected_edges,
        )

    logger.info(
        "Intersection produced %d results across %d/%d zones (edge hits=%d)",
        len(results),
        hit_zones,
        len(polys_proj),
        intersected_edges,
    )

    if not results:
        logger.warning("No intersections found — returning empty GeoDataFrame")
        # Column order must match the record dict built in the non-empty path:
        # proportion, geometry, [edge_length_m], [proportional_length_m], edge_*, proportional_*, zone_*
        empty_cols = ["proportion", "geometry"]
        if include_edge_length_m:
            empty_cols.append("edge_length_m")
        if include_proportional_length_m:
            empty_cols.append("proportional_length_m")
        return gpd.GeoDataFrame(
            columns=empty_cols,
            geometry="geometry",
            crs=output_epsg,
        )

    result_gdf = gpd.GeoDataFrame(results, geometry="geometry", crs=road_network_epsg)
    if output_epsg != road_network_epsg:
        result_gdf = result_gdf.to_crs(epsg=output_epsg)

    if output_path:
        from osm_chordify.utils.io import save_geodataframe
        save_geodataframe(result_gdf, output_path)
        logger.info("Saved intersection to %s", output_path)

    return result_gdf
