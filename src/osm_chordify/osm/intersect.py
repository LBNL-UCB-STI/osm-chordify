"""Road-network / zone-polygon intersection.

* :func:`load_osm_edges` — loads edges from an OSM GPKG file.
* :func:`intersect_road_network_with_zones` — intersects road-network edges
  with zone polygons (TAZ, census tracts, ISRM cells, …) and returns
  proportional lengths per zone.
"""

import logging
import os
import json
import warnings

import geopandas as gpd
import pandas as pd
from pyproj import CRS
from tqdm import tqdm

logger = logging.getLogger(__name__)

LINE_GEOMETRY_TYPES = {"LineString", "MultiLineString"}
INTERSECTION_CHUNK_SIZE = 50000


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
    if path.endswith(".parquet"):
        return gpd.read_parquet(path)
    return gpd.read_file(path)


def _load_zones(zones):
    """Resolve *zones* to a GeoDataFrame of polygons."""
    if isinstance(zones, gpd.GeoDataFrame):
        return zones
    if not isinstance(zones, (str, os.PathLike)):
        raise TypeError(
            f"zones must be a GeoDataFrame or a file path, got {type(zones).__name__}"
        )
    path = str(zones)
    if path.endswith(".parquet"):
        return gpd.read_parquet(path)
    return gpd.read_file(path)


def _project_if_needed(gdf, epsg):
    """Project a GeoDataFrame only when its CRS differs from the target EPSG."""
    if gdf.crs is not None:
        current_epsg = gdf.crs.to_epsg()
        if current_epsg == epsg:
            return gdf
    return gdf.to_crs(epsg=epsg)


def _is_geographic_crs(crs):
    """Return True when the CRS is geographic (degrees), else False."""
    return crs is not None and getattr(crs, "is_geographic", False)


def _require_projected_epsg(epsg):
    """Validate that the working intersection CRS is projected."""
    crs = CRS.from_user_input(epsg)
    if crs.is_geographic:
        raise ValueError(
            f"Intersection EPSG:{epsg} is geographic. "
            "Use a projected CRS for length-based intersection outputs."
        )
    return crs


def _load_saved_geodataframe(path):
    """Load a previously saved GeoDataFrame from a supported output format."""
    ext = os.path.splitext(str(path))[1].lower()
    if ext == ".parquet":
        return gpd.read_parquet(path)
    if ext == ".gpkg":
        return gpd.read_file(path)
    return gpd.read_file(path)


def _cache_sidecar_path(output_path):
    """Return the metadata sidecar path for an intersection output."""
    return f"{output_path}.cache.json"


def _fingerprint_source(source):
    """Build a lightweight fingerprint for a path-based intersection input."""
    if isinstance(source, (str, os.PathLike)):
        path = os.path.abspath(os.path.expanduser(str(source)))
        if os.path.exists(path):
            stat = os.stat(path)
            return {
                "kind": "path",
                "path": path,
                "size": stat.st_size,
                "mtime_ns": stat.st_mtime_ns,
            }
        return {
            "kind": "path",
            "path": path,
            "missing": True,
        }

    if isinstance(source, gpd.GeoDataFrame):
        bounds = [round(v, 6) for v in source.total_bounds] if len(source) else []
        return {
            "kind": "gdf",
            "rows": int(len(source)),
            "crs": source.crs.to_string() if source.crs is not None else None,
            "bounds": bounds,
            "columns": list(source.columns),
        }

    return {"kind": type(source).__name__}


def _intersection_cache_metadata(
    road_network,
    road_network_epsg,
    zones,
    output_epsg,
):
    """Return metadata used to validate a reusable intersection cache."""
    return {
        "road_network": _fingerprint_source(road_network),
        "road_network_epsg": road_network_epsg,
        "zones": _fingerprint_source(zones),
        "output_epsg": output_epsg,
        "schema_version": 1,
    }


def _try_load_cached_intersection(output_path, metadata):
    """Load an existing cached intersection when the metadata matches."""
    if not output_path or not os.path.exists(output_path):
        return None

    cache_path = _cache_sidecar_path(output_path)
    if not os.path.exists(cache_path):
        return None

    try:
        with open(cache_path, "r", encoding="utf-8") as f:
            cached = json.load(f)
    except Exception:
        return None

    if cached != metadata:
        return None

    logger.info("Reusing cached intersection from %s", output_path)
    return _load_saved_geodataframe(output_path)


def _write_intersection_cache_metadata(output_path, metadata):
    """Persist intersection cache metadata next to the saved output."""
    if not output_path:
        return
    cache_path = _cache_sidecar_path(output_path)
    with open(cache_path, "w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2, sort_keys=True)


def _edge_output_name(col):
    """Preserve already-prefixed chained columns instead of stacking edge_ repeatedly."""
    if col.startswith(("edge_", "zone_")) or col in {
        "zone_edge_proportion",
        "edge_link_length_m",
        "zone_link_length_m",
    }:
        return col
    return f"edge_{col}"


def _zone_output_name(col, existing_keys):
    """Prefix zone attributes and avoid collisions with carried-through columns."""
    candidate = f"zone_{col}"
    if candidate not in existing_keys:
        return candidate
    candidate = f"zone2_{col}"
    if candidate not in existing_keys:
        return candidate
    suffix = 3
    while f"zone{suffix}_{col}" in existing_keys:
        suffix += 1
    return f"zone{suffix}_{col}"


def _build_result_gdf(matched, edge_attr_cols, zone_attr_cols, crs):
    """Build the standardized intersection GeoDataFrame from matched rows."""
    result_cols = ["zone_edge_proportion", "edge_link_length_m", "zone_link_length_m", "geometry"]
    result_data = {
        "zone_edge_proportion": matched["zone_edge_proportion"],
        "edge_link_length_m": matched["edge_link_length_m"],
        "zone_link_length_m": matched["zone_link_length_m"],
        "geometry": matched["geometry"],
    }

    existing_keys = set(result_cols)
    for col in edge_attr_cols:
        out_col = _edge_output_name(col)
        result_data[out_col] = matched[col]
        if out_col not in result_cols:
            result_cols.append(out_col)
            existing_keys.add(out_col)

    for col in zone_attr_cols:
        out_col = _zone_output_name(col, existing_keys)
        result_data[out_col] = matched[col]
        if out_col not in result_cols:
            result_cols.append(out_col)
            existing_keys.add(out_col)

    return gpd.GeoDataFrame(result_data, geometry="geometry", crs=crs)[result_cols]


def _compute_intersection_chunk(chunk, crs):
    """Compute exact intersections for one matched edge/zone chunk."""
    edge_geom_series = gpd.GeoSeries(chunk["edge_geometry"], crs=crs)
    zone_geom_series = gpd.GeoSeries(chunk["zone_geometry"], crs=crs)
    intersection_series = gpd.GeoSeries(
        edge_geom_series.intersection(zone_geom_series),
        crs=crs,
    )
    nonempty_mask = ~intersection_series.is_empty
    line_mask = intersection_series.geom_type.isin(LINE_GEOMETRY_TYPES)
    chunk = chunk.loc[nonempty_mask & line_mask].copy()
    if chunk.empty:
        return chunk

    intersection_series = intersection_series.loc[nonempty_mask & line_mask]
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="Geometry is in a geographic CRS.*",
            category=UserWarning,
        )
        edge_lengths = edge_geom_series.loc[chunk.index].length.round(2)
        zone_lengths = intersection_series.length.round(2)

    chunk["geometry"] = intersection_series
    chunk["edge_link_length_m"] = edge_lengths
    chunk["zone_link_length_m"] = zone_lengths
    chunk["zone_edge_proportion"] = (
        (zone_lengths / edge_lengths).where(edge_lengths > 0, 0).round(4)
    )
    return chunk


def _is_county_like_zones(polys_proj):
    """Heuristic to enable the fast path for small county-style polygon sets."""
    return len(polys_proj) <= 100 and (
        "COUNTYFP" in polys_proj.columns or "GEOID" in polys_proj.columns
    )


def intersect_road_network_with_zones(
    road_network,
    road_network_epsg,
    zones,
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

    The ``zone_edge_proportion`` column indicates what fraction of each edge's
    geometry falls within the intersecting zone. Two length columns are
    always included:
    - ``edge_link_length_m``: full original edge length in meters.
    - ``zone_link_length_m``: intersected segment length inside the zone.
    Parameters
    ----------
    road_network : str, os.PathLike, or gpd.GeoDataFrame
        Road-network edges with line geometries.
    road_network_epsg : int
        EPSG code for the road network's CRS (e.g. ``26910`` for UTM 10N).
    zones : str, os.PathLike, or gpd.GeoDataFrame
        Zone polygons.
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
        attributes (prefixed ``zone_``), plus ``zone_edge_proportion``,
        ``edge_link_length_m``, ``zone_link_length_m``, and ``geometry``.
    """
    _require_projected_epsg(road_network_epsg)
    edges_gdf = _load_edges(road_network)
    polygons_gdf = _load_zones(zones)

    if output_epsg is None:
        output_epsg = road_network_epsg

    cache_metadata = _intersection_cache_metadata(
        road_network=road_network,
        road_network_epsg=road_network_epsg,
        zones=zones,
        output_epsg=output_epsg,
    )
    cached_result = _try_load_cached_intersection(output_path, cache_metadata)
    if cached_result is not None:
        return cached_result

    # Project both inputs to road_network_epsg for intersection
    logger.info("Projecting geometries to EPSG:%d", road_network_epsg)
    edges_proj = _project_if_needed(edges_gdf, road_network_epsg)
    polys_proj = _project_if_needed(polygons_gdf, road_network_epsg)

    edge_attr_cols = [c for c in edges_gdf.columns if c != "geometry"]
    zone_attr_cols = [c for c in polygons_gdf.columns if c != "geometry"]

    logger.info("Intersecting %d zones with %d edges", len(polys_proj), len(edges_proj))

    edges_indexed = edges_proj.reset_index(drop=True).reset_index(names="__edge_idx")
    polys_indexed = polys_proj.reset_index(drop=True).reset_index(names="__zone_idx")
    direct_matches = None
    crossing_edges = edges_indexed

    if _is_county_like_zones(polys_indexed):
        logger.info("Using county fast path for non-boundary-crossing edges")
        boundaries = polys_indexed.geometry.boundary
        boundary_union = boundaries.union_all() if hasattr(boundaries, "union_all") else boundaries.unary_union
        crosses_boundary = edges_indexed.geometry.intersects(boundary_union)
        contained_edges = edges_indexed.loc[~crosses_boundary].copy()
        crossing_edges = edges_indexed.loc[crosses_boundary].copy()

        if not contained_edges.empty:
            reps = contained_edges.copy()
            reps["geometry"] = reps.geometry.representative_point()
            direct_matches = gpd.sjoin(
                reps[["__edge_idx", "geometry"]],
                polys_indexed[["__zone_idx", "geometry"]],
                how="inner",
                predicate="within",
            ).drop(columns=["index_right"])
            direct_matches = direct_matches.merge(
                edges_indexed.drop(columns=["geometry"]),
                on="__edge_idx",
                how="left",
            ).merge(
                polys_indexed.rename(columns={"geometry": "zone_geometry"}),
                on="__zone_idx",
                how="left",
            )
            direct_matches["geometry"] = contained_edges.set_index("__edge_idx").loc[
                direct_matches["__edge_idx"], "geometry"
            ].values
            edge_lengths = gpd.GeoSeries(direct_matches["geometry"], crs=road_network_epsg).length.round(2)
            direct_matches["edge_link_length_m"] = edge_lengths
            direct_matches["zone_link_length_m"] = edge_lengths
            direct_matches["zone_edge_proportion"] = 1.0

    pairs = gpd.sjoin(
        crossing_edges[["__edge_idx", "geometry"]],
        polys_indexed[["__zone_idx", "geometry"]],
        how="inner",
        predicate="intersects",
    )
    pairs = pairs.rename(columns={"geometry": "edge_geometry"}).drop(columns=["index_right"])

    edge_hits = len(pairs) + (0 if direct_matches is None else len(direct_matches))
    hit_zones = (
        pairs["__zone_idx"].nunique() if len(pairs) else 0
    ) + (
        0 if direct_matches is None else direct_matches["__zone_idx"].nunique()
    )

    logger.info(
        "Spatial join produced %d candidate edge-zone pairs across %d/%d zones",
        edge_hits,
        hit_zones,
        len(polys_proj),
    )

    if len(pairs) == 0 and (direct_matches is None or direct_matches.empty):
        logger.warning("No intersections found — returning empty GeoDataFrame")
        empty_cols = ["zone_edge_proportion", "edge_link_length_m", "zone_link_length_m", "geometry"]
        return gpd.GeoDataFrame(
            columns=empty_cols,
            geometry="geometry",
            crs=output_epsg,
        )

    exact_parts = []
    if len(pairs):
        matched = pairs.merge(
            edges_indexed.drop(columns=["geometry"]),
            on="__edge_idx",
            how="left",
        ).merge(
            polys_indexed.rename(columns={"geometry": "zone_geometry"}),
            on="__zone_idx",
            how="left",
        )

        with tqdm(
            total=len(matched),
            desc="Computing intersections",
            unit="pair",
            dynamic_ncols=True,
            leave=True,
            mininterval=0.5,
        ) as pbar:
            for start in range(0, len(matched), INTERSECTION_CHUNK_SIZE):
                chunk = matched.iloc[start:start + INTERSECTION_CHUNK_SIZE].copy()
                chunk_result = _compute_intersection_chunk(chunk, road_network_epsg)
                if not chunk_result.empty:
                    exact_parts.append(chunk_result)
                pbar.update(len(chunk))
                pbar.set_postfix_str(
                    f"hit_zones={hit_zones}, pieces={sum(len(p) for p in exact_parts)}, edge_hits={edge_hits}"
                )

    frames = []
    if direct_matches is not None and not direct_matches.empty:
        frames.append(direct_matches)
    if exact_parts:
        frames.extend(exact_parts)

    if not frames:
        logger.warning("No line intersections found — returning empty GeoDataFrame")
        empty_cols = ["zone_edge_proportion", "edge_link_length_m", "zone_link_length_m", "geometry"]
        return gpd.GeoDataFrame(
            columns=empty_cols,
            geometry="geometry",
            crs=output_epsg,
        )

    matched_all = gpd.GeoDataFrame(
        frames[0] if len(frames) == 1 else pd.concat(frames, ignore_index=True),
        geometry="geometry",
        crs=road_network_epsg,
    )
    result_gdf = _build_result_gdf(matched_all, edge_attr_cols, zone_attr_cols, road_network_epsg)

    logger.info(
        "Intersection produced %d results across %d/%d zones (edge hits=%d)",
        len(result_gdf),
        hit_zones,
        len(polys_proj),
        edge_hits,
    )

    if output_epsg != road_network_epsg:
        result_gdf = result_gdf.to_crs(epsg=output_epsg)

    if output_path:
        from osm_chordify.utils.io import save_geodataframe
        save_geodataframe(result_gdf, output_path)
        _write_intersection_cache_metadata(output_path, cache_metadata)
        logger.info("Saved intersection to %s", output_path)

    return result_gdf
