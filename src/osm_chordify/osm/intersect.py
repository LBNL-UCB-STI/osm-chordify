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
from concurrent.futures import ThreadPoolExecutor, as_completed

import geopandas as gpd
import pandas as pd
from pyproj import CRS
from shapely.geometry import GeometryCollection, box
from tqdm import tqdm

logger = logging.getLogger(__name__)

LINE_GEOMETRY_TYPES = {"LineString", "MultiLineString"}
POLYGON_GEOMETRY_TYPES = {"Polygon", "MultiPolygon"}
INTERSECTION_CHUNK_SIZE = 50000
COUNTY_BOUNDARY_SCREEN_CHUNK_SIZE = 100000
PARALLEL_INTERSECTION_MIN_PAIRS = 200000
MAX_INTERSECTION_WORKERS = min(8, os.cpu_count() or 1)


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
    prefilter_zones_to_network_bbox,
    zone_label=None,
):
    """Return metadata used to validate a reusable intersection cache."""
    return {
        "road_network": _fingerprint_source(road_network),
        "road_network_epsg": road_network_epsg,
        "zones": _fingerprint_source(zones),
        "output_epsg": output_epsg,
        "prefilter_zones_to_network_bbox": prefilter_zones_to_network_bbox,
        "zone_label": _label_prefix(zone_label),
        "schema_version": 3,
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
        "edge_surface_m2",
        "zone_surface_m2",
        "zone_piece_proportion",
        "piece_link_length_m",
        "zone_piece_length_m",
        "piece_surface_m2",
    }:
        return col
    return f"edge_{col}"


def _zone_output_name(col, existing_keys, prefix="zone"):
    """Prefix zone attributes and avoid collisions with carried-through columns."""
    candidate = f"{prefix}_{col}"
    if candidate not in existing_keys:
        return candidate
    candidate = f"{prefix}2_{col}"
    if candidate not in existing_keys:
        return candidate
    suffix = 3
    while f"{prefix}{suffix}_{col}" in existing_keys:
        suffix += 1
    return f"{prefix}{suffix}_{col}"


def _label_prefix(zone_label):
    return zone_label.strip().lower() if zone_label else None


def _line_metric_names(zone_label=None):
    prefix = _label_prefix(zone_label)
    if prefix is None:
        return {
            "proportion": "zone_edge_proportion",
            "edge_length": "edge_link_length_m",
            "zone_length": "zone_link_length_m",
        }
    return {
        "proportion": f"{prefix}_zone_edge_proportion",
        "edge_length": f"{prefix}_edge_link_length_m",
        "zone_length": f"{prefix}_zone_link_length_m",
    }


def _polygon_metric_names(zone_label=None):
    prefix = _label_prefix(zone_label)
    if prefix is None:
        return {
            "proportion": "zone_edge_proportion",
            "edge_length": "edge_link_length_m",
            "zone_length": "zone_link_length_m",
            "edge_surface": "edge_surface_m2",
            "zone_surface": "zone_surface_m2",
        }
    return {
        "proportion": f"{prefix}_zone_edge_proportion",
        "edge_length": f"{prefix}_edge_link_length_m",
        "zone_length": f"{prefix}_zone_link_length_m",
        "edge_surface": f"{prefix}_edge_surface_m2",
        "zone_surface": f"{prefix}_zone_surface_m2",
    }


def _cascade_metric_names(zone_label=None):
    prefix = _label_prefix(zone_label)
    if prefix is None:
        return {
            "proportion": "zone_piece_proportion",
            "piece_length": "piece_link_length_m",
            "zone_length": "zone_piece_length_m",
            "piece_surface": "piece_surface_m2",
            "zone_surface": "zone_surface_m2",
        }
    return {
        "proportion": f"{prefix}_zone_piece_proportion",
        "piece_length": f"{prefix}_piece_link_length_m",
        "zone_length": f"{prefix}_zone_piece_length_m",
        "piece_surface": f"{prefix}_piece_surface_m2",
        "zone_surface": f"{prefix}_zone_surface_m2",
    }


def _build_result_gdf(matched, edge_attr_cols, zone_attr_cols, crs, metric_names=None, zone_label=None):
    """Build the standardized intersection GeoDataFrame from matched rows."""
    metric_names = metric_names or _line_metric_names()
    result_cols = [
        metric_names["proportion"],
        metric_names["edge_length"],
        metric_names["zone_length"],
        "geometry",
    ]
    result_data = {
        metric_names["proportion"]: matched[metric_names["proportion"]],
        metric_names["edge_length"]: matched[metric_names["edge_length"]],
        metric_names["zone_length"]: matched[metric_names["zone_length"]],
        "geometry": matched["geometry"],
    }
    for optional_key in ("edge_surface", "zone_surface"):
        if optional_key in metric_names and metric_names[optional_key] in matched.columns:
            result_data[metric_names[optional_key]] = matched[metric_names[optional_key]]
            result_cols.insert(-1, metric_names[optional_key])

    existing_keys = set(result_cols)
    for col in edge_attr_cols:
        out_col = _edge_output_name(col)
        result_data[out_col] = matched[col]
        if out_col not in result_cols:
            result_cols.append(out_col)
            existing_keys.add(out_col)

    for col in zone_attr_cols:
        out_col = _zone_output_name(col, existing_keys, prefix=_label_prefix(zone_label) or "zone")
        result_data[out_col] = matched[col]
        if out_col not in result_cols:
            result_cols.append(out_col)
            existing_keys.add(out_col)

    return gpd.GeoDataFrame(result_data, geometry="geometry", crs=crs)[result_cols]


def _build_void_rows(zones, edge_attr_cols, zone_attr_cols, crs, metric_names=None, zone_label=None):
    """Build placeholder rows for kept zones with no intersecting link pieces."""
    metric_names = metric_names or _line_metric_names()
    result_cols = [
        metric_names["proportion"],
        metric_names["edge_length"],
        metric_names["zone_length"],
        "geometry",
    ]
    result_data = {
        metric_names["proportion"]: [pd.NA] * len(zones),
        metric_names["edge_length"]: [pd.NA] * len(zones),
        metric_names["zone_length"]: [pd.NA] * len(zones),
        "geometry": [GeometryCollection() for _ in range(len(zones))],
    }
    for optional_key in ("edge_surface", "zone_surface"):
        if optional_key in metric_names:
            result_data[metric_names[optional_key]] = [pd.NA] * len(zones)
            result_cols.insert(-1, metric_names[optional_key])

    existing_keys = set(result_cols)
    for col in edge_attr_cols:
        out_col = _edge_output_name(col)
        result_data[out_col] = [pd.NA] * len(zones)
        if out_col not in result_cols:
            result_cols.append(out_col)
            existing_keys.add(out_col)

    for col in zone_attr_cols:
        out_col = _zone_output_name(col, existing_keys, prefix=_label_prefix(zone_label) or "zone")
        result_data[out_col] = zones[col].tolist()
        if out_col not in result_cols:
            result_cols.append(out_col)
            existing_keys.add(out_col)

    return gpd.GeoDataFrame(result_data, geometry="geometry", crs=crs)[result_cols]


def _compute_intersection_chunk(chunk, crs, metric_names=None):
    """Compute exact intersections for one matched edge/zone chunk."""
    metric_names = metric_names or _line_metric_names()
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
    chunk[metric_names["edge_length"]] = edge_lengths
    chunk[metric_names["zone_length"]] = zone_lengths
    chunk[metric_names["proportion"]] = (
        (zone_lengths / edge_lengths).where(edge_lengths > 0, 0).round(4)
    )
    return chunk


def _is_county_like_zones(polys_proj):
    """Heuristic to enable the fast path for small county-style polygon sets."""
    return len(polys_proj) <= 100 and (
        "COUNTYFP" in polys_proj.columns or "GEOID" in polys_proj.columns
    )


def _prefilter_zones_to_network_bbox(polys_proj, edges_proj):
    """Keep only zones that intersect the overall road-network bounding box."""
    if edges_proj.empty or polys_proj.empty:
        return polys_proj.iloc[0:0].copy()

    minx, miny, maxx, maxy = edges_proj.total_bounds
    network_bbox = box(minx, miny, maxx, maxy)
    logger.info(
        "Filtering %d zones against the road-network bounding box before exact intersection",
        len(polys_proj),
    )
    with tqdm(
        total=2,
        desc="Filtering zones",
        unit="step",
        dynamic_ncols=True,
        leave=True,
        mininterval=0.5,
    ) as pbar:
        pbar.set_postfix_str("building network bbox")
        search_bounds = network_bbox.bounds
        pbar.update(1)

        pbar.set_postfix_str("screening zones")
        candidate_zone_indices = polys_proj.sindex.intersection(search_bounds)
        if len(candidate_zone_indices):
            candidate_grid = polys_proj.iloc[sorted(candidate_zone_indices)].copy()
            exact_keep = candidate_grid.geometry.intersects(network_bbox)
            filtered = candidate_grid.loc[exact_keep].copy()
        else:
            filtered = polys_proj.iloc[0:0].copy()
        pbar.update(1)
        pbar.set_postfix_str(f"kept={len(filtered)}, dropped={len(polys_proj) - len(filtered)}")

    logger.info(
        "Zone prefilter kept %d/%d zones after network-bbox screening",
        len(filtered),
        len(polys_proj),
    )
    return filtered


def _chunk_slices(total, chunk_size):
    """Yield slice boundaries for chunked processing."""
    for start in range(0, total, chunk_size):
        yield start, min(start + chunk_size, total)


def _screen_crosses_boundary_chunk(geoms, boundary_union):
    """Return boolean intersects(boundary_union) for one geometry chunk."""
    return geoms.intersects(boundary_union)


def _compute_chunk_worker(func, chunk, epsg, metric_names, length_col):
    """Top-level worker wrapper for parallel chunk execution."""
    if length_col is None:
        return func(chunk, epsg, metric_names=metric_names)
    return func(chunk, epsg, length_col, metric_names=metric_names)


def _compute_chunked_exact_parts(
    matched,
    epsg,
    compute_func,
    desc,
    postfix_template,
    metric_names=None,
    length_col=None,
    enable_parallel=False,
):
    """Run chunked exact geometry processing, optionally in parallel."""
    exact_parts = []
    if not len(matched):
        return exact_parts

    total = len(matched)
    chunk_specs = [(start, matched.iloc[start:end].copy()) for start, end in _chunk_slices(total, INTERSECTION_CHUNK_SIZE)]
    use_parallel = enable_parallel and MAX_INTERSECTION_WORKERS > 1 and total >= PARALLEL_INTERSECTION_MIN_PAIRS

    with tqdm(
        total=total,
        desc=desc,
        unit="pair",
        dynamic_ncols=True,
        leave=True,
        mininterval=0.5,
    ) as pbar:
        if use_parallel:
            with ThreadPoolExecutor(max_workers=MAX_INTERSECTION_WORKERS) as executor:
                futures = {
                    executor.submit(
                        _compute_chunk_worker,
                        compute_func,
                        chunk,
                        epsg,
                        metric_names,
                        length_col,
                    ): len(chunk)
                    for _, chunk in chunk_specs
                }
                for future in as_completed(futures):
                    chunk_result = future.result()
                    if not chunk_result.empty:
                        exact_parts.append(chunk_result)
                    pbar.update(futures[future])
                    pbar.set_postfix_str(postfix_template(sum(len(p) for p in exact_parts)))
        else:
            for _, chunk in chunk_specs:
                chunk_result = _compute_chunk_worker(
                    compute_func,
                    chunk,
                    epsg,
                    metric_names,
                    length_col,
                )
                if not chunk_result.empty:
                    exact_parts.append(chunk_result)
                pbar.update(len(chunk))
                pbar.set_postfix_str(postfix_template(sum(len(p) for p in exact_parts)))

    return exact_parts


def intersect_road_network_with_zones(
    road_network,
    road_network_epsg,
    zones,
    output_path=None,
    output_epsg=None,
    prefilter_zones_to_network_bbox=False,
    zone_label=None,
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
    prefilter_zones_to_network_bbox : bool, optional
        If ``True``, keep only zones that intersect the overall road-network
        bounding box before exact intersection.
    zone_label : str, optional
        Prefix used for the current-step zone attributes and metric columns.
        For example, ``zone_label="aermod"`` produces columns such as
        ``aermod_zone_edge_proportion`` and ``aermod_COUNTYFP``.
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
        prefilter_zones_to_network_bbox=prefilter_zones_to_network_bbox,
        zone_label=zone_label,
    )
    cached_result = _try_load_cached_intersection(output_path, cache_metadata)
    if cached_result is not None:
        return cached_result

    # Project both inputs to road_network_epsg for intersection
    logger.info("Projecting geometries to EPSG:%d", road_network_epsg)
    edges_proj = _project_if_needed(edges_gdf, road_network_epsg)
    polys_proj = _project_if_needed(polygons_gdf, road_network_epsg)

    if prefilter_zones_to_network_bbox and len(polys_proj):
        polys_proj = _prefilter_zones_to_network_bbox(polys_proj, edges_proj)

    edge_attr_cols = [c for c in edges_gdf.columns if c != "geometry"]
    zone_attr_cols = [c for c in polygons_gdf.columns if c != "geometry"]
    metric_names = _line_metric_names(zone_label)

    logger.info("Intersecting %d zones with %d edges", len(polys_proj), len(edges_proj))

    edges_indexed = edges_proj.reset_index(drop=True).reset_index(names="__edge_idx")
    polys_indexed = polys_proj.reset_index(drop=True).reset_index(names="__zone_idx")
    direct_matches = None
    crossing_edges = edges_indexed

    if _is_county_like_zones(polys_indexed):
        logger.info("Using county fast path for non-boundary-crossing edges")
        boundaries = polys_indexed.geometry.boundary
        boundary_union = boundaries.union_all() if hasattr(boundaries, "union_all") else boundaries.unary_union
        cross_masks = []
        for start, end in _chunk_slices(len(edges_indexed), COUNTY_BOUNDARY_SCREEN_CHUNK_SIZE):
            cross_masks.append(_screen_crosses_boundary_chunk(edges_indexed.geometry.iloc[start:end], boundary_union))
        crosses_boundary = pd.concat(cross_masks, ignore_index=True) if cross_masks else pd.Series(dtype=bool)
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
            direct_matches[metric_names["edge_length"]] = edge_lengths
            direct_matches[metric_names["zone_length"]] = edge_lengths
            direct_matches[metric_names["proportion"]] = 1.0

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
        if prefilter_zones_to_network_bbox and len(polys_proj):
            logger.warning("No intersections found inside bbox-filtered zones — returning voided zone rows")
            voided = _build_void_rows(
                polys_proj,
                edge_attr_cols,
                zone_attr_cols,
                output_epsg,
                metric_names=metric_names,
                zone_label=zone_label,
            )
            if output_path:
                from osm_chordify.utils.io import save_geodataframe
                save_geodataframe(voided, output_path)
                _write_intersection_cache_metadata(output_path, cache_metadata)
                logger.info("Saved intersection to %s", output_path)
            return voided
        logger.warning("No intersections found — returning empty GeoDataFrame")
        empty_cols = [metric_names["proportion"], metric_names["edge_length"], metric_names["zone_length"], "geometry"]
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
        exact_parts = _compute_chunked_exact_parts(
            matched=matched,
            epsg=road_network_epsg,
            compute_func=_compute_intersection_chunk,
            desc="Computing intersections",
            postfix_template=lambda pieces: f"hit_zones={hit_zones}, pieces={pieces}, edge_hits={edge_hits}",
            metric_names=metric_names,
            enable_parallel=_is_county_like_zones(polys_indexed),
        )

    frames = []
    if direct_matches is not None and not direct_matches.empty:
        frames.append(direct_matches)
    if exact_parts:
        frames.extend(exact_parts)

    if not frames:
        logger.warning("No line intersections found — returning empty GeoDataFrame")
        empty_cols = [metric_names["proportion"], metric_names["edge_length"], metric_names["zone_length"], "geometry"]
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
    result_gdf = _build_result_gdf(
        matched_all,
        edge_attr_cols,
        zone_attr_cols,
        road_network_epsg,
        metric_names=metric_names,
        zone_label=zone_label,
    )
    if prefilter_zones_to_network_bbox and "__zone_idx" in matched_all.columns:
        matched_zone_indices = set(matched_all["__zone_idx"].unique())
        missing_zone_rows = polys_indexed.loc[~polys_indexed["__zone_idx"].isin(matched_zone_indices)]
        if len(missing_zone_rows):
            voided = _build_void_rows(
                missing_zone_rows,
                edge_attr_cols,
                zone_attr_cols,
                road_network_epsg,
                metric_names=metric_names,
                zone_label=zone_label,
            )
            result_gdf = pd.concat([result_gdf, voided], ignore_index=True)
            result_gdf = gpd.GeoDataFrame(result_gdf, geometry="geometry", crs=road_network_epsg)

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


def _infer_edge_length_col(gdf, edge_length_col=None):
    """Resolve the source column used for full-link length values."""
    if edge_length_col is not None:
        if edge_length_col not in gdf.columns:
            raise ValueError(f"road_network is missing requested edge length column '{edge_length_col}'")
        return edge_length_col

    for candidate in ("edge_link_length_m", "edge_length_m", "edge_length", "length"):
        if candidate in gdf.columns:
            return candidate

    raise ValueError(
        "road_network polygon inputs must include a length column. "
        "Pass edge_length_col explicitly or provide one of: "
        "'edge_link_length_m', 'edge_length_m', 'edge_length', 'length'."
    )


def _infer_polygon_piece_length_col(gdf, piece_length_col=None):
    """Resolve the source column used for cascade polygon-piece lengths."""
    if piece_length_col is not None:
        if piece_length_col not in gdf.columns:
            raise ValueError(f"polygons are missing requested piece length column '{piece_length_col}'")
        return piece_length_col

    for candidate in ("zone_link_length_m", "edge_link_length_m", "edge_length_m", "edge_length", "length"):
        if candidate in gdf.columns:
            return candidate
    for suffix in ("_zone_piece_length_m", "_piece_link_length_m", "_zone_link_length_m", "_edge_link_length_m"):
        matches = [col for col in gdf.columns if col.endswith(suffix)]
        if matches:
            return matches[0]

    raise ValueError(
        "polygon cascade inputs must include a length column. "
        "Pass piece_length_col explicitly or provide one of: "
        "'zone_link_length_m', 'edge_link_length_m', 'edge_length_m', 'edge_length', 'length'."
    )


def _compute_polygon_intersection_chunk(chunk, crs, edge_length_col, metric_names=None):
    """Compute polygon-road/zone overlap pieces using area-based proportions."""
    metric_names = metric_names or _polygon_metric_names()
    edge_geom_series = gpd.GeoSeries(chunk["edge_geometry"], crs=crs)
    zone_geom_series = gpd.GeoSeries(chunk["zone_geometry"], crs=crs)
    intersection_series = gpd.GeoSeries(edge_geom_series.intersection(zone_geom_series), crs=crs)
    nonempty_mask = ~intersection_series.is_empty
    polygon_mask = intersection_series.geom_type.isin(POLYGON_GEOMETRY_TYPES)
    chunk = chunk.loc[nonempty_mask & polygon_mask].copy()
    if chunk.empty:
        return chunk

    intersection_series = intersection_series.loc[nonempty_mask & polygon_mask]
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="Geometry is in a geographic CRS.*",
            category=UserWarning,
        )
        edge_areas = edge_geom_series.loc[chunk.index].area.round(2)
        overlap_areas = intersection_series.area.round(2)

    proportions = (overlap_areas / edge_areas).where(edge_areas > 0, 0).round(4)
    full_lengths = pd.to_numeric(chunk[edge_length_col], errors="coerce")

    chunk["geometry"] = intersection_series
    chunk[metric_names["proportion"]] = proportions
    chunk[metric_names["edge_length"]] = full_lengths.round(2)
    chunk[metric_names["zone_length"]] = (full_lengths * proportions).round(2)
    chunk[metric_names["edge_surface"]] = edge_areas
    chunk[metric_names["zone_surface"]] = overlap_areas
    return chunk


def _build_cascade_result_gdf(matched, polygon_attr_cols, zone_attr_cols, crs, metric_names=None, zone_label=None):
    """Build a polygon-cascade result while preserving existing polygon columns."""
    metric_names = metric_names or _cascade_metric_names()
    result_data = {col: matched[col] for col in polygon_attr_cols}
    result_data[metric_names["proportion"]] = matched[metric_names["proportion"]]
    result_data[metric_names["piece_length"]] = matched[metric_names["piece_length"]]
    result_data[metric_names["zone_length"]] = matched[metric_names["zone_length"]]
    result_data[metric_names["piece_surface"]] = matched[metric_names["piece_surface"]]
    result_data[metric_names["zone_surface"]] = matched[metric_names["zone_surface"]]
    result_data["geometry"] = matched["geometry"]

    result_cols = list(polygon_attr_cols) + [
        metric_names["proportion"],
        metric_names["piece_length"],
        metric_names["zone_length"],
        metric_names["piece_surface"],
        metric_names["zone_surface"],
        "geometry",
    ]
    existing_keys = set(result_cols)

    for col in zone_attr_cols:
        out_col = _zone_output_name(col, existing_keys, prefix=_label_prefix(zone_label) or "zone")
        result_data[out_col] = matched[col]
        if out_col not in result_cols:
            result_cols.append(out_col)
            existing_keys.add(out_col)

    return gpd.GeoDataFrame(result_data, geometry="geometry", crs=crs)[result_cols]


def _build_cascade_void_rows(zones, polygon_attr_cols, zone_attr_cols, crs, metric_names=None, zone_label=None):
    """Build placeholder rows for zones kept by bbox prefilter with no polygon matches."""
    metric_names = metric_names or _cascade_metric_names()
    result_data = {col: [pd.NA] * len(zones) for col in polygon_attr_cols}
    result_data[metric_names["proportion"]] = [pd.NA] * len(zones)
    result_data[metric_names["piece_length"]] = [pd.NA] * len(zones)
    result_data[metric_names["zone_length"]] = [pd.NA] * len(zones)
    result_data[metric_names["piece_surface"]] = [pd.NA] * len(zones)
    result_data[metric_names["zone_surface"]] = [pd.NA] * len(zones)
    result_data["geometry"] = [GeometryCollection() for _ in range(len(zones))]

    result_cols = list(polygon_attr_cols) + [
        metric_names["proportion"],
        metric_names["piece_length"],
        metric_names["zone_length"],
        metric_names["piece_surface"],
        metric_names["zone_surface"],
        "geometry",
    ]
    existing_keys = set(result_cols)

    for col in zone_attr_cols:
        out_col = _zone_output_name(col, existing_keys, prefix=_label_prefix(zone_label) or "zone")
        result_data[out_col] = zones[col].tolist()
        if out_col not in result_cols:
            result_cols.append(out_col)
            existing_keys.add(out_col)

    return gpd.GeoDataFrame(result_data, geometry="geometry", crs=crs)[result_cols]


def _compute_polygon_cascade_chunk(chunk, crs, piece_length_col, metric_names=None):
    """Compute exact polygon/polygon overlaps for cascading zone intersections."""
    metric_names = metric_names or _cascade_metric_names()
    polygon_geom_series = gpd.GeoSeries(chunk["polygon_geometry"], crs=crs)
    zone_geom_series = gpd.GeoSeries(chunk["zone_geometry"], crs=crs)
    intersection_series = gpd.GeoSeries(polygon_geom_series.intersection(zone_geom_series), crs=crs)
    nonempty_mask = ~intersection_series.is_empty
    polygon_mask = intersection_series.geom_type.isin(POLYGON_GEOMETRY_TYPES)
    chunk = chunk.loc[nonempty_mask & polygon_mask].copy()
    if chunk.empty:
        return chunk

    intersection_series = intersection_series.loc[nonempty_mask & polygon_mask]
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="Geometry is in a geographic CRS.*",
            category=UserWarning,
        )
        polygon_areas = polygon_geom_series.loc[chunk.index].area.round(2)
        overlap_areas = intersection_series.area.round(2)

    proportions = (overlap_areas / polygon_areas).where(polygon_areas > 0, 0).round(4)
    piece_lengths = pd.to_numeric(chunk[piece_length_col], errors="coerce").round(2)

    chunk["geometry"] = intersection_series
    chunk[metric_names["proportion"]] = proportions
    chunk[metric_names["piece_length"]] = piece_lengths
    chunk[metric_names["zone_length"]] = (piece_lengths * proportions).round(2)
    chunk[metric_names["piece_surface"]] = polygon_areas
    chunk[metric_names["zone_surface"]] = overlap_areas
    return chunk


def intersect_road_polygons_with_zones(
    road_network,
    road_network_epsg,
    zones,
    output_path=None,
    output_epsg=None,
    prefilter_zones_to_network_bbox=True,
    edge_length_col=None,
    zone_label=None,
):
    """Intersect rectangular/polygon road links with zone polygons.

    This variant is for area-based road links, for example buffered or
    rectangular road segments. ``zone_edge_proportion`` is computed as the
    overlap area divided by the full road-link polygon area. That proportion is
    then applied to the full link length to derive ``zone_link_length_m``.

    By default it first filters zones to the overall road-network bounding box
    before exact intersection. That prefilter uses the same progress-bar path
    as the line-based workflow.
    When ``zone_label`` is provided, current-step metrics and zone attributes
    are prefixed with that label.
    """
    _require_projected_epsg(road_network_epsg)
    edges_gdf = _load_edges(road_network)
    polygons_gdf = _load_zones(zones)

    if output_epsg is None:
        output_epsg = road_network_epsg

    length_col = _infer_edge_length_col(edges_gdf, edge_length_col=edge_length_col)
    cache_metadata = _intersection_cache_metadata(
        road_network=road_network,
        road_network_epsg=road_network_epsg,
        zones=zones,
        output_epsg=output_epsg,
        prefilter_zones_to_network_bbox=prefilter_zones_to_network_bbox,
        zone_label=zone_label,
    )
    cached_result = _try_load_cached_intersection(output_path, {**cache_metadata, "polygon_mode": True, "edge_length_col": length_col})
    if cached_result is not None:
        return cached_result

    logger.info("Projecting geometries to EPSG:%d", road_network_epsg)
    edges_proj = _project_if_needed(edges_gdf, road_network_epsg)
    polys_proj = _project_if_needed(polygons_gdf, road_network_epsg)

    if prefilter_zones_to_network_bbox and len(polys_proj):
        polys_proj = _prefilter_zones_to_network_bbox(polys_proj, edges_proj)

    edge_attr_cols = [c for c in edges_gdf.columns if c != "geometry"]
    zone_attr_cols = [c for c in polygons_gdf.columns if c != "geometry"]
    metric_names = _polygon_metric_names(zone_label)

    logger.info("Intersecting %d zones with %d polygon road links", len(polys_proj), len(edges_proj))
    edges_indexed = edges_proj.reset_index(drop=True).reset_index(names="__edge_idx")
    polys_indexed = polys_proj.reset_index(drop=True).reset_index(names="__zone_idx")

    pairs = gpd.sjoin(
        edges_indexed[["__edge_idx", "geometry"]],
        polys_indexed[["__zone_idx", "geometry"]],
        how="inner",
        predicate="intersects",
    )
    pairs = pairs.rename(columns={"geometry": "edge_geometry"}).drop(columns=["index_right"])

    if len(pairs) == 0:
        if prefilter_zones_to_network_bbox and len(polys_proj):
            voided = _build_void_rows(
                polys_proj,
                edge_attr_cols,
                zone_attr_cols,
                output_epsg,
                metric_names=metric_names,
                zone_label=zone_label,
            )
            if output_path:
                from osm_chordify.utils.io import save_geodataframe
                save_geodataframe(voided, output_path)
                _write_intersection_cache_metadata(output_path, {**cache_metadata, "polygon_mode": True, "edge_length_col": length_col})
            return voided
        return gpd.GeoDataFrame(
            columns=[
                metric_names["proportion"],
                metric_names["edge_length"],
                metric_names["zone_length"],
                metric_names["edge_surface"],
                metric_names["zone_surface"],
                "geometry",
            ],
            geometry="geometry",
            crs=output_epsg,
        )

    matched = pairs.merge(
        edges_indexed.drop(columns=["geometry"]),
        on="__edge_idx",
        how="left",
    ).merge(
        polys_indexed.rename(columns={"geometry": "zone_geometry"}),
        on="__zone_idx",
        how="left",
    )

    exact_parts = []
    with tqdm(
        total=len(matched),
        desc="Computing polygon intersections",
        unit="pair",
        dynamic_ncols=True,
        leave=True,
        mininterval=0.5,
    ) as pbar:
        for start in range(0, len(matched), INTERSECTION_CHUNK_SIZE):
            chunk = matched.iloc[start:start + INTERSECTION_CHUNK_SIZE].copy()
            chunk_result = _compute_polygon_intersection_chunk(
                chunk,
                road_network_epsg,
                length_col,
                metric_names=metric_names,
            )
            if not chunk_result.empty:
                exact_parts.append(chunk_result)
            pbar.update(len(chunk))
            pbar.set_postfix_str(
                f"hit_zones={pairs['__zone_idx'].nunique()}, pieces={sum(len(p) for p in exact_parts)}, edge_hits={len(pairs)}"
            )

    if not exact_parts:
        return gpd.GeoDataFrame(
            columns=[
                metric_names["proportion"],
                metric_names["edge_length"],
                metric_names["zone_length"],
                metric_names["edge_surface"],
                metric_names["zone_surface"],
                "geometry",
            ],
            geometry="geometry",
            crs=output_epsg,
        )

    matched_all = gpd.GeoDataFrame(pd.concat(exact_parts, ignore_index=True), geometry="geometry", crs=road_network_epsg)
    result_gdf = _build_result_gdf(
        matched_all,
        edge_attr_cols,
        zone_attr_cols,
        road_network_epsg,
        metric_names=metric_names,
        zone_label=zone_label,
    )
    if prefilter_zones_to_network_bbox and "__zone_idx" in matched_all.columns:
        matched_zone_indices = set(matched_all["__zone_idx"].unique())
        missing_zone_rows = polys_indexed.loc[~polys_indexed["__zone_idx"].isin(matched_zone_indices)]
        if len(missing_zone_rows):
            voided = _build_void_rows(
                missing_zone_rows,
                edge_attr_cols,
                zone_attr_cols,
                road_network_epsg,
                metric_names=metric_names,
                zone_label=zone_label,
            )
            result_gdf = gpd.GeoDataFrame(pd.concat([result_gdf, voided], ignore_index=True), geometry="geometry", crs=road_network_epsg)

    if output_epsg != road_network_epsg:
        result_gdf = result_gdf.to_crs(epsg=output_epsg)

    if output_path:
        from osm_chordify.utils.io import save_geodataframe
        save_geodataframe(result_gdf, output_path)
        _write_intersection_cache_metadata(output_path, {**cache_metadata, "polygon_mode": True, "edge_length_col": length_col})
        logger.info("Saved polygon-road intersection to %s", output_path)

    return result_gdf


def intersect_polygons_with_zones(
    polygons,
    polygons_epsg,
    zones,
    output_path=None,
    output_epsg=None,
    prefilter_zones_to_network_bbox=True,
    piece_length_col=None,
    zone_label=None,
):
    """Intersect an already-processed polygon layer with a new zone layer.

    Existing polygon-layer columns are preserved as-is. New zone attributes are
    added with prefixed names. The current-step metrics are:
    - ``zone_piece_proportion``: overlap area / full source polygon area
    - ``piece_link_length_m``: source polygon piece length
    - ``zone_piece_length_m``: derived in-zone length based on the proportion
    """
    _require_projected_epsg(polygons_epsg)
    polygons_gdf = _load_zones(polygons)
    zones_gdf = _load_zones(zones)

    if output_epsg is None:
        output_epsg = polygons_epsg

    piece_length_col = _infer_polygon_piece_length_col(polygons_gdf, piece_length_col=piece_length_col)
    cache_metadata = _intersection_cache_metadata(
        road_network=polygons,
        road_network_epsg=polygons_epsg,
        zones=zones,
        output_epsg=output_epsg,
        prefilter_zones_to_network_bbox=prefilter_zones_to_network_bbox,
        zone_label=zone_label,
    )
    cascade_cache_metadata = {**cache_metadata, "cascade_polygon_mode": True, "piece_length_col": piece_length_col}
    cached_result = _try_load_cached_intersection(output_path, cascade_cache_metadata)
    if cached_result is not None:
        return cached_result

    logger.info("Projecting geometries to EPSG:%d", polygons_epsg)
    polygons_proj = _project_if_needed(polygons_gdf, polygons_epsg)
    zones_proj = _project_if_needed(zones_gdf, polygons_epsg)

    if prefilter_zones_to_network_bbox and len(zones_proj):
        zones_proj = _prefilter_zones_to_network_bbox(zones_proj, polygons_proj)

    polygon_attr_cols = [c for c in polygons_gdf.columns if c != "geometry"]
    zone_attr_cols = [c for c in zones_gdf.columns if c != "geometry"]
    metric_names = _cascade_metric_names(zone_label)

    logger.info("Intersecting %d zones with %d polygon pieces", len(zones_proj), len(polygons_proj))
    polygons_indexed = polygons_proj.reset_index(drop=True).reset_index(names="__poly_idx")
    zones_indexed = zones_proj.reset_index(drop=True).reset_index(names="__zone_idx")
    direct_matches = None
    crossing_polygons = polygons_indexed

    if _is_county_like_zones(zones_indexed):
        logger.info("Using county fast path for non-boundary-crossing polygon pieces")
        boundaries = zones_indexed.geometry.boundary
        boundary_union = boundaries.union_all() if hasattr(boundaries, "union_all") else boundaries.unary_union
        cross_masks = []
        for start, end in _chunk_slices(len(polygons_indexed), COUNTY_BOUNDARY_SCREEN_CHUNK_SIZE):
            cross_masks.append(_screen_crosses_boundary_chunk(polygons_indexed.geometry.iloc[start:end], boundary_union))
        crosses_boundary = pd.concat(cross_masks, ignore_index=True) if cross_masks else pd.Series(dtype=bool)
        contained_polygons = polygons_indexed.loc[~crosses_boundary].copy()
        crossing_polygons = polygons_indexed.loc[crosses_boundary].copy()

        if not contained_polygons.empty:
            reps = contained_polygons.copy()
            reps["geometry"] = reps.geometry.representative_point()
            direct_matches = gpd.sjoin(
                reps[["__poly_idx", "geometry"]],
                zones_indexed[["__zone_idx", "geometry"]],
                how="inner",
                predicate="within",
            ).drop(columns=["index_right"])
            direct_matches = direct_matches.merge(
                polygons_indexed.drop(columns=["geometry"]),
                on="__poly_idx",
                how="left",
            ).merge(
                zones_indexed.rename(columns={"geometry": "zone_geometry"}),
                on="__zone_idx",
                how="left",
            )
            direct_matches["geometry"] = contained_polygons.set_index("__poly_idx").loc[
                direct_matches["__poly_idx"], "geometry"
            ].values
            with warnings.catch_warnings():
                warnings.filterwarnings(
                    "ignore",
                    message="Geometry is in a geographic CRS.*",
                    category=UserWarning,
                )
                polygon_areas = gpd.GeoSeries(direct_matches["geometry"], crs=polygons_epsg).area.round(2)
            piece_lengths = pd.to_numeric(direct_matches[piece_length_col], errors="coerce").round(2)
            direct_matches[metric_names["proportion"]] = 1.0
            direct_matches[metric_names["piece_length"]] = piece_lengths
            direct_matches[metric_names["zone_length"]] = piece_lengths
            direct_matches[metric_names["piece_surface"]] = polygon_areas
            direct_matches[metric_names["zone_surface"]] = polygon_areas

    pairs = gpd.sjoin(
        crossing_polygons[["__poly_idx", "geometry"]],
        zones_indexed[["__zone_idx", "geometry"]],
        how="inner",
        predicate="intersects",
    )
    pairs = pairs.rename(columns={"geometry": "polygon_geometry"}).drop(columns=["index_right"])

    if len(pairs) == 0 and (direct_matches is None or direct_matches.empty):
        if prefilter_zones_to_network_bbox and len(zones_proj):
            voided = _build_cascade_void_rows(
                zones_indexed,
                polygon_attr_cols,
                zone_attr_cols,
                output_epsg,
                metric_names=metric_names,
                zone_label=zone_label,
            )
            if output_path:
                from osm_chordify.utils.io import save_geodataframe
                save_geodataframe(voided, output_path)
                _write_intersection_cache_metadata(output_path, cascade_cache_metadata)
            return voided
        return gpd.GeoDataFrame(
            columns=list(polygon_attr_cols) + [
                metric_names["proportion"],
                metric_names["piece_length"],
                metric_names["zone_length"],
                metric_names["piece_surface"],
                metric_names["zone_surface"],
                "geometry",
            ],
            geometry="geometry",
            crs=output_epsg,
        )

    matched = pairs.merge(
        polygons_indexed.drop(columns=["geometry"]),
        on="__poly_idx",
        how="left",
    ).merge(
        zones_indexed.rename(columns={"geometry": "zone_geometry"}),
        on="__zone_idx",
        how="left",
    )

    exact_parts = []
    if len(pairs):
        exact_parts = _compute_chunked_exact_parts(
            matched=matched,
            epsg=polygons_epsg,
            compute_func=_compute_polygon_cascade_chunk,
            desc="Computing polygon intersections",
            postfix_template=lambda pieces: f"hit_zones={pairs['__zone_idx'].nunique()}, pieces={pieces}, poly_hits={len(pairs)}",
            metric_names=metric_names,
            length_col=piece_length_col,
            enable_parallel=_is_county_like_zones(zones_indexed),
        )

    frames = []
    if direct_matches is not None and not direct_matches.empty:
        frames.append(direct_matches)
    if exact_parts:
        frames.extend(exact_parts)

    if not frames:
        return gpd.GeoDataFrame(
            columns=list(polygon_attr_cols) + [
                metric_names["proportion"],
                metric_names["piece_length"],
                metric_names["zone_length"],
                metric_names["piece_surface"],
                metric_names["zone_surface"],
                "geometry",
            ],
            geometry="geometry",
            crs=output_epsg,
        )

    matched_all = gpd.GeoDataFrame(
        frames[0] if len(frames) == 1 else pd.concat(frames, ignore_index=True),
        geometry="geometry",
        crs=polygons_epsg,
    )
    result_gdf = _build_cascade_result_gdf(
        matched_all,
        polygon_attr_cols,
        zone_attr_cols,
        polygons_epsg,
        metric_names=metric_names,
        zone_label=zone_label,
    )
    if prefilter_zones_to_network_bbox and "__zone_idx" in matched_all.columns:
        matched_zone_indices = set(matched_all["__zone_idx"].unique())
        missing_zone_rows = zones_indexed.loc[~zones_indexed["__zone_idx"].isin(matched_zone_indices)]
        if len(missing_zone_rows):
            voided = _build_cascade_void_rows(
                missing_zone_rows,
                polygon_attr_cols,
                zone_attr_cols,
                polygons_epsg,
                metric_names=metric_names,
                zone_label=zone_label,
            )
            result_gdf = gpd.GeoDataFrame(pd.concat([result_gdf, voided], ignore_index=True), geometry="geometry", crs=polygons_epsg)

    if output_epsg != polygons_epsg:
        result_gdf = result_gdf.to_crs(epsg=output_epsg)

    if output_path:
        from osm_chordify.utils.io import save_geodataframe
        save_geodataframe(result_gdf, output_path)
        _write_intersection_cache_metadata(output_path, cascade_cache_metadata)
        logger.info("Saved polygon cascade intersection to %s", output_path)

    return result_gdf


def spatial_left_join_with_zones(
    gdf,
    gdf_epsg,
    zones,
    output_path=None,
    output_epsg=None,
    zone_label=None,
):
    """Spatially left-join a geometry layer with zone polygons.

    All rows from *gdf* are preserved. Zone attributes are appended with
    prefixed names and remain null where no zone intersects the input geometry.
    """
    _require_projected_epsg(gdf_epsg)
    input_gdf = _load_zones(gdf)
    zones_gdf = _load_zones(zones)

    if output_epsg is None:
        output_epsg = gdf_epsg

    logger.info("Projecting geometries to EPSG:%d", gdf_epsg)
    input_proj = _project_if_needed(input_gdf, gdf_epsg)
    zones_proj = _project_if_needed(zones_gdf, gdf_epsg)

    input_attr_cols = [c for c in input_gdf.columns if c != "geometry"]
    zone_attr_cols = [c for c in zones_gdf.columns if c != "geometry"]
    input_indexed = input_proj.reset_index(drop=True).reset_index(names="__input_idx")
    zones_indexed = zones_proj.reset_index(drop=True).reset_index(names="__zone_idx")

    joined = gpd.sjoin(
        input_indexed,
        zones_indexed[["__zone_idx", "geometry"] + zone_attr_cols],
        how="left",
        predicate="intersects",
    ).drop(columns=["index_right"])

    result_data = {col: joined[col] for col in input_attr_cols}
    result_data["geometry"] = joined["geometry"]
    result_cols = list(input_attr_cols) + ["geometry"]
    existing_keys = set(result_cols)

    for col in zone_attr_cols:
        out_col = _zone_output_name(col, existing_keys, prefix=_label_prefix(zone_label) or "zone")
        result_data[out_col] = joined[col]
        if out_col not in result_cols:
            result_cols.append(out_col)
            existing_keys.add(out_col)

    result_gdf = gpd.GeoDataFrame(result_data, geometry="geometry", crs=gdf_epsg)[result_cols]
    if output_epsg != gdf_epsg:
        result_gdf = result_gdf.to_crs(epsg=output_epsg)

    if output_path:
        from osm_chordify.utils.io import save_geodataframe
        save_geodataframe(result_gdf, output_path)
        logger.info("Saved spatial left join to %s", output_path)

    return result_gdf


def intersect_zones_with_zones(
    zones_a,
    zones_a_epsg,
    zones_b,
    zones_b_epsg,
    output_path=None,
    output_epsg=None,
):
    """Intersect one polygon zone set with another.

    Parameters
    ----------
    zones_a : str, os.PathLike, or gpd.GeoDataFrame
        First zone layer.
    zones_a_epsg : int
        EPSG code for the first zone layer CRS.
    zones_b : str, os.PathLike, or gpd.GeoDataFrame
        Second zone layer.
    zones_b_epsg : int
        EPSG code for the second zone layer CRS.
    output_path : str, optional
        If provided, save the result to this file.
    output_epsg : int, optional
        EPSG code for the output CRS. Defaults to ``zones_a_epsg``.

    Returns
    -------
    gpd.GeoDataFrame
        Polygon intersections with prefixed attributes from both inputs.
    """
    zones_a_gdf = _load_zones(zones_a)
    zones_b_gdf = _load_zones(zones_b)

    if output_epsg is None:
        output_epsg = zones_a_epsg

    zones_a_proj = _project_if_needed(zones_a_gdf, output_epsg)
    zones_b_proj = _project_if_needed(zones_b_gdf, output_epsg if zones_b_epsg != output_epsg else zones_b_epsg)

    a_cols = {col: f"zone_a_{col}" for col in zones_a_proj.columns if col != "geometry"}
    b_cols = {col: f"zone_b_{col}" for col in zones_b_proj.columns if col != "geometry"}
    zones_a_prefixed = zones_a_proj.rename(columns=a_cols)
    zones_b_prefixed = zones_b_proj.rename(columns=b_cols)

    result = gpd.overlay(zones_a_prefixed, zones_b_prefixed, how="intersection", keep_geom_type=False)
    if len(result):
        result = result.loc[result.geometry.geom_type.isin(POLYGON_GEOMETRY_TYPES)].copy()
    else:
        result = gpd.GeoDataFrame(columns=["geometry"], geometry="geometry", crs=f"EPSG:{output_epsg}")

    if output_path:
        from osm_chordify.utils.io import save_geodataframe
        save_geodataframe(result, output_path)

    return result
