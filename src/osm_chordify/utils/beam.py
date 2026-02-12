"""BEAM network ↔ polygon-OSM intersection mapping."""

import logging

import pandas as pd

logger = logging.getLogger(__name__)


def map_beam_network_to_polygon_intersection(
    network_df,
    intersection_gdf,
    polygon_id_col="polygon_id",
    beam_osm_id_col="attributeOrigId",
    beam_length_col="linkLength",
    beam_link_id_col="linkId",
):
    """Join a BEAM network DataFrame to a polygon-OSM intersection result.

    The join is performed on the OSM ID present in both datasets.  For each
    matched row the proportional network length is computed as
    ``beam_length * proportion``.

    Parameters
    ----------
    network_df : pd.DataFrame
        BEAM network table.  Must contain *beam_osm_id_col*,
        *beam_length_col*, and *beam_link_id_col*.
    intersection_gdf : pd.DataFrame or gpd.GeoDataFrame
        Output of :func:`osm_chordify.osm.intersect.intersect_edges_with_polygons`.
        Must contain ``osm_id``, ``proportion``, and *polygon_id_col*.
    polygon_id_col : str
        Name of the polygon-ID column in *intersection_gdf* (default
        ``"polygon_id"``).
    beam_osm_id_col : str
        Column in *network_df* that holds the OSM ID (default
        ``"attributeOrigId"``).
    beam_length_col : str
        Column in *network_df* that holds the link length (default
        ``"linkLength"``).
    beam_link_id_col : str
        Column in *network_df* that holds the link ID (default ``"linkId"``).

    Returns
    -------
    pd.DataFrame
        Merged result with an added ``proportional_network_length`` column
        and a ``network_polygon_id`` composite key.

    Raises
    ------
    ValueError
        If required columns are missing from the input DataFrames.
    """
    # Validate inputs
    for col in (beam_osm_id_col, beam_length_col, beam_link_id_col):
        if col not in network_df.columns:
            raise ValueError(f"Network DataFrame is missing '{col}' column")
    if "osm_id" not in intersection_gdf.columns:
        raise ValueError("Intersection DataFrame is missing 'osm_id' column")
    if polygon_id_col not in intersection_gdf.columns:
        raise ValueError(
            f"Intersection DataFrame is missing '{polygon_id_col}' column"
        )

    # Drop rows without an OSM ID
    network_df = network_df.dropna(subset=[beam_osm_id_col])

    # Ensure compatible join types
    network_df = network_df.copy()
    network_df[beam_osm_id_col] = network_df[beam_osm_id_col].astype(int)
    intersection_gdf = intersection_gdf.copy()
    intersection_gdf["osm_id"] = intersection_gdf["osm_id"].astype(int)

    logger.info("Merging network data with polygon-OSM intersection")
    merged = pd.merge(
        network_df,
        intersection_gdf,
        left_on=beam_osm_id_col,
        right_on="osm_id",
        how="inner",
    )
    logger.info("Merged result has %d rows", len(merged))

    merged["proportional_network_length"] = (
        merged[beam_length_col] * merged["proportion"]
    )
    merged["network_polygon_id"] = (
        merged[beam_link_id_col].astype(str)
        + "-"
        + merged[polygon_id_col].astype(str)
    )

    return merged
