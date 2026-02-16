"""Map a tabular network to a polygon-OSM intersection via a shared ID."""

import logging

import pandas as pd

logger = logging.getLogger(__name__)


def map_network_to_intersection(
    network_df,
    intersection_gdf,
    polygon_id_col="polygon_id",
    network_osm_id_col="attributeOrigId",
    network_length_col="linkLength",
    network_link_id_col="linkId",
):
    """Join a tabular network to a polygon-OSM intersection result.

    The join is performed on the OSM ID present in both datasets.  For each
    matched row the proportional network length is computed as
    ``link_length * proportion``.

    Parameters
    ----------
    network_df : pd.DataFrame
        Network table (CSV/parquet).  Must contain *network_osm_id_col*,
        *network_length_col*, and *network_link_id_col*.
    intersection_gdf : pd.DataFrame or gpd.GeoDataFrame
        Output of :func:`osm_chordify.osm.intersect.intersect_osm_with_zones`.
        Must contain ``osm_id``, ``proportion``, and *polygon_id_col*.
    polygon_id_col : str
        Name of the polygon-ID column in *intersection_gdf*.
    network_osm_id_col : str
        Column in *network_df* that holds the OSM ID for joining.
    network_length_col : str
        Column in *network_df* that holds the link length.
    network_link_id_col : str
        Column in *network_df* that holds the link ID.

    Returns
    -------
    pd.DataFrame
        Merged result with ``proportional_network_length`` and a
        ``network_polygon_id`` composite key.

    Raises
    ------
    ValueError
        If required columns are missing from the input DataFrames.
    """
    for col in (network_osm_id_col, network_length_col, network_link_id_col):
        if col not in network_df.columns:
            raise ValueError(f"Network DataFrame is missing '{col}' column")
    if "osm_id" not in intersection_gdf.columns:
        raise ValueError("Intersection DataFrame is missing 'osm_id' column")
    if polygon_id_col not in intersection_gdf.columns:
        raise ValueError(
            f"Intersection DataFrame is missing '{polygon_id_col}' column"
        )

    network_df = network_df.dropna(subset=[network_osm_id_col])

    network_df = network_df.copy()
    network_df[network_osm_id_col] = network_df[network_osm_id_col].astype(int)
    intersection_gdf = intersection_gdf.copy()
    intersection_gdf["osm_id"] = intersection_gdf["osm_id"].astype(int)

    logger.info("Merging network data with polygon-OSM intersection")
    merged = pd.merge(
        network_df,
        intersection_gdf,
        left_on=network_osm_id_col,
        right_on="osm_id",
        how="inner",
    )
    logger.info("Merged result has %d rows", len(merged))

    merged["proportional_network_length"] = (
        merged[network_length_col] * merged["proportion"]
    )
    merged["network_polygon_id"] = (
        merged[network_link_id_col].astype(str)
        + "-"
        + merged[polygon_id_col].astype(str)
    )

    return merged
