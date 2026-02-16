"""Map a tabular network to a road-network/zone intersection via a shared ID."""

import logging

import pandas as pd

logger = logging.getLogger(__name__)


def map_network_to_intersection(
    network_df,
    intersection_gdf,
    network_osm_id_col="attributeOrigId",
    intersection_osm_id_col="edge_osm_id",
):
    """Join a tabular network to an intersection result on OSM ID.

    All columns from both inputs are included in the output.

    Parameters
    ----------
    network_df : pd.DataFrame
        Network table (e.g. BEAM network CSV).
    intersection_gdf : pd.DataFrame or gpd.GeoDataFrame
        Output of :func:`osm_chordify.osm.intersect.intersect_road_network_with_zones`.
    network_osm_id_col : str
        Column in *network_df* that holds the OSM ID for joining.
    intersection_osm_id_col : str
        Column in *intersection_gdf* that holds the OSM ID for joining.

    Returns
    -------
    pd.DataFrame
        Merged result with all columns from both inputs.
    """
    if network_osm_id_col not in network_df.columns:
        raise ValueError(f"Network DataFrame is missing '{network_osm_id_col}' column")
    if intersection_osm_id_col not in intersection_gdf.columns:
        raise ValueError(
            f"Intersection DataFrame is missing '{intersection_osm_id_col}' column"
        )

    network_df = network_df.dropna(subset=[network_osm_id_col]).copy()
    network_df[network_osm_id_col] = network_df[network_osm_id_col].astype(int)

    intersection_gdf = intersection_gdf.copy()
    intersection_gdf[intersection_osm_id_col] = intersection_gdf[intersection_osm_id_col].astype(int)

    logger.info("Merging network data with intersection result")
    merged = pd.merge(
        network_df,
        intersection_gdf,
        left_on=network_osm_id_col,
        right_on=intersection_osm_id_col,
        how="inner",
    )
    logger.info("Merged result has %d rows", len(merged))

    return merged
