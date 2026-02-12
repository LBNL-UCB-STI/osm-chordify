"""Shared file-saving helpers for GeoDataFrames and DataFrames."""

import logging
import os

import geopandas as gpd
import pandas as pd

logger = logging.getLogger(__name__)


def save_geodataframe(gdf, output_path):
    """Save a GeoDataFrame, inferring the format from the file extension.

    Supported formats: ``.gpkg``, ``.geojson``, ``.shp``, ``.csv``.
    For CSV output the geometry is exported as a ``geometry_wkt`` column.

    Parameters
    ----------
    gdf : gpd.GeoDataFrame
        GeoDataFrame to save.
    output_path : str
        Destination file path.
    """
    output_dir = os.path.dirname(output_path)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    ext = os.path.splitext(output_path)[1].lower()
    if ext == ".gpkg":
        gdf.to_file(output_path, driver="GPKG")
    elif ext == ".geojson":
        gdf.to_file(output_path, driver="GeoJSON")
    elif ext == ".shp":
        gdf.to_file(output_path)
    elif ext == ".csv":
        gdf = gdf.copy()
        gdf["geometry_wkt"] = gdf.geometry.apply(lambda g: g.wkt)
        pd.DataFrame(gdf.drop(columns="geometry")).to_csv(output_path, index=False)
    else:
        logger.info("Unrecognized output format: %s, defaulting to GPKG", ext)
        gdf.to_file(output_path, driver="GPKG")

    logger.info("Saved to %s", output_path)


def save_dataframe(df, output_path):
    """Save a DataFrame (or GeoDataFrame) with optional geometry handling.

    For spatial formats (``.gpkg``, ``.geojson``) the DataFrame is promoted to
    a GeoDataFrame when a ``geometry`` column is present.  For ``.csv`` the
    geometry is exported as WKT and dropped.

    Parameters
    ----------
    df : pd.DataFrame or gpd.GeoDataFrame
        Data to save.
    output_path : str
        Destination file path.
    """
    output_dir = os.path.dirname(output_path)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    ext = os.path.splitext(output_path)[1].lower()

    if ext in (".gpkg", ".geojson"):
        if not isinstance(df, gpd.GeoDataFrame) and "geometry" in df.columns:
            df = gpd.GeoDataFrame(df, geometry="geometry")
        driver = "GPKG" if ext == ".gpkg" else "GeoJSON"
        df.to_file(output_path, driver=driver)
    elif ext == ".csv":
        df = df.copy()
        if "geometry" in df.columns:
            df["geometry_wkt"] = df["geometry"].apply(
                lambda g: g.wkt if g else None
            )
            df = df.drop(columns="geometry")
        df.to_csv(output_path, index=False)
    else:
        logger.info("Unrecognized output format: %s, defaulting to CSV", ext)
        df.to_csv(output_path, index=False)

    logger.info("Saved to %s", output_path)
