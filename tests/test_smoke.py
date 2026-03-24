import importlib.util
import tomllib
from pathlib import Path

import geopandas as gpd
import pandas as pd
from shapely.geometry import LineString

from osm_chordify import (
    build_osm_by_pop_density,
    export_network,
    intersect_road_network_with_zones,
)
from osm_chordify.osm.graph import download_and_prepare_osm_network
from osm_chordify.osm.intersect import load_osm_edges
from osm_chordify.utils.geo import name_osm_network, create_osm_highway_filter
from osm_chordify.utils.network import map_network_to_intersection
from osm_chordify.utils.io import save_geodataframe, save_dataframe


def test_public_api_imports():
    """Verify the public API is importable and callable."""
    assert callable(build_osm_by_pop_density)
    assert callable(export_network)
    assert callable(intersect_road_network_with_zones)
    assert callable(download_and_prepare_osm_network)


def test_module_imports():
    """Verify the modules are importable and callable."""
    assert callable(load_osm_edges)
    assert callable(intersect_road_network_with_zones)
    assert callable(map_network_to_intersection)
    assert callable(save_geodataframe)
    assert callable(save_dataframe)
    assert callable(name_osm_network)
    assert callable(create_osm_highway_filter)


def test_osmium_dependency_declared_in_project_metadata():
    """Verify the Python PBF export dependency is declared for installation."""
    pyproject = Path(__file__).resolve().parents[1] / "pyproject.toml"
    data = tomllib.loads(pyproject.read_text())
    dependencies = data["project"]["dependencies"]
    assert "osmium" in dependencies


def test_name_osm_network():
    """Verify name_osm_network produces expected output."""
    layers = {
        "main": {"geo_level": "county"},
        "residential": {"min_density_per_km2": 5500, "geo_level": "cbg"},
    }
    name = name_osm_network("sfbay", layers, True)
    assert name == "sfbay-cbg5500-strongConn-network"

    name_weak = name_osm_network("sfbay", layers, False)
    assert "weakConn" in name_weak

    layers_ferry = {**layers, "ferry": {"geo_level": "county"}}
    name_ferry = name_osm_network("seattle", layers_ferry, False)
    assert "-ferry-" in name_ferry


def test_save_helpers_support_parquet(tmp_path):
    if importlib.util.find_spec("pyarrow") is None:
        return

    geodf = gpd.GeoDataFrame(
        {
            "osm_id": [1],
            "value": [2.5],
            "geometry": [LineString([(0, 0), (1, 1)])],
        },
        geometry="geometry",
        crs="EPSG:4326",
    )

    geo_path = tmp_path / "edges.parquet"
    save_geodataframe(geodf, geo_path)
    loaded_geo = gpd.read_parquet(geo_path)
    assert len(loaded_geo) == 1
    assert "geometry" in loaded_geo.columns

    df_path = tmp_path / "mapping.parquet"
    save_dataframe(pd.DataFrame(geodf), df_path)
    loaded_df = gpd.read_parquet(df_path)
    assert len(loaded_df) == 1
    assert "geometry" in loaded_df.columns
