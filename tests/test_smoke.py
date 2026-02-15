from osm_chordify import download_osm_network, export_network, intersect_edges_with_polygons, download_and_build_osm
from osm_chordify.osm.intersect import load_osm_edges
from osm_chordify.utils.geo import name_osm_network, create_osm_highway_filter
from osm_chordify.utils.network import map_network_to_intersection
from osm_chordify.utils.io import save_geodataframe, save_dataframe


def test_public_api_imports():
    """Verify the public API is importable and callable."""
    assert callable(download_osm_network)
    assert callable(export_network)
    assert callable(intersect_edges_with_polygons)
    assert callable(download_and_build_osm)


def test_module_imports():
    """Verify the modules are importable and callable."""
    assert callable(load_osm_edges)
    assert callable(intersect_edges_with_polygons)
    assert callable(map_network_to_intersection)
    assert callable(save_geodataframe)
    assert callable(save_dataframe)
    assert callable(name_osm_network)
    assert callable(create_osm_highway_filter)


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
