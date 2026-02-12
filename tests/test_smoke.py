from osm_chordify import download_osm_network, export_network, get_area_config, intersect_edges_with_polygons
from osm_chordify.osm.intersect import load_osm_edges
from osm_chordify.utils.beam import map_beam_network_to_polygon_intersection
from osm_chordify.utils.io import save_geodataframe, save_dataframe


def test_public_api_imports():
    """Verify the public API is importable and callable."""
    assert callable(download_osm_network)
    assert callable(export_network)
    assert callable(get_area_config)
    assert callable(intersect_edges_with_polygons)


def test_new_module_imports():
    """Verify the new modules are importable and callable."""
    assert callable(load_osm_edges)
    assert callable(intersect_edges_with_polygons)
    assert callable(map_beam_network_to_polygon_intersection)
    assert callable(save_geodataframe)
    assert callable(save_dataframe)


def test_get_area_config():
    """Verify get_area_config returns a valid config dict."""
    config = get_area_config("sfbay")
    assert config["area"]["name"] == "sfbay"
    assert "network" in config
    assert "geo" in config
