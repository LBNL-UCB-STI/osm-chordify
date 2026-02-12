from osm_chordify import download_osm_network, export_network, get_area_config


def test_public_api_imports():
    """Verify the public API is importable and callable."""
    assert callable(download_osm_network)
    assert callable(export_network)
    assert callable(get_area_config)


def test_get_area_config():
    """Verify get_area_config returns a valid config dict."""
    config = get_area_config("sfbay")
    assert config["area"]["name"] == "sfbay"
    assert "network" in config
    assert "geo" in config
