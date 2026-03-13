"""Run diagnostics on a downloaded OSM PBF file."""

import importlib.util
import os
from pathlib import Path

from osm_chordify import diagnose_osm
from osm_chordify.utils.geo import name_osm_network


def _load_build_sfbay():
    path = Path(__file__).with_name("build_sfbay.py")
    spec = importlib.util.spec_from_file_location("build_sfbay_example", path)
    module = importlib.util.module_from_spec(spec)
    assert spec is not None and spec.loader is not None
    spec.loader.exec_module(module)
    return module


build_sfbay = _load_build_sfbay()
work_dir = build_sfbay.work_dir
area_config = build_sfbay.area_config
geo_config = build_sfbay.geo_config
osm_config = build_sfbay.osm_config

osm_name = name_osm_network(
    area_config["name"],
    osm_config["graph_layers"],
    osm_config["strongly_connected_components"],
)
osm_dir = f"{work_dir}/network/{osm_name}"

pbf_path = f"{osm_dir}/{osm_name}.osm.pbf"
utm_epsg = geo_config["utm_epsg"]

if __name__ == "__main__":
    if os.environ.get("OSM_CHORDIFY_SKIP_DIAGNOSE") != "1":
        diagnose_osm(pbf_path, utm_epsg)
