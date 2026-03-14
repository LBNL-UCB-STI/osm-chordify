"""Intersect an OSM network with zone polygons and map a BEAM network."""

import importlib.util
import os
from pathlib import Path

from osm_chordify import intersect_road_network_with_zones, map_osm_with_beam_network
from osm_chordify.utils.geo import name_osm_network


def _load_build_sfbay():
    path = Path(__file__).with_name("build_sfbay.py")
    spec = importlib.util.spec_from_file_location("build_sfbay_example", path)
    module = importlib.util.module_from_spec(spec)
    assert spec is not None and spec.loader is not None
    spec.loader.exec_module(module)
    return module


build_sfbay = _load_build_sfbay()

work_dir = os.getcwd()
utm_epsg = build_sfbay.geo_config["utm_epsg"]
osm_name = name_osm_network(
    build_sfbay.area_config["name"],
    build_sfbay.osm_config["graph_layers"],
    build_sfbay.osm_config["strongly_connected_components"],
)
osm_dir = f"{work_dir}/network/{osm_name}"

# --- Step 1: Intersect OSM edges with zone polygons ---

osm_gpkg = f"{osm_dir}/{osm_name}.gpkg"
zone_polygon_path = f"{work_dir}/inmap/ISRM/isrm_polygon.shp"

out_dir = f"{work_dir}/polygon-{osm_name}"
os.makedirs(out_dir, exist_ok=True)
intersection_path = f"{out_dir}/polygon-{osm_name}.geojson"

if not os.path.exists(intersection_path):
    intersect_road_network_with_zones(
        road_network=osm_gpkg,
        road_network_epsg=utm_epsg,
        zones=zone_polygon_path,
        zones_epsg=utm_epsg,
        output_path=intersection_path,
    )

# --- Step 2: Map a BEAM network to the intersection ---

net_csv = f"{osm_dir}/network.csv.gz"
mapping_output = f"{out_dir}/polygon-network-mapping.geojson"

if not os.path.exists(mapping_output):
    map_osm_with_beam_network(
        osm_path=osm_gpkg,
        network_path=net_csv,
        network_osm_id_col="attributeOrigId",
        output_path=mapping_output,
    )
