"""Intersect an OSM network with zone polygons and map a BEAM network."""

import os
import sys
from pathlib import Path

EXAMPLES_DIR = Path(__file__).resolve().parent
if str(EXAMPLES_DIR) not in sys.path:
    sys.path.insert(0, str(EXAMPLES_DIR))

from _bootstrap import bootstrap_example_paths

bootstrap_example_paths(__file__)

from osm_chordify import intersect_road_network_with_zones, map_osm_with_beam_network
from osm_chordify.utils.geo import name_osm_network

work_dir = os.getcwd()
utm_epsg = 26910  # NAD83 / UTM zone 10N
area_name = "sfbay"
graph_layers = {
    "backbone": {
        "layer_role": "backbone",
        "geo_level": "county",
    },
    "connector": {
        "layer_role": "connector",
        "geo_level": "county",
    },
    "residential": {
        "layer_role": "residential",
        "min_density_per_km2": 5500,
        "geo_level": "cbg",
    },
}
strongly_connected_components = True
osm_name = name_osm_network(
    area_name,
    graph_layers,
    strongly_connected_components,
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
