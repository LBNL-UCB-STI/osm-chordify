"""Intersect an OSM network with zone polygons and map a BEAM network."""

import os

from osm_chordify import intersect_road_network_with_zones, map_osm_with_beam_network

work_dir = os.path.expanduser("~/Workspace/Simulation/sfbay")
utm_epsg = 26910
osm_name = "sfbay-main-residential"
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
