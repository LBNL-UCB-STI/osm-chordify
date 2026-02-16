"""Intersect an OSM network with a polygon grid and map a BEAM network."""

import os

from osm_chordify import intersect_network_geom_with_zones, map_osm_with_beam_network

work_dir = os.path.expanduser("~/Workspace/Simulation/sfbay")
utm_epsg = 26910
osm_name = "sfbay-main-residential"
osm_dir = f"{work_dir}/network/{osm_name}"

# --- Step 1: Intersect OSM edges with polygon grid ---

osm_geojson = f"{osm_dir}/{osm_name}.osm.geojson"
osm_gpkg = f"{osm_dir}/{osm_name}.gpkg"
grid_path = f"{work_dir}/inmap/ISRM/isrm_polygon.shp"
id_col = "isrm"

out_dir = f"{work_dir}/polygon-{osm_name}"
os.makedirs(out_dir, exist_ok=True)
intersection_path = f"{out_dir}/polygon-{osm_name}.geojson"

if not os.path.exists(intersection_path):
    intersect_network_geom_with_zones(
        grid_path=grid_path,
        id_col=id_col,
        osm_geojson_path=osm_geojson,
        osm_gpkg_path=osm_gpkg,
        epsg_utm=utm_epsg,
        output_path=intersection_path,
    )

# --- Step 2: Map a BEAM network to the intersection ---

net_csv = f"{osm_dir}/network.csv.gz"
mapping_output = f"{out_dir}/polygon-network-mapping.geojson"

if not os.path.exists(mapping_output):
    map_osm_with_beam_network(
        network_path=net_csv,
        intersection_path=intersection_path,
        id_col="polygon_id",
        osm_id_col="attributeOrigId",
        length_col="linkLength",
        link_id_col="linkId",
        output_path=mapping_output,
    )
