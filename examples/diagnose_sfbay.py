"""Run diagnostics on a downloaded OSM PBF file."""

import os

from osm_chordify import diagnose_osm

work_dir = os.path.expanduser("~/Workspace/Simulation/sfbay")
osm_name = "sfbay-main-residential"
osm_dir = f"{work_dir}/network/{osm_name}"

pbf_path = f"{osm_dir}/{osm_name}.osm.pbf"
utm_epsg = 26910

if __name__ == "__main__":
    diagnose_osm(pbf_path, utm_epsg)
