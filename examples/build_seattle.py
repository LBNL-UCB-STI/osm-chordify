"""Build an OSM network for the Seattle metro area (King, Kitsap, Pierce, Snohomish)."""

import sys
from pathlib import Path

EXAMPLES_DIR = Path(__file__).resolve().parent
if str(EXAMPLES_DIR) not in sys.path:
    sys.path.insert(0, str(EXAMPLES_DIR))

from _bootstrap import bootstrap_example_paths

bootstrap_example_paths(__file__)

from common import (
    backbone_highway_filter,
    base_osm_config,
    connector_highway_filter,
    highway_filter,
    run_example_build,
)

work_dir = str(Path.cwd() / "output")


area_config = {
    "name": "seattle",
    "state_fips": "53",
    "county_fips": ["033", "035", "053", "061"],  # King, Kitsap, Pierce, Snohomish
    "census_year": 2018,
}

geo_config = {
    "utm_epsg": 32048,  # NAD83 / Washington North
}

osm_config = base_osm_config(tolerance=2, strongly_connected_components=False)
osm_config["graph_layers"] = {
    "backbone": {
        "layer_role": "backbone",
        "geo_level": "county",
        "custom_filter": backbone_highway_filter(),
        "buffer_zone_in_meters": 5000,
        "protected_backbone": True,
    },
    "connector": {
        "layer_role": "connector",
        "geo_level": "county",
        "custom_filter": connector_highway_filter(),
        "buffer_zone_in_meters": 180,
    },
    "ferry": {
        "layer_role": "ferry",
        "geo_level": "county",
        "custom_filter": '["route"="ferry"]',
        "buffer_zone_in_meters": 10000,
    },
    "residential": {
        "layer_role": "residential",
        "min_density_per_km2": 120,
        "geo_level": "cbg",
        "custom_filter": highway_filter(include_residential=True),
        "buffer_zone_in_meters": 20,
    },
}

if __name__ == "__main__":
    run_example_build(
        default_output_dir=work_dir,
        osm_config=osm_config,
        area_config=area_config,
        geo_config=geo_config,
    )
