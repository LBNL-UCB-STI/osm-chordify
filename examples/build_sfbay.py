"""Build an OSM network for the San Francisco Bay Area (9 counties)."""

import sys
from pathlib import Path

EXAMPLES_DIR = Path(__file__).resolve().parent
if str(EXAMPLES_DIR) not in sys.path:
    sys.path.insert(0, str(EXAMPLES_DIR))

from common import (
    base_osm_config,
    highway_filter,
    run_example_build,
)

work_dir = str(Path.cwd() / "output")


area_config = {
    "name": "sfbay",
    "state_fips": "06",
    "county_fips": [
        "001", "013", "041", "055", "075", "081", "085", "095", "097",
    ],
    "census_year": 2018,
}

geo_config = {
    "utm_epsg": 26910,  # NAD83 / UTM zone 10N
}

osm_config = base_osm_config(tolerance=2, strongly_connected_components=True)
osm_config["graph_layers"] = {
    "main": {
        "geo_level": "county",
        "custom_filter": highway_filter(),
        "buffer_zone_in_meters": 200,
    },
    "residential": {
        "min_density_per_km2": 5500,
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
