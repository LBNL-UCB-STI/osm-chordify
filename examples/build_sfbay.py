"""Build an OSM network for the San Francisco Bay Area (9 counties)."""

import os
import sys
from pathlib import Path

from osm_chordify import build_osm_by_pop_density

EXAMPLES_DIR = Path(__file__).resolve().parent
if str(EXAMPLES_DIR) not in sys.path:
    sys.path.insert(0, str(EXAMPLES_DIR))

from common import (
    base_osm_config,
    expand_work_dir,
    highway_filter,
    parse_build_args,
    validate_built_network,
)

work_dir = expand_work_dir("sfbay")


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
    "taz_shp": "geo/shp/sfbay-tazs-epsg-26910.shp",
    "taz_id": "taz1454",
    "cbg_id": "GEOID",
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
    args = parse_build_args(work_dir, geo_config["taz_shp"])
    if args.census_api_key:
        os.environ["CENSUS_API_KEY"] = args.census_api_key

    geo_config["taz_shp"] = args.taz_shp

    result = build_osm_by_pop_density(
        work_dir=args.work_dir,
        osm_config=osm_config,
        area_config=area_config,
        geo_config=geo_config,
    )

    if not args.skip_validation:
        summary = validate_built_network(
            result["graph"],
            osm_path=result["exported"].get("osm"),
        )
        print(summary)
