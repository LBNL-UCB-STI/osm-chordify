"""Build an OSM network for the San Francisco Bay Area (9 counties)."""

import os

import osmnx as ox

from osm_chordify import build_osm_by_pop_density

work_dir = os.path.expanduser("~/Workspace/Simulation/sfbay")

osm_highways = [
    "motorway", "motorway_link", "trunk", "trunk_link",
    "primary", "primary_link", "secondary", "secondary_link",
    "tertiary", "tertiary_link", "unclassified",
]
osm_residential = ["residential"]


def _highway_filter(types):
    return f'["highway"~"{"|".join(types)}"]'


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

osm_config = {
    "osmnx_settings": {
        "log_console": True,
        "use_cache": True,
        "cache_only_mode": False,
        "all_oneway": True,
        "requests_timeout": 180,
        "overpass_memory": None,
        "max_query_area_size": 50 * 1000 * 50 * 1000,
        "overpass_rate_limit": False,
        "overpass_max_attempts": 3,
        "useful_tags_way": list(ox.settings.useful_tags_way) + [
            "maxweight", "hgv", "maxweight:hgv", "maxlength",
            "motorcar", "motor_vehicle", "goods", "truck",
        ],
        "overpass_url": "https://overpass-api.de/api",
    },
    "weight_limits": {"unit": "lbs", "mdv_max": 26000, "hdv_max": 80000},
    "download_enabled": True,
    "tolerance": 2,
    "strongly_connected_components": True,
    "graph_layers": {
        "main": {
            "geo_level": "county",
            "custom_filter": _highway_filter(list(set(osm_highways))),
            "buffer_zone_in_meters": 200,
        },
        "residential": {
            "min_density_per_km2": 5500,
            "geo_level": "cbg",
            "custom_filter": _highway_filter(
                list(set(osm_highways) | set(osm_residential))
            ),
            "buffer_zone_in_meters": 20,
        },
    },
}

if __name__ == "__main__":
    build_osm_by_pop_density(
        work_dir=work_dir,
        area_name=area_config["name"],
        osm_config=osm_config,
        area_config=area_config,
        geo_config=geo_config,
    )
