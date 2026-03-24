"""Public package exports for osm-chordify."""

from importlib import import_module

__all__ = [
    "build_osm_by_pop_density",
    "create_osm_highway_filter",
    "export_network",
    "intersect_road_network_with_zones",
    "map_osm_with_beam_network",
    "match_road_network_geometries",
    "diagnose_osm",
    "intersect_road_network_with_county_zones",
]

_EXPORTS = {
    "build_osm_by_pop_density": ("osm_chordify.main", "build_osm_by_pop_density"),
    "map_osm_with_beam_network": ("osm_chordify.main", "map_osm_with_beam_network"),
    "match_road_network_geometries": ("osm_chordify.main", "match_road_network_geometries"),
    "diagnose_osm": ("osm_chordify.main", "diagnose_osm"),
    "intersect_road_network_with_county_zones": (
        "osm_chordify.main",
        "intersect_road_network_with_county_zones",
    ),
    "intersect_road_network_with_zones": (
        "osm_chordify.osm.intersect",
        "intersect_road_network_with_zones",
    ),
    "export_network": ("osm_chordify.osm.export", "export_network"),
    "create_osm_highway_filter": ("osm_chordify.utils.geo", "create_osm_highway_filter"),
}


def __getattr__(name):
    """Lazily resolve public exports to avoid eager submodule imports."""
    if name not in _EXPORTS:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

    module_name, attr_name = _EXPORTS[name]
    value = getattr(import_module(module_name), attr_name)
    globals()[name] = value
    return value
