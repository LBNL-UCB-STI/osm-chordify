"""Export a network graph to one or more file formats."""

import importlib
import logging
import os
import pickle

import osmnx as ox

from osm_chordify.osm.xml import save_graph_xml

logger = logging.getLogger(__name__)

DEFAULT_EDGE_TAGS = [
    'highway', 'lanes', 'maxspeed', 'name', 'oneway', 'length',
    'tunnel', 'bridge', 'junction', 'edge_id', 'access', 'osm_id',
    'motor_vehicle', 'vehicle', 'motorcar',
    'access:car', 'access:vehicle', 'access:motor_vehicle',
    'cycleway', 'cycleway:left', 'cycleway:right',
    'sidewalk', 'foot', 'bicycle',
    'lts',  # Level of Traffic Stress if available
]

ALL_FORMATS = {"graphml", "pkl", "gpkg", "osm", "pbf", "geojson"}


def export_network(graph, output_dir, name, edge_tags=None, edge_tag_aggs=None,
                   formats=None):
    """Export a network graph to one or more file formats.

    Parameters
    ----------
    graph : nx.MultiDiGraph
        The network graph to export.
    output_dir : str
        Directory to write files into.
    name : str
        Base filename (without extension).
    edge_tags : list, optional
        OSM tags to include in XML export. Defaults to DEFAULT_EDGE_TAGS.
    edge_tag_aggs : list of tuples, optional
        Aggregation rules for XML export (e.g. ``[('length', 'sum')]``).
    formats : set of str, optional
        Which formats to export. Defaults to all:
        ``{"graphml", "pkl", "gpkg", "osm", "pbf", "geojson"}``.

    Returns
    -------
    dict
        Mapping of format name to output file path.
    """
    if edge_tags is None:
        edge_tags = DEFAULT_EDGE_TAGS
    if formats is None:
        formats = ALL_FORMATS

    os.makedirs(output_dir, exist_ok=True)

    # Define output file paths
    paths = {
        "graphml": os.path.join(output_dir, f'{name}.graphml'),
        "pkl":     os.path.join(output_dir, f'{name}.pkl'),
        "gpkg":    os.path.join(output_dir, f'{name}.gpkg'),
        "osm":     os.path.join(output_dir, f'{name}.osm'),
        "pbf":     os.path.join(output_dir, f'{name}.osm.pbf'),
        "geojson": os.path.join(output_dir, f'{name}.osm.geojson'),
    }

    exported = {}

    # GraphML
    if "graphml" in formats:
        ox.save_graphml(graph, filepath=paths["graphml"])
        logger.info("GRAPHML Network saved to '%s'.", paths["graphml"])
        exported["graphml"] = paths["graphml"]

    # Pickle
    if "pkl" in formats:
        with open(paths["pkl"], 'wb') as f:
            pickle.dump(graph, f)
        logger.info("PKL Network saved to '%s'.", paths["pkl"])
        exported["pkl"] = paths["pkl"]

    # GeoPackage
    if "gpkg" in formats:
        logger.info("Converting GraphML Network to GPKG Network...")
        ox.save_graph_geopackage(graph, filepath=paths["gpkg"])
        logger.info("GPKG Network saved to '%s'.", paths["gpkg"])
        exported["gpkg"] = paths["gpkg"]

    g_osm = None

    # OSM XML
    if "osm" in formats or "pbf" in formats or "geojson" in formats:
        logger.info("Creating OSM Network...")
        g_osm = _normalize_graph_for_osm_export(graph)
        save_graph_xml(
            g_osm,
            filepath=paths["osm"],
            edge_tags=edge_tags,
            edge_tag_aggs=edge_tag_aggs,
        )
        logger.info("OSM Network saved to '%s'.", paths["osm"])
        exported["osm"] = paths["osm"]

    # PBF (via pyosmium)
    if "pbf" in formats:
        _export_pbf_from_osm_xml(paths["osm"], paths["pbf"])
        logger.info("OSM PBF File saved to '%s'", paths["pbf"])
        exported["pbf"] = paths["pbf"]

    # GeoJSON (written directly from the graph edges)
    if "geojson" in formats:
        if g_osm is None:
            g_osm = _normalize_graph_for_osm_export(graph)
        _export_geojson_from_graph(g_osm, paths["geojson"])
        logger.info("OSM GEOJSON File saved to '%s'", paths["geojson"])
        exported["geojson"] = paths["geojson"]

    return exported


def _normalize_graph_for_osm_export(graph):
    """Normalize graph edge columns so they can be serialized safely."""
    nodes, edges = ox.graph_to_gdfs(graph)

    # simplify_graph(track_merged=True) can produce list-valued columns
    # (most notably osmid). File-based OSM/GeoJSON exports require scalars.
    for col in edges.columns:
        if edges[col].apply(lambda x: isinstance(x, list)).any():
            edges[col] = edges[col].apply(
                lambda x: min(x) if isinstance(x, list) else x
            )

    return ox.graph_from_gdfs(nodes, edges, graph_attrs=graph.graph)


def _import_pyosmium():
    try:
        return importlib.import_module("osmium")
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            "PBF export requires the Python package 'osmium' (pyosmium). "
            "Install it with `pip install osmium`."
        ) from exc


def _export_pbf_from_osm_xml(osm_path, pbf_path):
    """Convert an OSM XML file to OSM PBF using pyosmium."""
    osmium = _import_pyosmium()
    if os.path.exists(pbf_path):
        os.remove(pbf_path)
    with osmium.SimpleWriter(pbf_path) as writer:
        for obj in osmium.FileProcessor(osm_path):
            writer.add(obj)


def _export_geojson_from_graph(graph, geojson_path):
    """Write graph edges directly to GeoJSON without external CLI tools."""
    _, edges = ox.graph_to_gdfs(graph)
    edges = edges.reset_index()
    edges.to_file(geojson_path, driver="GeoJSON")
