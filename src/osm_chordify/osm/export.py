"""Export a network graph to one or more file formats."""

import os
import pickle
import subprocess

import osmnx as ox

from osm_chordify.osm.xml import save_graph_xml

# Absolute path to the package directory (for locating _osm_conf.ini)
_PACKAGE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

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
                   formats=None, osm_conf_path=None):
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
    osm_conf_path : str, optional
        Path to ``_osm_conf.ini`` for ogr2ogr. Defaults to the one shipped
        with the package.

    Returns
    -------
    dict
        Mapping of format name to output file path.
    """
    if edge_tags is None:
        edge_tags = DEFAULT_EDGE_TAGS
    if formats is None:
        formats = ALL_FORMATS
    if osm_conf_path is None:
        osm_conf_path = os.path.join(_PACKAGE_DIR, '_osm_conf.ini')

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
        print(f"GRAPHML Network saved to '{paths['graphml']}'.")
        exported["graphml"] = paths["graphml"]

    # Pickle
    if "pkl" in formats:
        with open(paths["pkl"], 'wb') as f:
            pickle.dump(graph, f)
        print(f"PKL Network saved to '{paths['pkl']}'.")
        exported["pkl"] = paths["pkl"]

    # GeoPackage
    if "gpkg" in formats:
        print("Converting GraphML Network to GPKG Network...")
        ox.save_graph_geopackage(graph, filepath=paths["gpkg"])
        print(f"GPKG Network saved to '{paths['gpkg']}'.")
        exported["gpkg"] = paths["gpkg"]

    # OSM XML
    if "osm" in formats or "pbf" in formats or "geojson" in formats:
        print("Creating OSM Network...")
        nodes, edges = ox.graph_to_gdfs(graph)
        g_osm = ox.graph_from_gdfs(nodes, edges, graph_attrs=graph.graph)
        save_graph_xml(
            g_osm,
            filepath=paths["osm"],
            edge_tags=edge_tags,
            edge_tag_aggs=edge_tag_aggs,
        )
        print(f"OSM Network saved to '{paths['osm']}'.")
        exported["osm"] = paths["osm"]

    # PBF (requires osmium CLI and an OSM file)
    if "pbf" in formats:
        cmd = (
            f"osmium cat {paths['osm']} -o - --output-format pbf,compression=zlib "
            f"| osmium sort -F pbf - -o {paths['pbf']} --overwrite"
        )
        subprocess.run(cmd, shell=True, check=True)
        print(f"OSM PBF File saved to '{paths['pbf']}'")
        exported["pbf"] = paths["pbf"]

    # GeoJSON (requires ogr2ogr and a PBF file)
    if "geojson" in formats:
        cmd = (
            f'ogr2ogr -f GeoJSON "{paths["geojson"]}" "{paths["pbf"]}" lines '
            f'--config OSM_CONFIG_FILE "{osm_conf_path}"'
        )
        subprocess.run(cmd, shell=True, check=True)
        print(f"OSM GEOJSON File saved to '{paths['geojson']}'")
        exported["geojson"] = paths["geojson"]

    return exported
