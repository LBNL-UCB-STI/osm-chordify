"""Geo/projection utilities."""

import logging

import geopandas as gpd
import networkx as nx
import osmnx as ox
import pyproj
from osmnx import settings

logger = logging.getLogger(__name__)


def name_osm_network(area_name, graph_layers, strongly_connected):
    """
    Generate an OSM product name from study-area parameters.

    Format: ``{area}[-{geo}{density}][-ferry]-{conn}Conn-network``

    Parameters
    ----------
    area_name : str
        Study area name (e.g. ``"sfbay"``).
    graph_layers : dict
        Graph layer definitions (must include ``"residential"`` key if used).
    strongly_connected : bool
        Whether the network uses strongly connected components.

    Returns
    -------
    str
    """
    if "residential" in graph_layers:
        density = str(graph_layers["residential"]["min_density_per_km2"])
        res_part = f"-{graph_layers['residential']['geo_level']}{density}"
    else:
        res_part = ""

    conn = "strong" if strongly_connected else "weak"
    ferry = "-ferry" if "ferry" in graph_layers else ""

    return f"{area_name}{res_part}{ferry}-{conn}Conn-network"


def create_osm_highway_filter(highway_types):
    """Convert a list of highway types to an OSM custom filter string.

    Parameters
    ----------
    highway_types : list[str]
        Highway type strings (e.g. ``["motorway", "trunk", "primary"]``).

    Returns
    -------
    str
        Filter in the format ``'["highway"~"type1|type2|..."]'``.
    """
    return f'["highway"~"{"|".join(highway_types)}"]'


def meters_to_degrees(lon, lat, utm_epsg, buffer_meters):
    """
    Calculate the equivalent buffer distance in degrees for a given buffer in meters,
    using a specified UTM projection for better precision.

    Parameters
    ----------
    lon : float
        Longitude coordinate (x) in WGS84
    lat : float
        Latitude coordinate (y) in WGS84
    utm_epsg : int
        The EPSG code for the UTM coordinate reference system (e.g., 26910 for UTM Zone 10N)
    buffer_meters : float
        Buffer distance in meters

    Returns
    -------
    float
        Equivalent buffer distance in degrees
    """
    # Create UTM CRS from EPSG code
    utm_crs = f"EPSG:{utm_epsg}"

    # Create transformers
    wgs84_to_utm = pyproj.Transformer.from_crs("EPSG:4326", utm_crs, always_xy=True)
    utm_to_wgs84 = pyproj.Transformer.from_crs(utm_crs, "EPSG:4326", always_xy=True)

    # Convert coordinates to UTM
    x_utm, y_utm = wgs84_to_utm.transform(lon, lat)

    # Calculate points at buffer distance in cardinal directions
    east_utm = (x_utm + buffer_meters, y_utm)
    north_utm = (x_utm, y_utm + buffer_meters)

    # Convert buffered points back to WGS84
    east_lon, east_lat = utm_to_wgs84.transform(*east_utm)
    north_lon, north_lat = utm_to_wgs84.transform(*north_utm)

    # Calculate degree differences
    lon_diff = abs(east_lon - lon)  # East-West difference (longitude)
    lat_diff = abs(north_lat - lat)  # North-South difference (latitude)

    # Return the average as an approximation
    # You could also return both separately if you need different buffers for lat/lon
    return (lon_diff + lat_diff) / 2


def to_convex_hull(input_data, utm_epsg, buffer_in_meters):
    """
    Create a buffered convex hull from input data.

    Parameters
    ----------
    input_data : GeoDataFrame, GeoSeries, or Shapely geometry
        The input geographic data
    utm_epsg : int
        EPSG code for the UTM projection to use for accurate distance calculations
    buffer_in_meters : float
        Buffer distance in meters

    Returns
    -------
    Shapely geometry
        The buffered convex hull
    """
    # Handle different input types
    if isinstance(input_data, gpd.GeoDataFrame):
        # GeoDataFrame: get the convex hull of all geometries
        convex_hull = input_data.geometry.unary_union.convex_hull
    elif isinstance(input_data, gpd.GeoSeries):
        # GeoSeries: get the convex hull of all geometries
        convex_hull = input_data.unary_union.convex_hull
    elif hasattr(input_data, 'geom_type'):
        # Shapely geometry: get its convex hull
        convex_hull = input_data.convex_hull
    else:
        raise TypeError("Input must be a GeoDataFrame, GeoSeries, or Shapely geometry")

    # Get centroid
    lon = convex_hull.centroid.x
    lat = convex_hull.centroid.y

    # Convert buffer distance
    buffer_in_degrees = meters_to_degrees(lon, lat, utm_epsg, buffer_in_meters)

    # Buffer in degrees
    buffered_convex_hull = convex_hull.buffer(buffer_in_degrees)

    return buffered_convex_hull


def project_graph(G: nx.MultiDiGraph, to_crs=None, to_latlong=False) -> nx.MultiDiGraph:
    """
    Project a graph from its current CRS to another.

    If `to_latlong` is True, this projects the graph to the coordinate
    reference system defined by `settings.default_crs`. Otherwise it projects
    it to the CRS defined by `to_crs`. If `to_crs` is `None`, it projects it
    to the CRS of an appropriate UTM zone given `geometry`'s bounds.

    Parameters
    ----------
    G
        The graph to be projected.
    to_crs
        If None, project to an appropriate UTM zone. Otherwise project to
        this CRS.
    to_latlong
        If True, project to `settings.default_crs` and ignore `to_crs`.

    Returns
    -------
    G_proj
        The projected graph.
    """
    if to_latlong:
        to_crs = settings.default_crs

    # STEP 1: PROJECT THE NODES
    gdf_nodes = ox.convert.graph_to_gdfs(G, edges=False)

    # project the nodes GeoDataFrame and extract the projected x/y values
    gdf_nodes_proj = ox.projection.project_gdf(gdf_nodes, to_crs=to_crs)
    gdf_nodes_proj["x"] = gdf_nodes_proj["geometry"].x
    gdf_nodes_proj["y"] = gdf_nodes_proj["geometry"].y
    to_crs = gdf_nodes_proj.crs

    # STEP 2: PROJECT THE EDGES
    # Always get edges with geometry, regardless of whether the graph is simplified
    gdf_edges = ox.convert.graph_to_gdfs(G, nodes=False, fill_edge_geometry=True)

    # If edges don't have a CRS but do have geometry, assign the source CRS
    if gdf_edges.crs is None and not gdf_edges.empty and 'geometry' in gdf_edges.columns:
        # If we're unsure about the source CRS, use what we know from the nodes
        source_crs = G.graph.get('crs', gdf_nodes.crs)
        if source_crs is not None:
            gdf_edges.crs = source_crs
            logger.debug("Setting edge CRS to %s before projection", source_crs)

    # Project the edges
    gdf_edges_proj = ox.projection.project_gdf(gdf_edges, to_crs=to_crs)

    # Debug output to verify projection worked
    if not gdf_edges_proj.empty and 'geometry' in gdf_edges_proj.columns:
        sample_geom = gdf_edges_proj.iloc[0]['geometry']
        if sample_geom is not None:
            logger.debug("Sample edge coordinate after projection: %s", next(iter(sample_geom.coords)))

    # STEP 3: REBUILD GRAPH
    # turn projected node/edge gdfs into a graph and update its CRS attribute
    G_proj = ox.convert.graph_from_gdfs(gdf_nodes_proj, gdf_edges_proj, graph_attrs=G.graph)
    G_proj.graph["crs"] = to_crs

    logger.info("Projected graph with %d nodes and %d edges", len(G), len(G.edges))

    # Final verification
    nodes_check, edges_check = ox.convert.graph_to_gdfs(G_proj)
    logger.debug("Verified: Nodes CRS: %s, Edges CRS: %s", nodes_check.crs, edges_check.crs)

    return G_proj
