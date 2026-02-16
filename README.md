# osm-chordify

Download, filter, and export OSM road networks using population-density-based boundaries. Intersect the resulting geometries with polygon grids (TAZ, census tracts, ISRM cells, etc.) and map them to simulation networks like BEAM.

## Installation

Requires Python 3.10+.

```bash
pip install git+https://github.com/LBNL-UCB-STI/osm-chordify.git
```

For development:

```bash
git clone https://github.com/LBNL-UCB-STI/osm-chordify.git
cd osm-chordify
pip install -e ".[dev]"
```

For diagnostics (connectivity checks, link-length histograms):

```bash
pip install -e ".[diagnostics]"
```

### Census API key

`build_osm_by_pop_density` uses the Census Bureau API to fetch population data for density-based filtering. Set your API key as an environment variable:

```bash
export CENSUS_API_KEY="your_key_here"
```

To make it persistent, add the line above to your `~/.bashrc`, `~/.zshrc`, or `.env` file.

Get a free key at https://api.census.gov/data/key_signup.html.

## API

```python
from osm_chordify import (
    build_osm_by_pop_density,
    intersect_road_network_with_zones,
    map_osm_with_beam_network,
    match_road_network_geometries,
    diagnose_osm,
)
```

### `build_osm_by_pop_density`

Download and build a multi-layer OSM network. Each layer can target a different geographic level (county, census block group, tract) and apply a population density threshold to limit residential road downloads to urban areas.

```python
build_osm_by_pop_density(
    work_dir="~/Workspace/Simulation/sfbay",
    osm_config=osm_config,   # osmnx settings, graph layers, filters
    area_config=area_config,  # name, state/county FIPS, census year
    geo_config=geo_config,    # UTM EPSG, TAZ shapefile path
)
```

The `osm_config` dict defines one or more `graph_layers`, each with a `custom_filter` (Overpass QL highway filter), `geo_level`, optional `min_density_per_km2`, and `buffer_zone_in_meters`. Layers are downloaded independently and merged into a single network exported as GeoJSON + GPKG.

### `intersect_road_network_with_zones`

Intersect road-network edges with zone polygons. Computes the proportion of each edge within each zone and the corresponding length in meters. All attributes from both inputs are carried through, prefixed with `edge_` and `zone_` to avoid collisions. Both `road_network` and `zones` accept either a file path (GPKG, GeoJSON, Shapefile) or a GeoDataFrame.

```python
# From files
intersect_road_network_with_zones(
    road_network="sfbay.gpkg",
    zones="isrm_polygon.shp",
    epsg_utm=26910,
    output_path="intersection.geojson",
)

# From GeoDataFrames
result_gdf = intersect_road_network_with_zones(
    road_network=edges_gdf,
    zones=polygons_gdf,
    epsg_utm=26910,
)
```

### `map_osm_with_beam_network`

Join a BEAM network CSV to a zone-network intersection result on the shared OSM ID. All columns from both inputs are included in the output.

```python
map_osm_with_beam_network(
    network_path="network.csv.gz",
    intersection_path="intersection.geojson",
    network_osm_id_col="attributeOrigId",
    output_path="mapping.geojson",
)
```

### `match_road_network_geometries`

Spatially match link geometries between two road networks. Accepts file paths or GeoDataFrames. Use `matching="strict"` for high-overlap matches only, or `matching="flexible"` (default) for partial and nearby matches.

```python
result_gdf = match_road_network_geometries(
    network_a="osm_network.gpkg",
    network_b="beam_network.geojson",
    epsg_utm=26910,
    matching="flexible",
    output_path="matched.geojson",
)
```

> **Note:** This function is a placeholder and not yet implemented.

### `diagnose_osm`

Run connectivity and link-length diagnostics on a `.osm.pbf` file. Reports connected components, short/long link anomalies, and saves a histogram. Requires the `diagnostics` extra.

```python
diagnose_osm("sfbay.osm.pbf", epsg_utm=26910)
```

## Examples

See the [`examples/`](examples/) directory for complete, runnable scripts:

| Script | Description |
|--------|-------------|
| [`build_sfbay.py`](examples/build_sfbay.py) | Build OSM network for the SF Bay Area (9 counties) |
| [`build_seattle.py`](examples/build_seattle.py) | Build OSM network for Seattle metro (4 counties + ferry) |
| [`intersect_and_map_sfbay.py`](examples/intersect_and_map_sfbay.py) | Intersect OSM with polygon grid and map BEAM network |
| [`diagnose_sfbay.py`](examples/diagnose_sfbay.py) | Run diagnostics on a downloaded PBF |

## Project structure

```
src/osm_chordify/
    __init__.py          # Public API
    main.py              # Orchestration pipeline
    osm/
        graph.py         # OSM network download and preparation
        export.py        # Network export (GeoJSON, GPKG)
        intersect.py     # Edge-polygon intersection
        analyze.py       # PBF analysis (osmium)
        diagnostics.py   # Network validation
        tags.py          # OSM tag parsing
    utils/
        geo.py           # Geographic utilities
        data_collection.py  # Census data and boundaries
        network.py       # Network-to-intersection mapping
        io.py            # File I/O helpers
```
