# osm-chordify

Download, filter, and export OSM road networks using population-density-based boundaries. Intersect the resulting geometries with polygon grids (TAZ, census tracts, ISRM cells, etc.) and map them to BEAM/MATSim simulation network.

## Installation

Requires Python 3.10+.

From PyPI:

```bash
pip install osm-chordify
```

With diagnostics extras:

```bash
pip install "osm-chordify[diagnostics]"
```

From GitHub:

```bash
pip install git+https://github.com/LBNL-UCB-STI/osm-chordify.git
```

For development:

```bash
git clone https://github.com/LBNL-UCB-STI/osm-chordify.git
cd osm-chordify
pip install -e ".[dev]"
```

For diagnostics in a local clone:

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

You can also pass the key directly to the example scripts with `--census-api-key`.

## API

```python
from osm_chordify import (
    build_osm_by_pop_density,
    create_osm_highway_filter,
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
    work_dir="./output",
    osm_config=osm_config,   # osmnx settings, graph layers, filters
    area_config=area_config,  # name, state/county FIPS, census year
    geo_config=geo_config,    # UTM EPSG and projection settings
)
```

The `osm_config` dict defines one or more `graph_layers`, each with a `custom_filter` (Overpass QL highway filter), `geo_level`, optional `min_density_per_km2`, and `buffer_zone_in_meters`. Layers are downloaded independently and merged into a single network that can be exported as GraphML, PKL, GPKG, OSM XML, OSM PBF, and GeoJSON.

### `intersect_road_network_with_zones`

Intersect road-network edges with zone polygons. Computes the proportion of each edge that sits within each zone. All attributes from both inputs are carried through, prefixed with `edge_` and `zone_`. Use `proportional_cols` to specify edge columns that should be scaled by the intersection proportion. Both `road_network` and `zones` accept either a file path (GPKG, GeoJSON, Shapefile) or a GeoDataFrame.

```python
# From files — scale edge_length and vmt by proportion
intersect_road_network_with_zones(
    road_network="sfbay.gpkg",
    road_network_epsg=26910,
    zones="isrm_polygon.shp",
    zones_epsg=26910,
    proportional_cols=["edge_length", "vmt"],
    output_path="intersection.geojson",
    output_epsg=26910,
)

# From GeoDataFrames
result_gdf = intersect_road_network_with_zones(
    road_network=edges_gdf,
    road_network_epsg=26910,
    zones=polygons_gdf,
    zones_epsg=26910,
    proportional_cols="edge_length",
)
```

### `map_osm_with_beam_network`

Join a BEAM network CSV to a zone-network intersection result on the shared OSM ID. All columns from both inputs are included in the output.

```python
map_osm_with_beam_network(
    osm_path="sfbay.gpkg",
    network_path="network.csv.gz",
    network_osm_id_col="attributeOrigId",
    output_path="mapping.geojson",
)
```

### `match_road_network_geometries`

Spatially match link geometries between two road networks. Accepts file paths or GeoDataFrames. Use `matching="strict"` for high-overlap matches only, or `matching="flexible"` (default) for partial and nearby matches.

```python
result_gdf = match_road_network_geometries(
    network_a="beam_osm.geojson",
    network_a_epsg=4326,
    network_b="hpms_network.shp",
    network_b_epsg=26910,
    matching="flexible",
    output_path="matched.geojson",
    output_epsg=26910,
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

### Running example builds

The build examples require:

- internet access for Census / pygris / Overpass downloads
- a Census API key, either in `CENSUS_API_KEY` or passed with `--census-api-key`

Seattle example:

```bash
./.venv/bin/python examples/build_seattle.py
```

With an explicit API key and output directory:

```bash
./.venv/bin/python examples/build_seattle.py \
  --census-api-key "$CENSUS_API_KEY" \
  --output-dir ./output
```

SF Bay example:

```bash
./.venv/bin/python examples/build_sfbay.py --output-dir ./output
```

By default, the example scripts write under `./output` relative to the directory you run them from, and they create the needed structure automatically:

- `output/geo/`
- `output/network/`

After the build completes, the scripts automatically run graph/XML validation unless you disable it:

```bash
./.venv/bin/python examples/build_seattle.py --skip-validation
```

The validation step checks topology, duplicate/near-duplicate coordinates, missing edge attributes, and XML node-reference integrity, then prints a readable summary of what it inspected.

To regression-test the resulting network against duplicate-node coordinate issues that can create R5 self-loops, run:

```bash
pytest tests/test_graph.py -m "not integration" -k "coordinate or integrity"
```

## Project structure

```
src/osm_chordify/
    __init__.py          # Public API
    main.py              # Orchestration pipeline
    osm/
        __init__.py      # OSM subpackage
        analyze.py       # PBF analysis helpers
        diagnostics.py   # Network validation
        export.py        # Network export (GraphML, GPKG, OSM XML, PBF, GeoJSON)
        graph.py         # OSM network download and preparation
        intersect.py     # Edge-polygon intersection
        simplify.py      # Graph simplification aggregation helpers
        tags.py          # OSM tag parsing
        xml.py           # OSM XML serialization
    utils/
        __init__.py      # Utility subpackage
        data_collection.py  # Census data and boundaries
        geo.py           # Geographic utilities
        io.py            # File I/O helpers
        network.py       # Network-to-intersection mapping
```
