# osm-chordify

Download, filter, and export OSM road networks using population-density-based boundaries. Intersect the resulting geometries with polygon grids (TAZ, census tracts, ISRM cells, etc.) and map them to BEAM/MATSim simulation network.

## Installation

Requires Python 3.10+.

From PyPI:

```bash
pip install osm-chordify
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

## U.S. Census API key

`build_osm_by_pop_density` uses the U.S. Census Bureau API to fetch population data for density-based filtering. Set your API key as an environment variable:

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
)
```

Core functions:

- `build_osm_by_pop_density`: build a multi-layer OSM network from county / tract / CBG boundaries and density filters
- `intersect_road_network_with_zones`: intersect network edges with zone polygons
- `map_osm_with_beam_network`: join BEAM network attributes onto OSM geometries
- `match_road_network_geometries`: spatially match geometries across two road networks

Networks can be exported as GraphML, PKL, GPKG, OSM XML, OSM PBF, and GeoJSON.

## Examples

See the [`examples/`](examples/) directory for complete, runnable scripts:

| Script | Description |
|--------|-------------|
| [`build_sfbay.py`](examples/build_sfbay.py) | Build OSM network for the SF Bay Area (9 counties) |
| [`build_seattle.py`](examples/build_seattle.py) | Build OSM network for Seattle metro (4 counties + ferry) |
| [`intersect_and_map_sfbay.py`](examples/intersect_and_map_sfbay.py) | Intersect OSM with polygon grid and map BEAM network |
| [`diagnose_sfbay.py`](examples/diagnose_sfbay.py) | Run diagnostics on a downloaded PBF |

For commands, output layout, validation behavior, and project structure, see [docs/USAGE.md](docs/USAGE.md).
