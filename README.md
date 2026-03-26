# osm-chordify

Download, filter, and export OSM road networks using population-density-based boundaries. Intersect the resulting geometries with polygon grids (TAZ, census tracts, ISRM cells, etc.) and map them to BEAM/MATSim simulation network.

No AI Training Use: The contents of this repository may not be used to train,
fine-tune, evaluate, benchmark, or improve machine learning or artificial
intelligence models without prior written permission. See [NOAI](./NOAI).

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
    build_area_mask_from_counties,
    build_osm_by_pop_density,
    create_osm_highway_filter,
    intersect_polygons_with_zones,
    intersect_road_polygons_with_zones,
    intersect_road_network_with_county_zones,
    intersect_road_network_with_zones,
    spatial_left_join_with_zones,
    intersect_zones_with_zones,
    map_osm_with_beam_network,
    match_road_network_geometries,
)
```

Core functions:

- `build_area_mask_from_counties`: build a fused county-based land mask or whole-area mask from FIPS codes
- `build_osm_by_pop_density`: build a multi-layer OSM network from county / tract / CBG boundaries and density filters
- `intersect_polygons_with_zones`: cascade an already-intersected polygon layer into a new zone layer
- `intersect_road_polygons_with_zones`: intersect polygon/rectangular road links with zones using area-based proportions
- `intersect_road_network_with_county_zones`: intersect a network with county polygons fetched by state/county FIPS codes
- `intersect_road_network_with_zones`: intersect network edges with zone polygons
- `spatial_left_join_with_zones`: spatially left-join any geometry layer with zone attributes while keeping unmatched rows
- `intersect_zones_with_zones`: intersect one polygon zone layer with another
- `map_osm_with_beam_network`: join BEAM network attributes onto OSM geometries
- `match_road_network_geometries`: spatially match geometries across two road networks

Intersection outputs always include:

- `zone_edge_proportion`
- `edge_link_length_m`
- `zone_link_length_m`

These three columns always describe the latest intersection step. In chained intersections, earlier prefixed `edge_...` / `zone_...` fields are carried forward, but the top-level fixed columns are recomputed for the current step.

For dense grids, `intersect_road_network_with_zones(...)` can optionally
prefilter zones to the overall road-network bounding box before exact
intersection:

- `prefilter_zones_to_network_bbox=False` disables the prefilter
- `prefilter_zones_to_network_bbox=True` keeps only zones intersecting the
  network bounding box

This reduces the zone set before exact intersection and shows a dedicated
`Filtering zones` progress bar. Zones kept by the bounding-box prefilter that
still contain no intersecting link pieces are preserved instead of being
dropped; the fixed numeric intersection columns remain null for those rows.

For rectangular/polygon road links, use
`intersect_road_polygons_with_zones(...)` instead. That variant computes
`zone_edge_proportion` from overlap area and derives `zone_link_length_m` by
applying that proportion to the full link length column. By default it first
filters zones to the overall road-network bounding box before exact
intersection and shows the same `Filtering zones` progress bar used by the
line-based workflow.

For cascading polygon workflows, use:

- `intersect_polygons_with_zones(...)`
  - preserves all existing columns from the input polygon layer
  - adds new zone attributes with `zone_*` prefixes
  - computes current-step metrics:
    - `zone_piece_proportion`
    - `piece_link_length_m`
    - `zone_piece_length_m`
- `spatial_left_join_with_zones(...)`
  - keeps all input rows
  - appends zone attributes where matched
  - leaves zone columns null where no zone intersects

Mask generation semantics:

- `include_water=False`
  - land mask
  - uses shoreline-clipped cartographic county boundaries
  - fuses counties without taking a convex hull
  - optional `buffer_m` is applied to that non-convex land geometry
- `include_water=True`
  - whole-area mask
  - uses full TIGER/Line county boundaries
  - takes the convex hull of the fused counties
  - optional `buffer_m` is applied to that convex whole-area geometry

Networks can be exported as GraphML, PKL, GPKG, OSM XML, OSM PBF, and GeoJSON.

You can also run the main workflows via CLI:

```bash
python -m osm_chordify.main --help
```

## Examples

See the [`examples/`](examples/) directory for complete, runnable scripts:

| Script | Description |
|--------|-------------|
| [`build_sfbay.py`](examples/build_sfbay.py) | Build OSM network for the SF Bay Area (9 counties) |
| [`build_seattle.py`](examples/build_seattle.py) | Build OSM network for Seattle metro (4 counties + ferry) |
| [`compare_osm_pbf.py`](examples/compare_osm_pbf.py) | Compare validation and diagnostics metrics across two built `.osm.pbf` artifacts |
| [`intersect_zones.py`](examples/intersect_zones.py) | Intersect one polygon zone layer with another |
| [`intersect_and_map_sfbay.py`](examples/intersect_and_map_sfbay.py) | Intersect OSM with polygon grid and map BEAM network |
| [`diagnose_osm_pbf.py`](examples/diagnose_osm_pbf.py) | Validate and diagnose a built `.osm.pbf` artifact |

For commands, output layout, validation behavior, and project structure, see [docs/USAGE.md](docs/USAGE.md).
