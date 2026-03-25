# Usage

## Build examples

The build examples require:

- internet access for U.S. Census Bureau / `pygris` / Overpass downloads
- a U.S. Census Bureau API key, either in `CENSUS_API_KEY` or passed with `--census-api-key`

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

The validation step checks:

- topology: self-loops, isolates, connectivity
- coordinate integrity: duplicate XML-rounded coordinates and near-duplicate nodes
- edge completeness: `edge_id`, `length`, `speed_kph`, `oneway`, geometry
- export integrity: duplicate XML coordinates and dangling `<nd ref>` node references

## Validation tests

To regression-test duplicate-coordinate issues that can create downstream self-loops:

```bash
pytest tests/test_graph.py -m "not integration" -k "coordinate or integrity"
```

To run the export integration tests with printed summaries:

```bash
./.venv/bin/pytest tests/test_export.py -m integration -s
```

This runs:

- a generic real-download smoke workflow test
- a Seattle example validation test

## PBF diagnostics

To validate and diagnose a built `.osm.pbf` artifact:

```bash
./.venv/bin/python examples/diagnose_osm_pbf.py \
./output/network/seattle-cbg0-ferry-weakConn-network/seattle-cbg0-ferry-weakConn-network.osm.pbf \
  --epsg-utm 32048
```

If a sibling `.pkl`, `.graphml`, or `.osm` file exists next to the PBF, the script will also print the same build-validation metrics used by the example build validator before running the PBF-specific diagnostics.

The PBF diagnostics step uses the optional `pyrosm` dependency. Histogram output also uses `matplotlib` when available; if `matplotlib` is missing, diagnostics still run and only the histogram is skipped.

To compare two built `.osm.pbf` artifacts:

```bash
./.venv/bin/python examples/compare_osm_pbf.py \
  ./output/network/network-a/network-a.osm.pbf \
  ./output/network/network-b/network-b.osm.pbf \
  --epsg-utm-a 32048 \
  --epsg-utm-b 32048
```

When sibling `.pkl`, `.graphml`, or `.osm` files exist, the comparison also includes the build-validation metrics and reports deltas for the two artifacts.

## County-zone intersection

To intersect a road network with county polygons selected by FIPS code, use
`intersect_road_network_with_county_zones(...)`. This uses the same cached
boundary collection path as the build workflow, so county geometries are
downloaded/reused automatically under the provided `work_dir`.

Example:

```python
from osm_chordify import intersect_road_network_with_county_zones

result = intersect_road_network_with_county_zones(
    road_network="./output/network/sfbay-cbg5500-strongConn-network/sfbay-cbg5500-strongConn-network.gpkg",
    road_network_epsg=26910,
    state_fips_code="06",
    county_fips_codes=["001", "013"],
    year=2020,
    work_dir="./output/geo",
    area_name="alameda-contra-costa",
    output_path="./output/geo/county-intersections.parquet",
)
```

`road_network_epsg` must be a projected CRS. Intersections now reject
geographic CRSs such as `4326`, because the fixed length outputs are defined in
meters.

The returned output always includes:

- `zone_edge_proportion`
- `edge_link_length_m`
- `zone_link_length_m`

These three columns always refer to the current intersection step. If you
chain intersections (for example network -> counties -> grid1 -> grid2), the
latest output's top-level fixed columns are recomputed for that latest step.
Older carried attributes remain available through the prefixed `edge_...` /
`zone_...` columns.

`area_name` only controls the cached boundary filename; the actual county
selection comes from `state_fips_code` and `county_fips_codes`. County
boundaries are cached in the requested intersection CRS, so
`road_network_epsg=26910` writes an `..._epsg26910.geojson` cache instead of
always materializing WGS84 first.

## Project structure

```text
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
        data_collection.py  # U.S. Census and boundary collection
        geo.py           # Geographic utilities
        io.py            # File I/O helpers
        network.py       # Network-to-intersection mapping
```
