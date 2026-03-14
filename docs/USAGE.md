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
