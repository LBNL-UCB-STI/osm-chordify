# Changelog

## 0.2.1

- Added area-agnostic PBF diagnostics and comparison workflows, including `diagnose_osm_pbf.py`, `compare_osm_pbf.py`, and corresponding `python -m osm_chordify.main` CLI commands.
- Added an example workflow for mapping `network.csv(.gz)` to `.osm.pbf` and saving a GeoJSON join, and hardened PBF reading for GDAL/pyogrio indexing quirks.
- Made the package CLI cleaner by adding lazy package exports and deferring heavy imports so `python -m osm_chordify.main --help` does not emit runtime/import noise.
- Reworked the example scripts to bootstrap the local `src/` tree automatically, so they run directly from a repo checkout without needing `PYTHONPATH=src`.
- Split the short project overview from the longer operational instructions by keeping the README compact and moving detailed how-to material into `docs/USAGE.md`.
- Simplified SF Bay and Seattle example wiring, aligned `main.py` with generic CLI-driven workflows, and improved example/test coverage for the new commands.

## 0.2.0

- Replaced external `osmium` / `ogr2ogr` export subprocesses with Python-backed `.osm.pbf` and GeoJSON export, and added the Python `osmium` dependency to package metadata.
- Added multi-endpoint Overpass fallback for live downloads instead of relying on a single Overpass host.
- Hardened graph preparation for current OSMnx behavior, including the second consolidation pass used to eliminate near-duplicate nodes that can collapse during XML/PBF export.
- Added stronger graph and export validation for built networks, including duplicate-coordinate checks, near-duplicate-node checks, edge attribute completeness checks, and XML node-reference integrity checks.
- Updated example build scripts to use `--output-dir`, `--census-api-key`, optional `--skip-validation`, and a default `./output` workspace with auto-created `geo/` and `network/` subdirectories.
- Aligned the example and integration validation output so Seattle validation and build-time validation use the same shared summary logic, including highway-type counts.
- Added integration coverage for live download plus PBF export/parse workflow, with a Seattle example validation test and a smaller real-download smoke test.
- Updated README guidance for example usage, output layout, validation behavior, and current export capabilities.

## 0.1.1

- Replaced `print(...)` statements with structured `logging` across OSM processing and export modules.
- Standardized logger messages to lazy formatting style for consistency and library-friendly output.
- Improved diagnostics and graph/export logging behavior to avoid direct stdout writes in library code.

## 0.1.0

- Initial scaffold
