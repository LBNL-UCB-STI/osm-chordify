# Changelog

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
