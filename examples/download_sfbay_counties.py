#!/usr/bin/env python3

"""Download SF Bay Area county boundaries using the shared boundary helper."""

from __future__ import annotations

import argparse
from pathlib import Path

from _bootstrap import bootstrap_example_paths

bootstrap_example_paths(__file__)

from osm_chordify.utils.data_collection import collect_geographic_boundaries
from osm_chordify.utils.io import save_geodataframe


SF_BAY_STATE_FIPS = "06"
SF_BAY_COUNTY_FIPS = ["001", "013", "041", "055", "075", "081", "085", "095", "097"]


def download_sfbay_counties(work_dir: str, output_path: str | None = None, epsg: int = 4326, year: int = 2020):
    """Download and optionally save SF Bay Area county boundaries."""
    counties = collect_geographic_boundaries(
        state_fips_code=SF_BAY_STATE_FIPS,
        county_fips_codes=SF_BAY_COUNTY_FIPS,
        year=year,
        area_name="sfbay",
        geo_level="county",
        work_dir=work_dir,
        target_epsg=epsg,
    )

    if output_path:
        save_geodataframe(counties, output_path)
    return counties


def main():
    parser = argparse.ArgumentParser(description="Download SF Bay Area county boundaries.")
    parser.add_argument("--work-dir", default="output/geo", help="Directory used for boundary cache files.")
    parser.add_argument("--output-path", help="Optional output file path (.geojson, .gpkg, .parquet, ...).")
    parser.add_argument("--epsg", type=int, default=4326, help="Target EPSG for the downloaded county boundaries.")
    parser.add_argument("--year", type=int, default=2020, help="Boundary reference year.")
    args = parser.parse_args()

    work_dir = str(Path(args.work_dir).expanduser())
    output_path = str(Path(args.output_path).expanduser()) if args.output_path else None
    counties = download_sfbay_counties(
        work_dir=work_dir,
        output_path=output_path,
        epsg=args.epsg,
        year=args.year,
    )
    print(f"Downloaded {len(counties)} SF Bay counties in EPSG:{args.epsg}")
    if output_path:
        print(f"Saved to {output_path}")


if __name__ == "__main__":
    main()
