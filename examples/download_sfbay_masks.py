#!/usr/bin/env python3

"""Generate both SF Bay area masks: land-only and whole-area."""

from __future__ import annotations

import argparse
from pathlib import Path

from _bootstrap import bootstrap_example_paths

bootstrap_example_paths(__file__)

from osm_chordify import build_area_mask_from_counties


SF_BAY_STATE_FIPS = "06"
SF_BAY_COUNTY_FIPS = ["001", "013", "041", "055", "075", "081", "085", "095", "097"]


def generate_sfbay_masks(work_dir: str, land_output_path: str, whole_area_output_path: str, epsg: int = 4326, year: int = 2020, buffer_m: float = 0.0):
    """Generate and save both land-only and whole-area masks for SF Bay."""
    land_mask = build_area_mask_from_counties(
        state_fips_code=SF_BAY_STATE_FIPS,
        county_fips_codes=SF_BAY_COUNTY_FIPS,
        year=year,
        work_dir=work_dir,
        area_name="sfbay",
        output_epsg=epsg,
        include_water=False,
        buffer_m=buffer_m,
        output_path=land_output_path,
    )
    whole_area_mask = build_area_mask_from_counties(
        state_fips_code=SF_BAY_STATE_FIPS,
        county_fips_codes=SF_BAY_COUNTY_FIPS,
        year=year,
        work_dir=work_dir,
        area_name="sfbay",
        output_epsg=epsg,
        include_water=True,
        buffer_m=buffer_m,
        output_path=whole_area_output_path,
    )
    return land_mask, whole_area_mask


def main():
    parser = argparse.ArgumentParser(description="Generate land-only and whole-area masks for the SF Bay counties.")
    parser.add_argument("--work-dir", default="output/geo", help="Directory used for boundary cache files.")
    parser.add_argument("--land-output-path", required=True, help="Output file path for the land-only mask.")
    parser.add_argument("--whole-area-output-path", required=True, help="Output file path for the whole-area mask.")
    parser.add_argument("--epsg", type=int, default=4326, help="Target EPSG for both output masks.")
    parser.add_argument("--year", type=int, default=2020, help="Boundary reference year.")
    parser.add_argument("--buffer-m", type=float, default=0.0, help="Meter buffer applied to both generated masks.")
    args = parser.parse_args()

    work_dir = str(Path(args.work_dir).expanduser())
    land_output_path = str(Path(args.land_output_path).expanduser())
    whole_area_output_path = str(Path(args.whole_area_output_path).expanduser())
    generate_sfbay_masks(
        work_dir=work_dir,
        land_output_path=land_output_path,
        whole_area_output_path=whole_area_output_path,
        epsg=args.epsg,
        year=args.year,
        buffer_m=args.buffer_m,
    )
    print(f"Saved land-only mask to {land_output_path}")
    print(f"Saved whole-area mask to {whole_area_output_path}")


if __name__ == "__main__":
    main()
