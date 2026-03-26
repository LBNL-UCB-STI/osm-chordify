#!/usr/bin/env python3

"""Intersect the SF Bay land mask with the SF Bay whole-area mask."""

from __future__ import annotations

import argparse
from pathlib import Path

from _bootstrap import bootstrap_example_paths

bootstrap_example_paths(__file__)

from osm_chordify import build_area_mask_from_counties, intersect_zones_with_zones


SF_BAY_STATE_FIPS = "06"
SF_BAY_COUNTY_FIPS = ["001", "013", "041", "055", "075", "081", "085", "095", "097"]


def intersect_sfbay_masks(output_epsg: int, output_path: str):
    """Build and intersect the SF Bay land and whole-area masks."""
    output_path_obj = Path(output_path).expanduser()
    work_dir = str(output_path_obj.parent)

    land_mask = build_area_mask_from_counties(
        state_fips_code=SF_BAY_STATE_FIPS,
        county_fips_codes=SF_BAY_COUNTY_FIPS,
        year=2020,
        work_dir=work_dir,
        area_name="sfbay",
        output_epsg=output_epsg,
        include_water=False,
    )
    whole_area_mask = build_area_mask_from_counties(
        state_fips_code=SF_BAY_STATE_FIPS,
        county_fips_codes=SF_BAY_COUNTY_FIPS,
        year=2020,
        work_dir=work_dir,
        area_name="sfbay",
        output_epsg=output_epsg,
        include_water=True,
    )

    return intersect_zones_with_zones(
        zones_a=land_mask,
        zones_a_epsg=output_epsg,
        zones_b=whole_area_mask,
        zones_b_epsg=output_epsg,
        output_path=str(output_path_obj),
        output_epsg=output_epsg,
    )


def main():
    parser = argparse.ArgumentParser(description="Intersect the SF Bay land mask with the SF Bay whole-area mask.")
    parser.add_argument("--output-epsg", type=int, default=26910, help="Target EPSG for the generated masks and overlap.")
    parser.add_argument("--output-path", required=True, help="Output file path (.geojson, .gpkg, .parquet, ...).")
    args = parser.parse_args()

    output_path = str(Path(args.output_path).expanduser())
    result = intersect_sfbay_masks(output_epsg=args.output_epsg, output_path=output_path)
    print(f"Created {len(result)} SF Bay mask overlap pieces")
    print(f"Saved to {output_path}")


if __name__ == "__main__":
    main()
