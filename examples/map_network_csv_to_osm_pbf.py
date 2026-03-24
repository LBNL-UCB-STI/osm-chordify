"""Map a network CSV/CSV.GZ to a built OSM PBF and save a spatial join."""

import argparse
import sys
from pathlib import Path

EXAMPLES_DIR = Path(__file__).resolve().parent
if str(EXAMPLES_DIR) not in sys.path:
    sys.path.insert(0, str(EXAMPLES_DIR))

from _bootstrap import bootstrap_example_paths

bootstrap_example_paths(__file__)

from osm_chordify import map_osm_with_beam_network


def map_network_csv_to_osm_pbf(
    network_csv_path,
    osm_pbf_path,
    output_path,
    network_osm_id_col="attributeOrigId",
):
    """Map a network CSV/CSV.GZ to an OSM PBF and save a spatial join."""
    if not str(osm_pbf_path).lower().endswith(".pbf"):
        raise ValueError(
            f"osm_pbf_path must point to a .pbf/.osm.pbf file, got: {osm_pbf_path}"
        )
    if Path(output_path).suffix.lower() not in {".parquet", ".gpkg", ".geojson"}:
        raise ValueError("output_path must end with .parquet, .gpkg, or .geojson")

    return map_osm_with_beam_network(
        osm_path=osm_pbf_path,
        network_path=network_csv_path,
        network_osm_id_col=network_osm_id_col,
        output_path=output_path,
    )


def _parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("network_csv_path")
    parser.add_argument("osm_pbf_path")
    parser.add_argument("output_path")
    parser.add_argument("--network-osm-id-col", default="attributeOrigId")
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    map_network_csv_to_osm_pbf(
        network_csv_path=args.network_csv_path,
        osm_pbf_path=args.osm_pbf_path,
        output_path=args.output_path,
        network_osm_id_col=args.network_osm_id_col,
    )
