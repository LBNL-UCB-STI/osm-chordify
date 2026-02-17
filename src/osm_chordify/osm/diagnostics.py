"""Validation & QA helpers."""

import csv
import logging
import os
import subprocess

import osmnx as ox
import pandas as pd

logger = logging.getLogger(__name__)


def find_long_tags_in_gdf(gdf, element_type="elements"):
    """
    Find columns and combinations of attributes that exceed 250 characters in a GeoDataFrame.

    Parameters
    ----------
    gdf : GeoDataFrame
        The input GeoDataFrame (can be either nodes or edges)
    element_type : str, optional
        The type of elements being analyzed ("nodes" or "edges") for output messages

    Returns
    -------
    tuple
        (long_tags, long_comb_tags) where:
        - long_tags: dict of individual columns with values >= 250 characters
        - long_comb_tags: dict of rows with combined attribute length >= 250 characters
    """
    logger.info("Analyzing %s...", element_type)

    # Find individual columns with values longer than 250 characters
    long_tags = {}
    for column in gdf.columns:
        # Convert all values to strings and check their lengths
        max_length = gdf[column].astype(str).str.len().max()
        if max_length >= 250:
            long_tags[column] = max_length

    # Print results for individual columns
    if long_tags:
        logger.info("Individual %s columns with values >= 250 characters:", element_type)
        for column, length in long_tags.items():
            logger.info("Column '%s': max length = %d characters", column, length)
            # Print an example of a long value
            long_value_idx = gdf[column].astype(str).str.len().idxmax()
            logger.info("Example long value: %s", gdf[column].iloc[long_value_idx])
    else:
        logger.info("No individual %s columns found with values >= 250 characters", element_type)

    # Find combinations of attributes that exceed 250 characters
    logger.info("Checking %s attribute combinations...", element_type)
    # Get all rows where any combination of attributes might be long
    long_comb_tags = {}
    for idx, row in gdf.iterrows():
        comb_length = 0
        contributing_cols = []

        for col in gdf.columns:
            value = str(row[col])
            if len(value) > 0 and value.lower() != 'nan':  # Skip empty or NaN values
                value_length = len(value)
                comb_length += value_length
                if value_length > 0:  # Only add if the value has length
                    contributing_cols.append({
                        'column': col,
                        'length': value_length,
                        'value': value
                    })

        if comb_length >= 250:
            long_comb_tags[idx] = {
                'total_length': comb_length,
                'contributing_columns': contributing_cols
            }

    # Print results for combinations
    if long_comb_tags:
        logger.info("%s rows with combined attribute length >= 250 characters:", element_type.capitalize())
        for idx, info in long_comb_tags.items():
            logger.info("Row %s:", idx)
            logger.info("Total combined length: %d characters", info['total_length'])
            logger.info("Contributing columns:")
            for col_info in info['contributing_columns']:
                logger.info("- %s: length=%d chars", col_info['column'], col_info['length'])
                if col_info['length'] > 50:  # Show value only if it's significantly long
                    logger.info("  Value: %s...", col_info['value'][:50])  # Show first 50 chars
    else:
        logger.info("No combinations of %s attributes found exceeding 250 characters", element_type)

    return long_tags, long_comb_tags


def check_duplicate_edge_ids(edges_gdf, id_column='edge_id'):
    """
    Check for duplicate edge IDs in an OSMnx edges GeoDataFrame.

    Parameters
    ----------
    edges_gdf : GeoDataFrame
        The edges GeoDataFrame from ox.graph_to_gdfs()
    id_column : str, default 'edge_id'
        The column name containing the edge IDs to check

    Returns
    -------
    tuple
        (has_duplicates, duplicate_info) where:
        - has_duplicates: Boolean indicating if duplicates were found
        - duplicate_info: DataFrame containing the duplicate IDs and their counts
    """
    # Count occurrences of each edge_id
    id_counts = edges_gdf[id_column].value_counts()

    # Filter to only those with count > 1 (duplicates)
    duplicates = id_counts[id_counts > 1]

    if len(duplicates) > 0:
        # Create a DataFrame with duplicate IDs and their counts
        duplicate_info = duplicates.reset_index()
        duplicate_info.columns = ['edge_id', 'count']

        # Get examples of each duplicate
        examples = []
        for dup_id in duplicate_info['edge_id']:
            # Get the first few examples of this duplicate ID
            example_edges = edges_gdf[edges_gdf[id_column] == dup_id].head(3)
            examples.append(example_edges)

        if examples:
            # Concatenate all example edges into one DataFrame
            examples_df = pd.concat(examples)
            duplicate_info = (duplicate_info, examples_df)

        logger.info("Found %d duplicate edge IDs out of %d total edges", len(duplicates), len(edges_gdf))
        return True, duplicate_info
    else:
        logger.info("No duplicate edge IDs found in %d edges", len(edges_gdf))
        return False, None


def check_invalid_coordinates(graph):
    """
    Check for invalid coordinates in the graph nodes.

    Parameters
    ----------
    graph : networkx.MultiDiGraph
        The graph to check

    Returns
    -------
    tuple
        (has_invalid, invalid_nodes) where:
        - has_invalid: boolean indicating if any invalid coordinates were found
        - invalid_nodes: list of node IDs with invalid coordinates
    """
    nodes, _ = ox.graph_to_gdfs(graph)

    # Check for NaN, infinite, or out-of-range coordinates
    invalid_x = ~nodes['x'].between(-180, 180) | nodes['x'].isna() | nodes['x'].abs().eq(float('inf'))
    invalid_y = ~nodes['y'].between(-90, 90) | nodes['y'].isna() | nodes['y'].abs().eq(float('inf'))

    # Combine invalid x or y
    invalid_nodes = nodes[invalid_x | invalid_y]

    if len(invalid_nodes) > 0:
        logger.warning("Found %d nodes with invalid coordinates:", len(invalid_nodes))
        for idx, node in invalid_nodes.iterrows():
            logger.warning("  Node ID: %s, x: %s, y: %s", idx, node['x'], node['y'])
        return True, invalid_nodes.index.tolist()

    return False, []


def scan_network_directories_for_ways(directory):
    def calculate_ways(osm_file):
        try:
            # Use osmium to get file info with summary
            result = subprocess.run(['osmium', 'fileinfo', '-e', osm_file],
                                    capture_output=True, text=True)
            # Initialize ways_count variable
            ways_count = 0

            # Extract the number of ways from the output
            for line in result.stdout.splitlines():
                if "Number of ways" in line:
                    ways_count = line.split(":")[1].strip()  # Get the number of ways
                    break  # Stop after finding the count

            return ways_count  # Return the number of ways
        except Exception as e:
            logger.error("Error processing %s: %s", osm_file, e)
        return 0

    output_file = os.path.join(directory, 'ways_count.csv')
    scanned_files = set()

    # Check if output file exists and load already processed files
    if os.path.exists(output_file):
        try:
            with open(output_file, 'r', newline='') as f:
                reader = csv.reader(f)
                next(reader, None)  # Skip header, safely
                for row in reader:
                    if len(row) >= 3:  # Ensure the row has enough columns
                        scanned_files.add(row[2])  # Add scanned file path to the set
        except Exception as e:
            logger.error("Error reading existing CSV: %s", e)
    else:
        # Create the output file and write the header
        with open(output_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['name', 'ways', 'path'])
            logger.info("Created output file: %s", output_file)

    logger.info("Scanning directory: %s", directory)  # Log current directory being scanned
    for root, dirs, files in os.walk(directory):
        # Skip archive directories
        if 'archive' in root.lower():
            logger.info("Ignoring archive directory: %s", root)
            continue

        # Look for the first osm.pbf file using next() with a generator expression
        osm_file_path = next((os.path.join(root, file) for file in files if file.endswith('.osm.pbf')), None)

        if osm_file_path is not None:
            if osm_file_path in scanned_files:
                logger.info("PBF file already processed: %s", osm_file_path)  # Log already processed directory
                continue
            else:
                # Extract network name from the file name or directory name
                network_name = os.path.basename(root)  # Use the directory name as the network name
                number_of_ways = calculate_ways(osm_file_path)

                # Ensure file ends with newline before appending
                if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
                    with open(output_file, 'rb+') as f:
                        f.seek(-1, os.SEEK_END)  # Go to the last byte
                        last_char = f.read(1)
                        if last_char != b'\n':
                            f.seek(0, os.SEEK_END)  # Go to the end of the file
                            f.write(b'\n')  # Add a newline if it doesn't end with one

                # Append result to the output CSV file
                with open(output_file, 'a', newline='') as f:
                    writer = csv.writer(f)
                    writer.writerow([network_name, number_of_ways, osm_file_path])
                    # Write network name, number of ways, and path
                    logger.info("Appended to CSV: %s, %s, %s", network_name, number_of_ways, osm_file_path)
        else:
            logger.info("No OSM file found in this directory: %s.", root)  # Log message if no file found
            continue  # Skip to the next directory if no file is found
