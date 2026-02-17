"""PBF analysis: OSMTagHandler, stats."""

import logging
import sys
from collections import Counter, defaultdict

import osmium
import pandas as pd

logger = logging.getLogger(__name__)


class OSMTagHandler(osmium.SimpleHandler):
    def __init__(self):
        osmium.SimpleHandler.__init__(self)
        # Create separate counters for different element types
        self.way_tag_counters = defaultdict(Counter)
        self.node_tag_counters = defaultdict(Counter)
        self.relation_tag_counters = defaultdict(Counter)
        self.other_tags_counters = defaultdict(Counter)
        self.records = []
        self.unique_tags = set()
        self.total_ways = 0
        self.total_nodes = 0
        self.total_relations = 0

    def process_other_tags(self, other_tags_str):
        """Parse hstore-formatted other_tags string into a dictionary"""
        if not other_tags_str:
            return {}

        parsed_tags = {}
        try:
            # Handle the format: "key"=>"value","key2"=>"value2",...
            current = ""
            in_quotes = False
            key = None
            parts = []

            # First split into key=>value parts
            for char in other_tags_str:
                if char == '"' and (not current or current[-1] != '\\'):
                    in_quotes = not in_quotes

                current += char

                if char == ',' and not in_quotes:
                    parts.append(current[:-1])  # Remove the trailing comma
                    current = ""

            if current:  # Add the last part if there is one
                parts.append(current)

            # Now process each part to extract key and value
            for part in parts:
                if "=>" in part:
                    key_val = part.split("=>")
                    if len(key_val) == 2:
                        k = key_val[0].strip().strip('"')
                        v = key_val[1].strip().strip('"')
                        parsed_tags[k] = v

                        # Update counter for this key-value pair
                        self.other_tags_counters[k][v] += 1
        except Exception as e:
            logger.error("Error parsing other_tags: %s, value: %s", e, other_tags_str[:100])

        return parsed_tags

    def way(self, w):
        """Process a way and its tags"""
        self.total_ways += 1

        # Extract all tags into a dictionary
        tags_dict = {}
        other_tags_dict = {}

        for tag in w.tags:
            tag_key = tag.k
            tag_value = tag.v

            # Add to our unique tags set
            self.unique_tags.add(tag_key)

            # Store the tag and update its counter for ways
            tags_dict[tag_key] = tag_value
            self.way_tag_counters[tag_key][tag_value] += 1

            # Check if this is an other_tags field that needs parsing
            if tag_key == 'other_tags':
                other_tags_dict = self.process_other_tags(tag_value)

        # Store the record with all its tags
        record = {
            'id': w.id,
            **tags_dict,
            'other_tags_parsed': other_tags_dict
        }

        self.records.append(record)

    # Also handle nodes and relations if needed
    def node(self, n):
        self.total_nodes += 1
        for tag in n.tags:
            self.unique_tags.add(tag.k)
            self.node_tag_counters[tag.k][tag.v] += 1

    def relation(self, r):
        self.total_relations += 1
        for tag in r.tags:
            self.unique_tags.add(tag.k)
            self.relation_tag_counters[tag.k][tag.v] += 1


def analyze_osm_pbf(file_path, num_top_values=10):
    """
    Analyze an OSM PBF file and return statistics about all tags,
    separated by element type (way, node, relation)

    Args:
        file_path: Path to the OSM PBF file
        num_top_values: Optional, Number of top values to report for each tag

    Returns:
        Dictionary of statistics and DataFrame of records
    """
    logger.info("Analyzing OSM PBF file: %s", file_path)
    handler = OSMTagHandler()

    # Process the file
    handler.apply_file(file_path)

    logger.info("Processed %d ways, %d nodes, and %d relations", handler.total_ways, handler.total_nodes, handler.total_relations)
    logger.info("Found %d unique tag keys", len(handler.unique_tags))

    # Create summary statistics for ways
    way_stats = {}
    for tag_key, counter in handler.way_tag_counters.items():
        total = sum(counter.values())
        way_stats[tag_key] = {
            'count': total,
            'unique_values': len(counter),
            'top_values': dict(counter.most_common(num_top_values)),
            'percent_present':
                round(total / handler.total_ways * 100, 2) if handler.total_ways > 0 else 0
        }

    # Sort way stats by frequency
    way_stats = {k: v for k, v in sorted(
        way_stats.items(),
        key=lambda item: item[1]['count'],
        reverse=True
    )}

    # Create summary statistics for nodes
    node_stats = {}
    for tag_key, counter in handler.node_tag_counters.items():
        total = sum(counter.values())
        node_stats[tag_key] = {
            'count': total,
            'unique_values': len(counter),
            'top_values': dict(counter.most_common(num_top_values)),
            'percent_present':
                round(total / handler.total_nodes * 100, 2) if handler.total_nodes > 0 else 0
        }

    # Sort node stats by frequency
    node_stats = {k: v for k, v in sorted(
        node_stats.items(),
        key=lambda item: item[1]['count'],
        reverse=True
    )}

    # Create summary statistics for relations
    relation_stats = {}
    for tag_key, counter in handler.relation_tag_counters.items():
        total = sum(counter.values())
        relation_stats[tag_key] = {
            'count': total,
            'unique_values': len(counter),
            'top_values': dict(counter.most_common(num_top_values)),
            'percent_present':
                round(total / handler.total_relations * 100, 2) if handler.total_relations > 0 else 0
        }

    # Sort relation stats by frequency
    relation_stats = {k: v for k, v in sorted(
        relation_stats.items(),
        key=lambda item: item[1]['count'],
        reverse=True
    )}

    # Create similar statistics for other_tags fields
    other_tags_stats = {}
    for tag_key, counter in handler.other_tags_counters.items():
        total = sum(counter.values())
        other_tags_stats[tag_key] = {
            'count': total,
            'unique_values': len(counter),
            'top_values': dict(counter.most_common(num_top_values)),
            'percent_present':
                round(total / (handler.total_ways + handler.total_nodes + handler.total_relations) * 100, 2)
        }

    # Sort other_tags stats by frequency
    other_tags_stats = {k: v for k, v in sorted(
        other_tags_stats.items(),
        key=lambda item: item[1]['count'],
        reverse=True
    )}

    # Create a DataFrame from the records
    records_df = pd.DataFrame(handler.records) if handler.records else pd.DataFrame()

    # Return the summary statistics and records
    return {
        'total_ways': handler.total_ways,
        'total_nodes': handler.total_nodes,
        'total_relations': handler.total_relations,
        'unique_tags': list(handler.unique_tags),
        'way_stats': way_stats,
        'node_stats': node_stats,
        'relation_stats': relation_stats,
        'other_tags_stats': other_tags_stats
    }, records_df


def print_tag_stats(stats, category_name="Tags", element_type="Elements", limit=None):
    """Print tag statistics in a formatted way"""
    logger.info("\n=== %s Statistics for %s ===", category_name, element_type)
    logger.info("Total unique %s: %d", category_name.lower(), len(stats))

    for i, (tag, data) in enumerate(stats.items()):
        if limit and i >= limit:
            logger.info("\n... and %d more %s.", len(stats) - limit, category_name.lower())
            break

        logger.info("\n%d. %s: %d instances (%s%% of %s)", i + 1, tag, data['count'], data['percent_present'], element_type.lower())
        logger.info("   Unique values: %d", data['unique_values'])
        logger.info("   Top values:")

        # Print top values with their counts
        for val, count in data['top_values'].items():
            # Truncate very long values
            display_val = val[:50] + "..." if len(val) > 50 else val
            logger.info("     - %s: %d", display_val, count)


def main(file_path=None):
    """Main function to analyze an OSM PBF file"""
    if not file_path:
        logger.info("No file provided. To analyze a file, run: python osm_analyzer.py <file.osm.pbf>")
        return

    # Analyze the PBF file
    stats, records_df = analyze_osm_pbf(file_path, 30)

    logger.info("\n=== OSM PBF Analysis Summary ===")
    logger.info("Total ways processed: %d", stats['total_ways'])
    logger.info("Total nodes processed: %d", stats['total_nodes'])
    logger.info("Total relations processed: %d", stats['total_relations'])
    logger.info("Total unique tags found: %d", len(stats['unique_tags']))

    # Print way tag statistics
    print_tag_stats(stats['way_stats'], "Way Tags", "Ways", limit=20)

    # Print node tag statistics
    print_tag_stats(stats['node_stats'], "Node Tags", "Nodes", limit=20)

    # Print relation tag statistics
    print_tag_stats(stats['relation_stats'], "Relation Tags", "Relations", limit=20)

    # Print other_tags statistics
    print_tag_stats(stats['other_tags_stats'], "other_tags Keys", "All Elements", limit=20)

    # Show column names in the data
    if not records_df.empty:
        logger.info("\n=== DataFrame Columns ===")
        columns = list(records_df.columns)
        for i, col in enumerate(columns):
            logger.info("%d. %s", i + 1, col)

    return stats, records_df


if __name__ == "__main__":
    if len(sys.argv) < 2:
        main()
    else:
        main(sys.argv[1])
