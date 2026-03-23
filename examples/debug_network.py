"""Diagnostic script: trace disconnected fragments and bridge corridor integrity.

Usage
-----
    python examples/debug_network.py \\
        --pkl ~/Workspace/Simulation/sfbay/network/sfbay-cbg5500-weakConn-network/sfbay-cbg5500-weakConn-network.pkl \\
        --raw-pkl ~/Workspace/Simulation/sfbay/network/raw_osm_graph_sfbay_060704d91c.pkl

What it does
------------
1.  Loads the *processed* .pkl (final graph after the full pipeline).
2.  Finds every weakly-connected component that is NOT the largest one
    (i.e. the fragments that were discarded by the LCC step).
    For each fragment it reports the node count, edge count, highway type
    breakdown, bridge/tunnel edge count, and geographic bounding box.
3.  Checks whether well-known Bay Area bridge corridors are reachable in
    the processed graph (directed AND undirected).
4.  If --raw-pkl is supplied, runs the same corridor check on the raw graph
    so you can compare which bridges exist in the raw download but not in
    the processed network.
"""

import argparse
import pickle
import sys
from collections import Counter
from pathlib import Path

EXAMPLES_DIR = Path(__file__).resolve().parent
if str(EXAMPLES_DIR) not in sys.path:
    sys.path.insert(0, str(EXAMPLES_DIR))

from _bootstrap import bootstrap_example_paths

bootstrap_example_paths(__file__)

import networkx as nx
import osmnx as ox

# ---------------------------------------------------------------------------
# Well-known Bay Area bridge corridors
# Each entry: (name, west_lon, south_lat, east_lon, north_lat)
# ---------------------------------------------------------------------------
# Each corridor entry:
#   (name, entry_lon, entry_lat, exit_lon, exit_lat, bbox_west, bbox_south, bbox_east, bbox_north)
# entry = where traffic enters from one side, exit = where it should emerge on the other side.
# The directed check finds the nearest node to (entry_lon, entry_lat) and checks if any node
# near (exit_lon, exit_lat) is reachable in the directed graph.
BAY_AREA_CORRIDORS = [
    # Bay Bridge west span: westbound traffic enters from YBI (east) and exits at SF (west)
    ("Bay Bridge (west span, westbound)", -122.369, 37.816, -122.395, 37.793,
     -122.405, 37.785, -122.355, 37.830),
    # Bay Bridge east span: westbound traffic enters from Oakland (east) and reaches YBI (west)
    ("Bay Bridge (east span, westbound)", -122.290, 37.814, -122.358, 37.816,
     -122.375, 37.800, -122.265, 37.835),
    # Golden Gate: southbound enters Marin (north), exits SF (south)
    ("Golden Gate Bridge (southbound)", -122.479, 37.829, -122.478, 37.810,
     -122.490, 37.808, -122.465, 37.835),
    # Richmond-San Rafael: westbound enters Richmond (east), exits San Rafael (west)
    # exit coords set to I-580/US-101 interchange area in San Rafael
    ("Richmond-San Rafael Bridge (westbound)", -122.360, 37.941, -122.530, 37.970,
     -122.460, 37.920, -122.340, 37.960),
    # San Mateo: westbound enters Hwy 92 near I-880 (east), exits Hwy 92 near US-101 (west)
    ("San Mateo Bridge (westbound)", -122.065, 37.608, -122.275, 37.556,
     -122.245, 37.560, -122.055, 37.620),
    # Dumbarton: westbound enters east side (Newark), exits west (East Palo Alto)
    # entry updated — original was ~3 km west of actual east approach; bbox extended east
    ("Dumbarton Bridge (westbound)", -121.990, 37.508, -122.115, 37.508,
     -122.130, 37.495, -121.975, 37.520),
    # Caldecott Tunnel: westbound enters Orinda side, exits Oakland/Berkeley
    ("Caldecott Tunnel (westbound)", -122.195, 37.845, -122.225, 37.845,
     -122.230, 37.835, -122.185, 37.860),
    # Broadway Tunnel: westbound enters Bush St side (east), exits Polk St (west)
    ("Broadway Tunnel (westbound)", -122.413, 37.800, -122.423, 37.800,
     -122.428, 37.795, -122.408, 37.806),
]


def _normalize_hw(value) -> str:
    if isinstance(value, list):
        return str(value[0]).strip() if value else "unknown"
    if value is None:
        return "unknown"
    return str(value).strip()


def _is_bridge_or_tunnel(data: dict) -> bool:
    for key in ("bridge", "tunnel"):
        v = data.get(key)
        if v and str(v).strip().lower() in {"yes", "true", "1"}:
            return True
    return False


def _edges_in_bbox(G, west, south, east, north):
    """Return list of (u, v, k, data) for edges whose geometry touches the bbox."""
    result = []
    for u, v, k, data in G.edges(data=True, keys=True):
        geom = data.get("geometry")
        if geom is not None:
            x1, y1, x2, y2 = geom.bounds
            if x1 <= east and x2 >= west and y1 <= north and y2 >= south:
                result.append((u, v, k, data))
        else:
            # fall back to node coordinates
            nu = G.nodes.get(u, {})
            nv = G.nodes.get(v, {})
            xu, yu = nu.get("x"), nu.get("y")
            xv, yv = nv.get("x"), nv.get("y")
            if all(c is not None for c in (xu, yu, xv, yv)):
                if (west <= xu <= east or west <= xv <= east) and (
                    south <= yu <= north or south <= yv <= north
                ):
                    result.append((u, v, k, data))
    return result


def _nodes_in_bbox(G, west, south, east, north):
    return [
        n
        for n, d in G.nodes(data=True)
        if west <= d.get("x", 999) <= east and south <= d.get("y", 999) <= north
    ]


def analyze_fragments(G):
    """Return info on all non-largest weakly-connected components."""
    components = sorted(nx.weakly_connected_components(G), key=len, reverse=True)
    largest = components[0]
    fragments = components[1:]
    print(f"\nMain component: {len(largest):,} nodes")
    print(f"Disconnected fragments: {len(fragments)} (total {sum(len(c) for c in fragments):,} nodes)\n")

    for i, comp in enumerate(fragments):
        sub = G.subgraph(comp)
        hw_counts: Counter = Counter()
        bt_count = 0
        lons, lats = [], []
        for n, d in sub.nodes(data=True):
            x, y = d.get("x"), d.get("y")
            if x is not None and y is not None:
                lons.append(x)
                lats.append(y)
        for u, v, k, data in sub.edges(data=True, keys=True):
            hw_counts[_normalize_hw(data.get("highway"))] += 1
            if _is_bridge_or_tunnel(data):
                bt_count += 1

        major = {k: v for k, v in hw_counts.items()
                 if k in {"motorway", "motorway_link", "trunk", "trunk_link",
                          "primary", "primary_link"}}
        bbox_str = ""
        if lons:
            bbox_str = (f"lon [{min(lons):.4f}, {max(lons):.4f}]"
                        f" lat [{min(lats):.4f}, {max(lats):.4f}]")

        if major or bt_count > 0:
            print(
                f"  Fragment {i+1}: {len(comp)} nodes, {sub.number_of_edges()} edges"
                f" | bridge/tunnel: {bt_count}"
                f" | major hw: {dict(major)}"
                f" | {bbox_str}"
            )


def _nearby_nodes(G, lon, lat, radius_deg=0.01, top_n=10, highway_types=None):
    """Return up to top_n node IDs closest to (lon, lat) within radius_deg.

    If *highway_types* is given (a set of OSM highway values), only nodes that
    have at least one outgoing edge matching those types are returned.  Falls
    back to all nodes if no highway-filtered candidates are found within radius.
    """
    candidates = []
    for n, d in G.nodes(data=True):
        x, y = d.get("x"), d.get("y")
        if x is None or y is None:
            continue
        dist2 = (x - lon) ** 2 + (y - lat) ** 2
        if dist2 <= radius_deg ** 2:
            candidates.append((dist2, n))
    candidates.sort()

    if highway_types is None:
        return [n for _, n in candidates[:top_n]]

    # Prefer nodes with a major-highway outgoing edge
    filtered = [
        (d2, n) for d2, n in candidates
        if any(
            _normalize_hw(edata.get("highway")) in highway_types
            for _, _, edata in G.out_edges(n, data=True)
        )
    ]
    if filtered:
        return [n for _, n in filtered[:top_n]]
    # Fall back to unfiltered if none found
    return [n for _, n in candidates[:top_n]]


def _any_path(G, sources, targets):
    """Return True if any node in sources can reach any node in targets."""
    if not sources or not targets:
        return False
    # Use multi-source BFS from all sources simultaneously for efficiency
    try:
        reachable = nx.descendants(G, sources[0])
        reachable.add(sources[0])
        for s in sources[1:]:
            reachable |= nx.descendants(G, s)
            reachable.add(s)
        return bool(reachable & set(targets))
    except (nx.NodeNotFound, nx.NetworkXError):
        return False


def check_corridors(G, label):
    """Check bridge/tunnel edge presence and directed routing per corridor.

    For each corridor we sample up to 10 nodes near the 'entry' coordinates
    and up to 10 nodes near the 'exit' coordinates, then test whether ANY
    entry node can reach ANY exit node via the directed graph.  This is more
    robust than nearest-single-node because it avoids landing on the wrong
    one-way carriageway.
    """
    print(f"\n--- Corridor check: {label} ---")

    for name, entry_lon, entry_lat, exit_lon, exit_lat, west, south, east, north in BAY_AREA_CORRIDORS:
        edges = _edges_in_bbox(G, west, south, east, north)
        nodes_in = _nodes_in_bbox(G, west, south, east, north)
        bt_edges = [e for e in edges if _is_bridge_or_tunnel(e[3])]

        entry_nodes = _nearby_nodes(G, entry_lon, entry_lat, highway_types=_MAJOR_HW)
        exit_nodes = _nearby_nodes(G, exit_lon, exit_lat, highway_types=_MAJOR_HW)
        has_directed_path = _any_path(G, entry_nodes, exit_nodes)

        status = "OK" if bt_edges else ("MISSING" if not edges else "no bridge/tunnel tag")
        routed = "ROUTABLE" if has_directed_path else "NOT ROUTABLE"
        print(
            f"  {name}:"
            f" {len(edges)} edges ({len(bt_edges)} bridge/tunnel)"
            f" | {len(nodes_in)} nodes in bbox"
            f" | entry→exit directed: {routed}"
            f" | entry nodes: {entry_nodes[:3]} exit nodes: {exit_nodes[:3]}"
            f" | {status}"
        )


def dump_corridor_edges(G, name, west, south, east, north):
    """Print full edge details for a corridor bbox to diagnose direction issues."""
    print(f"\n=== Edge dump: {name} ===")
    edges = _edges_in_bbox(G, west, south, east, north)
    print(f"  {len(edges)} edges in bbox")
    for u, v, k, data in sorted(edges, key=lambda e: (e[0], e[1])):
        u_data = G.nodes.get(u, {})
        v_data = G.nodes.get(v, {})
        print(
            f"  {u} ({u_data.get('x', '?'):.4f},{u_data.get('y', '?'):.4f})"
            f" → {v} ({v_data.get('x', '?'):.4f},{v_data.get('y', '?'):.4f})"
            f"  hw={data.get('highway','?')}"
            f"  bridge={data.get('bridge','?')}"
            f"  oneway={data.get('oneway','?')}"
            f"  reversed={data.get('reversed','?')}"
            f"  length={data.get('length','?')}"
        )


_MAJOR_HW = {
    "motorway", "motorway_link", "trunk", "trunk_link", "primary", "primary_link",
}


def _major_hw_nodes_in_bbox(G, west, south, east, north):
    """Return node IDs that have at least one outgoing edge on a major highway within bbox."""
    result = []
    for n, d in G.nodes(data=True):
        x, y = d.get("x"), d.get("y")
        if x is None or y is None:
            continue
        if not (west <= x <= east and south <= y <= north):
            continue
        for _, _, edata in G.out_edges(n, data=True):
            hw = _normalize_hw(edata.get("highway"))
            if hw in _MAJOR_HW:
                result.append(n)
                break
    return result


def _multi_source_bfs(G, sources):
    """Return the set of all nodes reachable (directed) from any node in *sources*."""
    from collections import deque
    visited = set(sources)
    queue = deque(sources)
    while queue:
        n = queue.popleft()
        for _, v in G.out_edges(n):
            if v not in visited:
                visited.add(v)
                queue.append(v)
    return visited


def trace_approach(G, name):
    """Trace directed routing from major-highway approach nodes into the bridge bbox.

    Searches an extended approach zone (0.02° padding beyond the entry point) for
    motorway/trunk/primary nodes, then checks which of those can reach nodes inside
    the bridge/tunnel bbox.  This avoids the false-negative where _nearby_nodes lands
    on an unclassified local road rather than the actual highway approach.
    """
    corridor = next((c for c in BAY_AREA_CORRIDORS if c[0] == name), None)
    if corridor is None:
        print(f"Unknown corridor: {name}")
        return

    cname, entry_lon, entry_lat, exit_lon, exit_lat, west, south, east, north = corridor

    # Extend the entry search zone by 0.02° beyond the entry point
    pad = 0.02
    ew = min(entry_lon, west) - pad
    ee = max(entry_lon, east) + pad
    es = min(entry_lat, south) - pad
    en = max(entry_lat, north) + pad

    approach_nodes = _major_hw_nodes_in_bbox(G, ew, es, ee, en)
    bridge_nodes = _nodes_in_bbox(G, west, south, east, north)
    bridge_node_set = set(bridge_nodes)

    print(f"\n=== Approach trace: {name} ===")
    print(f"  Approach search zone: lon [{ew:.4f}, {ee:.4f}] lat [{es:.4f}, {en:.4f}]")
    print(f"  Major-highway approach nodes found: {len(approach_nodes)}")
    print(f"  Bridge bbox nodes: {len(bridge_nodes)}")

    if not approach_nodes:
        print("  !! No major-highway nodes found in approach zone — approach is missing")
        return
    if not bridge_nodes:
        print("  !! No nodes inside bridge bbox — bridge may be missing entirely")
        return

    # Forward BFS from all approach nodes — which bridge bbox nodes are reachable?
    reachable_from_approach = _multi_source_bfs(G, approach_nodes)
    reachable_bridge_nodes = reachable_from_approach & bridge_node_set

    print(
        f"  Bridge bbox nodes reachable from approach: "
        f"{len(reachable_bridge_nodes)} / {len(bridge_nodes)}"
    )

    if reachable_bridge_nodes:
        # Show 3 approach nodes closest to entry that can reach the bridge
        reachable_from = [
            n for n in approach_nodes
            if any(v in reachable_bridge_nodes for _, v in G.out_edges(n))
            or n in bridge_node_set
        ]
        closest = sorted(
            reachable_from or approach_nodes,
            key=lambda n: (G.nodes[n].get("x", 999) - entry_lon) ** 2
                         + (G.nodes[n].get("y", 999) - entry_lat) ** 2,
        )[:3]
        print("  Sample approach nodes near entry that feed into bridge:")
        for n in closest:
            nd = G.nodes[n]
            out_hw = [_normalize_hw(d.get("highway")) for _, _, d in G.out_edges(n, data=True)]
            print(
                f"    node {n} ({nd.get('x', '?'):.4f},{nd.get('y', '?'):.4f})"
                f"  out_highways={out_hw[:5]}"
            )
    else:
        # No approach node can reach bridge — show closest for manual inspection
        closest = sorted(
            approach_nodes,
            key=lambda n: (G.nodes[n].get("x", 999) - entry_lon) ** 2
                         + (G.nodes[n].get("y", 999) - entry_lat) ** 2,
        )[:5]
        print("  !! No approach node can reach bridge bbox. Closest approach nodes:")
        for n in closest:
            nd = G.nodes[n]
            out_edges = list(G.out_edges(n, data=True))
            print(
                f"    node {n} ({nd.get('x', '?'):.4f},{nd.get('y', '?'):.4f})"
                f"  out_edges={len(out_edges)}"
            )
            for _, v, d in out_edges[:3]:
                vd = G.nodes.get(v, {})
                print(
                    f"      → {v} ({vd.get('x', '?'):.4f},{vd.get('y', '?'):.4f})"
                    f"  hw={_normalize_hw(d.get('highway'))}"
                    f"  bridge={d.get('bridge','?')}"
                    f"  oneway={d.get('oneway','?')}"
                )

    # Exit side: reverse BFS from exit approach nodes — which bridge bbox nodes can reach exit?
    exit_pad = 0.02
    xw = min(exit_lon, west) - exit_pad
    xe = max(exit_lon, east) + exit_pad
    xs = min(exit_lat, south) - exit_pad
    xn = max(exit_lat, north) + exit_pad
    exit_approach_nodes = _major_hw_nodes_in_bbox(G, xw, xs, xe, xn)

    # Reverse BFS: predecessors of exit nodes — finds what can reach them
    G_rev = G.reverse(copy=False)
    reachable_to_exit = _multi_source_bfs(G_rev, exit_approach_nodes)
    bridge_nodes_reaching_exit = bridge_node_set & reachable_to_exit

    print(
        f"  Bridge bbox nodes that can reach exit approach: "
        f"{len(bridge_nodes_reaching_exit)} / {len(bridge_nodes)}"
    )


def geometry_gaps(G, west, south, east, north, threshold_m=5.0):
    """Find edges in bbox where geometry endpoints don't match node coordinates.

    After ``ox.consolidate_intersections`` the node (x, y) is updated to the
    centroid of the merged cluster, but the edge ``geometry`` LineString may
    still store the old pre-merge endpoint coordinates.  This creates a visual
    gap in QGIS even though the graph is topologically intact.
    """
    DEG_TO_M = 111_000.0

    def _dist_m(ax, ay, bx, by):
        dx = (ax - bx) * DEG_TO_M
        dy = (ay - by) * DEG_TO_M
        return (dx * dx + dy * dy) ** 0.5

    print(f"\n=== Geometry-gap check (threshold {threshold_m} m) ===")
    print(f"  bbox: lon [{west:.4f}, {east:.4f}] lat [{south:.4f}, {north:.4f}]")

    edges = _edges_in_bbox(G, west, south, east, north)
    gaps = []
    for u, v, k, data in edges:
        geom = data.get("geometry")
        if geom is None:
            continue
        coords = list(geom.coords)
        if len(coords) < 2:
            continue
        gx0, gy0 = coords[0]
        gx1, gy1 = coords[-1]
        nu = G.nodes.get(u, {})
        nv = G.nodes.get(v, {})
        nx_u, ny_u = nu.get("x"), nu.get("y")
        nx_v, ny_v = nv.get("x"), nv.get("y")
        if None in (nx_u, ny_u, nx_v, ny_v):
            continue
        gap_start = _dist_m(gx0, gy0, nx_u, ny_u)
        gap_end   = _dist_m(gx1, gy1, nx_v, ny_v)
        if gap_start > threshold_m or gap_end > threshold_m:
            gaps.append((max(gap_start, gap_end), u, v, k, data,
                         gap_start, gap_end,
                         gx0, gy0, gx1, gy1,
                         nx_u, ny_u, nx_v, ny_v))

    if not gaps:
        print(f"  No geometry gaps > {threshold_m} m found.")
        return

    gaps.sort(reverse=True)
    print(f"  {len(gaps)} edges with geometry gap > {threshold_m} m:\n")
    for max_gap, u, v, k, data, gs, ge, gx0, gy0, gx1, gy1, nxu, nyu, nxv, nyv in gaps[:30]:
        print(f"  {u}→{v} k={k}  hw={data.get('highway','?')}  bridge={data.get('bridge','?')}  length={data.get('length','?')}")
        if gs > threshold_m:
            print(f"    start gap {gs:.1f} m: geom=({gx0:.6f},{gy0:.6f})  node_u=({nxu:.6f},{nyu:.6f})")
        if ge > threshold_m:
            print(f"    end   gap {ge:.1f} m: geom=({gx1:.6f},{gy1:.6f})  node_v=({nxv:.6f},{nyv:.6f})")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pkl", required=True, help="Path to processed .pkl graph")
    parser.add_argument("--raw-pkl", help="Path to raw cached .pkl graph (optional)")
    _corridor_names = [c[0] for c in BAY_AREA_CORRIDORS]
    parser.add_argument(
        "--dump-corridor",
        choices=_corridor_names,
        metavar="NAME",
        help=(
            "Print full edge list for a specific corridor. "
            f"Choices: {_corridor_names}"
        ),
    )
    parser.add_argument(
        "--trace-approach",
        choices=_corridor_names,
        metavar="NAME",
        help=(
            "Trace directed routing from major-highway approach nodes into a corridor bbox. "
            "More reliable than --dump-corridor for diagnosing severed approach connections. "
            f"Choices: {_corridor_names}"
        ),
    )
    parser.add_argument(
        "--geometry-gaps",
        choices=_corridor_names,
        metavar="NAME",
        help="Scan edges in a corridor bbox for geometry-node coordinate mismatches.",
    )
    parser.add_argument(
        "--gap-threshold",
        type=float,
        default=5.0,
        help="Metres threshold for --geometry-gaps (default: 5.0)",
    )
    args = parser.parse_args()

    print(f"Loading processed graph: {args.pkl}")
    with open(Path(args.pkl).expanduser(), "rb") as f:
        G = pickle.load(f)
    print(f"  {G.number_of_nodes():,} nodes, {G.number_of_edges():,} edges")

    analyze_fragments(G)
    check_corridors(G, "processed graph")

    if args.raw_pkl:
        print(f"\nLoading raw graph: {args.raw_pkl}")
        with open(Path(args.raw_pkl).expanduser(), "rb") as f:
            G_raw = pickle.load(f)
        print(f"  {G_raw.number_of_nodes():,} nodes, {G_raw.number_of_edges():,} edges")
        analyze_fragments(G_raw)
        check_corridors(G_raw, "raw graph")

    if args.dump_corridor:
        corridor = next(c for c in BAY_AREA_CORRIDORS if c[0] == args.dump_corridor)
        _, _, _, _, _, west, south, east, north = corridor
        dump_corridor_edges(G, args.dump_corridor, west, south, east, north)
        if args.raw_pkl:
            print("\n--- same corridor in raw graph ---")
            dump_corridor_edges(G_raw, args.dump_corridor, west, south, east, north)

    if args.trace_approach:
        trace_approach(G, args.trace_approach)
        if args.raw_pkl:
            print("\n--- same approach trace in raw graph ---")
            trace_approach(G_raw, args.trace_approach)

    if args.geometry_gaps:
        corridor = next(c for c in BAY_AREA_CORRIDORS if c[0] == args.geometry_gaps)
        _, _, _, _, _, west, south, east, north = corridor
        geometry_gaps(G, west, south, east, north, threshold_m=args.gap_threshold)
        if args.raw_pkl:
            print("\n--- same geometry-gap check in raw graph ---")
            geometry_gaps(G_raw, west, south, east, north, threshold_m=args.gap_threshold)



if __name__ == "__main__":
    main()