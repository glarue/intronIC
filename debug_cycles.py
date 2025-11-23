#!/usr/bin/env python3
"""
Debug script to analyze cycle generation in annotation parsing.

Usage:
    python debug_cycles.py <annotation_file.gff.gz>
"""

import sys
from pathlib import Path
from collections import defaultdict
from networkx import DiGraph, find_cycle, is_directed_acyclic_graph

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent))

from extraction.annotator import AnnotationHierarchyBuilder
from file_io.parsers import BioGLAnnotationParser


def analyze_cycles(annotation_file: str, clean_names: bool = True):
    """Analyze cycles in annotation graph without removing them."""

    print(f"Analyzing: {annotation_file}")
    print("=" * 80)

    # Parse annotation file
    parser = BioGLAnnotationParser(clean_names=clean_names)
    annotations = list(parser.parse_file(annotation_file))
    print(f"Parsed {len(annotations):,} annotation lines")

    # Build graph (copied from AnnotationHierarchyBuilder)
    builder = AnnotationHierarchyBuilder(
        child_features=['cds', 'exon'],
        clean_names=clean_names
    )

    # Build graph WITHOUT cycle removal
    feat_graph = DiGraph()
    feat_index = {}
    unique_coords = defaultdict(set)

    for ann in annotations:
        if ann.start > ann.stop:
            continue

        features = builder._create_features_from_annotation(ann)

        for feat in features:
            feat_type = feat.attributes.get('_orig_feat_type', '')
            if not feat_type and hasattr(feat, 'feature_type'):
                feat_type = feat.feature_type
            feat_type = feat_type.lower() if feat_type else 'unknown'

            parent_name = getattr(feat, 'parent_id', None) or feat.attributes.get('_parent_name')
            grandparent_name = feat.attributes.get('_grandparent_name')

            # Create name
            if feat_type in builder.child_features:
                name = f"{feat_type}_Parent={parent_name}:{feat.start}_{feat.stop}"
            else:
                name = feat.feature_id

            # Check duplicates
            if feat_type in builder.child_features and parent_name:
                check_coords = (feat_type, feat.start, feat.stop)
                if check_coords in unique_coords[parent_name]:
                    continue
                unique_coords[parent_name].add(check_coords)

            # Add to index
            feat_index[name] = feat

            # Build graph edges
            if parent_name is not None:
                feat_graph.add_edge(parent_name, name)

            if grandparent_name is not None:
                feat_graph.add_edge(grandparent_name, parent_name)

    print(f"Built graph: {feat_graph.number_of_nodes():,} nodes, {feat_graph.number_of_edges():,} edges")

    # Check for cycles
    print("\nChecking for cycles...")

    if is_directed_acyclic_graph(feat_graph):
        print("✓ Graph is acyclic (no cycles)")
        return

    print("✗ Graph contains cycles!")
    print("\nFinding cycles...")

    all_cycles = []
    max_cycles_to_show = 20

    # Find up to max_cycles_to_show cycles
    for i in range(1000):
        try:
            cycle_edges = list(find_cycle(feat_graph))
            if cycle_edges:
                all_cycles.append(cycle_edges)
                # Remove edges to find next cycle
                feat_graph.remove_edges_from(cycle_edges)
            else:
                break
        except:
            # No more cycles
            break

    print(f"\nFound {len(all_cycles)} cycle(s)")
    print("\n" + "=" * 80)

    # Show first few cycles in detail
    for i, cycle_edges in enumerate(all_cycles[:max_cycles_to_show]):
        print(f"\nCycle {i+1}:")
        print("-" * 80)

        # Show the cycle path
        nodes_in_cycle = set()
        for parent, child in cycle_edges:
            nodes_in_cycle.add(parent)
            nodes_in_cycle.add(child)
            print(f"  {parent}")
            print(f"    ↓")
            print(f"  {child}")

        print(f"\nNodes in cycle: {len(nodes_in_cycle)}")

        # Try to identify what types of features are involved
        for node in sorted(nodes_in_cycle):
            if node in feat_index:
                feat = feat_index[node]
                feat_type = feat.attributes.get('_orig_feat_type', 'unknown')
                print(f"  - {node[:80]} ({feat_type})")
            else:
                print(f"  - {node[:80]} (not in feat_index - placeholder?)")

    if len(all_cycles) > max_cycles_to_show:
        print(f"\n... and {len(all_cycles) - max_cycles_to_show} more cycles")

    # Summary statistics
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)

    all_cycle_edges = [edge for cycle in all_cycles for edge in cycle]
    unique_nodes = set()
    for parent, child in all_cycle_edges:
        unique_nodes.add(parent)
        unique_nodes.add(child)

    print(f"Total cycles found: {len(all_cycles)}")
    print(f"Total edges in cycles: {len(all_cycle_edges)}")
    print(f"Unique nodes involved: {len(unique_nodes)}")

    # Check for self-loops
    self_loops = [(p, c) for cycle in all_cycles for p, c in cycle if p == c]
    if self_loops:
        print(f"\nSelf-loops found: {len(self_loops)}")
        for parent, child in self_loops[:10]:
            print(f"  - {parent} -> {child}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python debug_cycles.py <annotation_file.gff.gz>")
        sys.exit(1)

    annotation_file = sys.argv[1]
    analyze_cycles(annotation_file)
