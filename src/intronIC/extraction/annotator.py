"""
Build hierarchical gene/transcript/exon structures from annotation files.

This module handles parsing annotation files (GFF3/GTF) and building the
parent-child relationships between genes, transcripts, and exons/CDS features.
"""

from collections import defaultdict
from typing import Dict, Iterable, List, Optional, Set, Tuple

from networkx import (
    DiGraph,
    find_cycle,
    is_directed_acyclic_graph,
    lexicographical_topological_sort,
)

from intronIC.core.models import Exon, Gene, Transcript
from intronIC.file_io.parsers import AnnotationLine, BioGLAnnotationParser
from intronIC.utils.coordinates import CoordinateSystem, GenomicCoordinate


class AnnotationHierarchyBuilder:
    """
    Builds hierarchical gene -> transcript -> exon structures from annotations.

    This class takes annotation data and constructs a directed acyclic graph (DAG)
    of parent-child relationships, then creates properly linked Gene, Transcript,
    and Exon objects.

    Examples:
        >>> builder = AnnotationHierarchyBuilder(['exon', 'cds'])
        >>> genes = builder.build_from_file('annotation.gff3')
        >>> print(f"Found {len(genes)} genes")

    Attributes:
        child_features: Feature types to treat as children (e.g., 'exon', 'cds')
        parent_class_map: Maps feature types to their parent/grandparent classes
    """

    def __init__(
        self, child_features: List[str], clean_names: bool = True, messenger=None
    ):
        """
        Initialize the hierarchy builder.

        Args:
            child_features: List of feature types to extract (e.g., ['exon', 'cds'])
            clean_names: If True, remove Ensembl-style 'gene:' and 'transcript:' prefixes.
                         This only removes colon-separated prefixes (Ensembl convention).
                         RefSeq-style hyphen/underscore prefixes are preserved to avoid collisions.
            messenger: Optional UnifiedMessenger for logging (if None, uses print)
        """
        self.child_features = [f.lower() for f in child_features]
        self.clean_names = clean_names
        self.messenger = messenger

        # Maps feature types to expected parent classes
        self.parent_class_map = {
            "parent": {
                "exon": Transcript,
                "cds": Transcript,
                "transcript": Gene,
                "gene": None,
            },
            "grandparent": {
                "exon": Gene,
                "cds": Gene,
                "transcript": None,
                "gene": None,
            },
        }

        # Feature index - populated after building hierarchy
        self.feature_index: Dict[str, Gene | Transcript | Exon] = {}

    @staticmethod
    def _create_coordinate(ann: AnnotationLine) -> GenomicCoordinate:
        """
        Convert annotation line to genomic coordinate.

        Args:
            ann: Parsed annotation line

        Returns:
            GenomicCoordinate object in 1-based system
        """
        return GenomicCoordinate(
            chromosome=ann.region,
            start=ann.start,
            stop=ann.stop,
            strand=ann.strand,
            system="1-based",  # GFF3/GTF are 1-based
        )

    def build_from_file(self, annotation_file: str) -> List[Gene]:
        """
        Build gene hierarchy from an annotation file.

        Args:
            annotation_file: Path to GFF3/GTF annotation file

        Returns:
            List of Gene objects with children populated

        Examples:
            >>> builder = AnnotationHierarchyBuilder(['exon'])
            >>> genes = builder.build_from_file('chr19.gff3')
            >>> gene = genes[0]
            >>> print(f"{gene.name}: {len(gene.children)} transcripts")
        """
        # Parse annotation file (streaming - don't materialize all at once!)
        parser = BioGLAnnotationParser(clean_names=self.clean_names)
        annotations_iter = parser.parse_file(annotation_file)

        # Build hierarchy from iterator (processes as it goes)
        return self.build_from_annotations(annotations_iter)

    def build_from_chromosome(
        self, annotation_file: str, chromosome: str
    ) -> List[Gene]:
        """
        Build gene hierarchy for a single contig only.

        This enables contig-by-contig processing for memory optimization.
        Only genes on the specified contig are returned.

        Args:
            annotation_file: Path to GFF3/GTF annotation file
            chromosome: Contig/chromosome name to filter for

        Returns:
            List of Gene objects for the specified contig

        Examples:
            >>> builder = AnnotationHierarchyBuilder(['exon'])
            >>> chr1_genes = builder.build_from_chromosome('genome.gff3', 'chr1')
            >>> print(f"Chr1 has {len(chr1_genes)} genes")

        Note:
            This is more memory-efficient than build_from_file() for large genomes
            because it only processes features from one contig at a time.
        """
        from intronIC.file_io.parsers import filter_annotation_lines_by_chromosome

        # Get contig-filtered annotation lines (streaming)
        chr_annotations = filter_annotation_lines_by_chromosome(
            annotation_file, chromosome
        )

        # Build hierarchy from filtered annotations
        return self.build_from_annotations(chr_annotations)

    def build_from_lines(self, lines: List[str]) -> List[Gene]:
        """
        Build gene hierarchy from annotation lines (streaming-friendly).

        This method enables memory-efficient streaming processing by accepting
        pre-extracted lines instead of requiring full file parsing.

        Args:
            lines: List of annotation lines (GFF3/GTF format)

        Returns:
            List of Gene objects with children populated

        Examples:
            >>> # Get lines for chr1 from index
            >>> from file_io.annotation_index import build_contig_index, extract_contig_lines
            >>> index = build_contig_index('annotations.gff3.gz')
            >>> chr1_lines = extract_contig_lines('annotations.gff3.gz', index.get_line_numbers('chr1'))
            >>>
            >>> # Build hierarchy from just chr1's lines
            >>> builder = AnnotationHierarchyBuilder(['exon'])
            >>> chr1_genes = builder.build_from_lines(chr1_lines)

        Note:
            This is the key method for streaming annotation processing.
            It allows building gene hierarchies for one contig at a time
            without loading the entire annotation file into memory.
        """
        # Parse lines using same parser
        parser = BioGLAnnotationParser(clean_names=self.clean_names)
        annotations = list(parser.parse_lines(lines))

        # Build hierarchy
        return self.build_from_annotations(annotations)

    def build_from_annotations(
        self, annotations: Iterable[AnnotationLine]
    ) -> List[Gene]:
        """
        Build gene hierarchy from parsed annotation lines.

        Args:
            annotations: Iterable of parsed AnnotationLine objects (can be generator)

        Returns:
            List of top-level Gene objects

        Note:
            Accepts an iterable (including generators) to enable streaming processing.
            This avoids materializing all AnnotationLine objects in memory at once.
        """
        # Step 1: Create features and build graph
        feat_graph = DiGraph()
        feat_index: Dict[str, Gene | Transcript | Exon] = {}
        unique_coords: Dict[str, Set[Tuple]] = defaultdict(set)

        for ann in annotations:
            # Skip invalid features where start > stop (annotation errors)
            # Note: Allow start == stop for 1bp features (now supported by GenomicCoordinate)
            if ann.start > ann.stop:
                continue

            features = self._create_features_from_annotation(ann)

            for feat in features:
                # Get feature type - use stored value from attributes
                feat_type = feat.attributes.get("_orig_feat_type", "")
                if not feat_type and hasattr(feat, "feature_type"):
                    # For Gene/Transcript which have feature_type property
                    feat_type = feat.feature_type
                feat_type = feat_type.lower() if feat_type else "unknown"

                parent_name = getattr(feat, "parent_id", None) or feat.attributes.get(
                    "_parent_name"
                )
                grandparent_name = feat.attributes.get("_grandparent_name")

                # For child features (exon/cds), create unique name per parent
                # This handles cases where the same exon is shared across multiple transcripts
                if feat_type in self.child_features:
                    # Create unique name: feature_type:parent:coords
                    name = f"{feat_type}_Parent={parent_name}:{feat.start}_{feat.stop}"
                else:
                    name = feat.feature_id

                # Check for duplicate coordinates per parent for child features
                # (same exon coordinates appearing twice in the same transcript)
                if feat_type in self.child_features and parent_name:
                    check_coords = (feat_type, feat.start, feat.stop)
                    if check_coords in unique_coords[parent_name]:
                        # Skip duplicate coordinates for this parent
                        continue
                    unique_coords[parent_name].add(check_coords)

                # Add feature to index
                feat_index[name] = feat

                # Build graph edges
                if parent_name is not None:
                    feat_graph.add_edge(parent_name, name)

                    # Ensure parent exists in index
                    parent_class = self.parent_class_map["parent"].get(feat_type)
                    if parent_class and parent_name not in feat_index:
                        # Create minimal placeholder parent with dummy coordinates
                        dummy_coord = GenomicCoordinate(
                            chromosome="unknown",
                            start=1,
                            stop=2,  # stop must be > start
                            strand="+",
                            system="1-based",
                        )
                        if parent_class == Gene:
                            feat_index[parent_name] = Gene(
                                feature_id=parent_name, coordinates=dummy_coord
                            )
                        else:
                            # For placeholder transcripts without a parent gene,
                            # set parent_id to own ID (matches original intronIC behavior)
                            feat_index[parent_name] = Transcript(
                                feature_id=parent_name,
                                coordinates=dummy_coord,
                                parent_id=parent_name,  # Use self as parent when gene is missing
                            )

                    # Ensure grandparent exists if specified
                    if grandparent_name is not None:
                        feat_graph.add_edge(grandparent_name, parent_name)
                        grandparent_class = self.parent_class_map["grandparent"].get(
                            feat_type
                        )
                        if grandparent_class and grandparent_name not in feat_index:
                            # Create minimal placeholder grandparent
                            dummy_coord = GenomicCoordinate(
                                chromosome="unknown",
                                start=1,
                                stop=2,  # stop must be > start
                                strand="+",
                                system="1-based",
                            )
                            feat_index[grandparent_name] = Gene(
                                feature_id=grandparent_name, coordinates=dummy_coord
                            )

        # Step 1b: Create placeholder genes for transcripts with missing parent genes
        # This handles cases like C_gene_segment -> gene where the gene feature
        # exists in the annotation but wasn't created (because parent=[])
        dummy_coord = GenomicCoordinate(
            chromosome="unknown",
            start=1,
            stop=2,
            strand="+",
            system="1-based",
        )
        for name, feat in list(feat_index.items()):
            if isinstance(feat, Transcript):
                parent_id = getattr(feat, "parent_id", None)
                if parent_id and parent_id not in feat_index:
                    # Create placeholder Gene for this transcript's parent
                    feat_index[parent_id] = Gene(
                        feature_id=parent_id, coordinates=dummy_coord
                    )
                    feat_graph.add_edge(parent_id, name)

        # Step 2: Remove cycles if any
        self._remove_cycles(feat_graph)

        # Step 3: Build parent-child relationships via topological sort
        top_level_names = self._build_relationships(feat_graph, feat_index)

        # Step 4: Extract top-level features (genes)
        top_level = [feat_index[name] for name in top_level_names if name in feat_index]

        # For per-contig processing, empty contigs are valid (e.g., scaffolds with no genes)
        # Only raise error if we received annotations but couldn't build any hierarchy
        if not top_level and feat_index:
            raise ValueError(
                "Could not establish parent-child relationships. "
                "Check annotation file format."
            )
        elif not top_level:
            # Empty contig (no features) - return empty list
            return []

        # Step 5: Wrap orphan transcripts in genes (for annotations missing gene features)
        # When gene features are missing, transcripts become top-level but need Gene wrappers
        from dataclasses import replace

        wrapped_top_level = []
        for feature in top_level:
            if isinstance(feature, Transcript):
                # Transcript has no parent gene - create a Gene wrapper with same ID
                # Use a unique key by prefixing with "gene_wrapper:"
                gene_id = f"gene_wrapper:{feature.feature_id}"
                if gene_id not in feat_index:
                    gene = Gene(
                        feature_id=gene_id,
                        coordinates=feature.coordinates,
                        attributes={"_orig_feat_type": "gene", "_synthetic": True},
                    )
                    # Create updated transcript with parent_id pointing to gene
                    updated_transcript = replace(feature, parent_id=gene_id)

                    # Update feature_index with new transcript
                    feat_index[feature.feature_id] = updated_transcript

                    # Add transcript to gene's children
                    gene.add_child(feature.feature_id)

                    # Add gene to index
                    feat_index[gene_id] = gene
                    wrapped_top_level.append(gene)
                else:
                    # Shouldn't happen, but handle gracefully
                    wrapped_top_level.append(feature)
            else:
                # Already a Gene
                wrapped_top_level.append(feature)

        top_level = wrapped_top_level

        # Step 6: Set family sizes and coding lengths
        # TODO: These methods don't exist in the new model yet
        # Will be implemented in a later phase if needed
        # for gene in top_level:
        #     if isinstance(gene, Gene):
        #         gene.set_family_size()
        #         gene.coding_length = gene.get_coding_length()

        # Store feature index for later use
        self.feature_index = feat_index

        return top_level

    def _create_features_from_annotation(
        self, ann: AnnotationLine
    ) -> List[Gene | Transcript | Exon]:
        """
        Create feature objects from an annotation line.

        An annotation line may have multiple parents, so this returns
        a list of features (one per parent).

        Args:
            ann: Parsed annotation line

        Returns:
            List of feature objects (Gene, Transcript, or Exon)
        """
        feat_type = ann.feat_type.lower()

        # Map feature types to classes
        feature_class_map = {
            "gene": Gene,
            "transcript": Transcript,
            "mrna": Transcript,  # Common alias
            "exon": Exon,
            "cds": Exon,
        }

        # Default to Transcript for unknown types
        feature_class = feature_class_map.get(feat_type, Transcript)

        # Create coordinate once
        coord = self._create_coordinate(ann)

        features = []
        # Handle multiple parents (some annotations have this)
        parents = ann.parent if isinstance(ann.parent, list) else [ann.parent]

        for parent in parents:
            # Skip if both parent and grandparent are the same (annotation error)
            grandparent = ann.grandparent
            if grandparent == parent:
                grandparent = None

            # Create attributes dict with metadata needed during hierarchy building
            attributes = {
                "_parent_name": parent,  # Temporary for graph building
                "_grandparent_name": grandparent,  # Temporary for graph building
                "_line_number": ann.line_number,  # Temporary for ordering
                "_orig_feat_type": feat_type,  # Store original feature type
            }

            # Create feature instance
            if feature_class == Exon:
                # Parse phase
                phase = None
                if ann.phase is not None and ann.phase != ".":
                    try:
                        phase = int(ann.phase)
                    except (ValueError, TypeError):
                        phase = None

                feat = Exon(
                    feature_id=ann.name,
                    coordinates=coord,
                    parent_id=parent,
                    attributes=attributes,
                    phase=phase,
                    is_coding=(feat_type == "cds"),
                )
            elif feature_class == Transcript:
                feat = Transcript(
                    feature_id=ann.name,
                    coordinates=coord,
                    parent_id=parent,
                    attributes=attributes,
                )
            else:  # Gene
                feat = Gene(
                    feature_id=ann.name, coordinates=coord, attributes=attributes
                )

            features.append(feat)

        return features

    def _remove_cycles(self, graph: DiGraph) -> None:
        """
        Detect and remove cycles from the feature graph.

        Optimized approach:
        1. First check if graph is acyclic (fast check)
        2. Only iterate if cycles exist
        3. find_cycle() returns one cycle at a time, so we loop until none remain

        Args:
            graph: NetworkX directed graph of feature relationships
        """
        # Fast check: if already acyclic, skip expensive cycle finding
        if is_directed_acyclic_graph(graph):
            return

        all_cycle_edges = []
        max_iterations = 1000  # Safety limit to prevent infinite loops

        for iteration in range(max_iterations):
            try:
                cycles = list(find_cycle(graph))
                if cycles:
                    # Store cycle information for debugging
                    all_cycle_edges.extend(cycles)
                    graph.remove_edges_from(cycles)
                else:
                    # No more cycles found
                    break
            except:
                # No cycles found (find_cycle raises exception when graph is acyclic)
                break

        if all_cycle_edges:
            # Log cycle warning at warning level (always visible)
            msg = f"Removed {len(all_cycle_edges)} cycle edge(s) from feature graph"
            if self.messenger:
                self.messenger.warning(msg)
            else:
                print(f"[!] Warning: {msg}")

            # Log detailed cycle information at debug level (only in debug mode)
            # Extract unique features involved in cycles
            features_in_cycles = set()
            cycle_details = []
            for parent, child in all_cycle_edges:
                features_in_cycles.add(parent)
                features_in_cycles.add(child)
                cycle_details.append(f"  Edge removed: {parent} -> {child}")

            debug_msg = "Cycle details:\n" + "\n".join(cycle_details)
            debug_msg += (
                f"\nFeatures involved in cycles ({len(features_in_cycles)} total):\n"
            )
            debug_msg += "\n".join(f"  {f}" for f in sorted(features_in_cycles))

            if self.messenger:
                self.messenger.log_only(debug_msg, level="debug")
            else:
                # Fallback to print if no messenger
                print(f"[!] {debug_msg}")

    def _build_relationships(
        self, graph: DiGraph, feat_index: Dict[str, Gene | Transcript | Exon]
    ) -> Set[str]:
        """
        Build parent-child relationships using topological sort.

        This ensures parents are processed before their children.

        Args:
            graph: NetworkX directed graph of features
            feat_index: Dictionary mapping names to feature objects

        Returns:
            Set of top-level feature names (usually genes)
        """
        top_level_names = set()

        # Process features in topological order (parents before children)
        for name in lexicographical_topological_sort(graph):
            if name not in feat_index:
                continue

            feat = feat_index[name]
            children = list(graph.successors(name))
            parents = list(graph.predecessors(name))

            # Track top-level features (have children but no parents)
            if children and not parents:
                top_level_names.add(name)
                continue

            # Fix missing grandparent for child features
            feat_type = feat.attributes.get("_orig_feat_type", "")
            if not feat_type and hasattr(feat, "feature_type"):
                feat_type = feat.feature_type
            feat_type = feat_type.lower() if feat_type else "unknown"

            grandparent_name = feat.attributes.get("_grandparent_name")
            parent_name = getattr(feat, "parent_id", None) or feat.attributes.get(
                "_parent_name"
            )

            if (
                feat_type in self.child_features
                and grandparent_name is None
                and parent_name
                and parent_name in feat_index
            ):
                parent_feat = feat_index[parent_name]
                gp_name = getattr(
                    parent_feat, "parent_id", None
                ) or parent_feat.attributes.get("_parent_name")
                if gp_name is None:
                    # When grandparent is missing (e.g., no gene features in annotation),
                    # use parent ID as grandparent ID (matches original intronIC behavior)
                    gp_name = parent_name
                # Update grandparent in attributes dict (which is mutable even though Exon is frozen)
                if hasattr(feat, "attributes") and isinstance(feat.attributes, dict):
                    feat.attributes["_grandparent_name"] = gp_name

            # Add this feature to its parents' children lists
            # Use 'name' (the unique name from the graph) not feat.feature_id
            for parent_name in parents:
                if parent_name in feat_index:
                    parent_feat = feat_index[parent_name]
                    if isinstance(parent_feat, (Gene, Transcript)):
                        parent_feat.add_child(name)  # Use graph name, not feature_id

        return top_level_names


def build_annotation_hierarchy(
    annotation_file: str, child_features: List[str]
) -> List[Gene]:
    """
    Convenience function to build gene hierarchy from an annotation file.

    This is a functional wrapper around AnnotationHierarchyBuilder for
    backwards compatibility with the original intronIC API.

    Args:
        annotation_file: Path to GFF3/GTF file
        child_features: List of feature types to extract (e.g., ['exon', 'cds'])

    Returns:
        List of Gene objects with complete hierarchy

    Examples:
        >>> genes = build_annotation_hierarchy('chr19.gff3', ['exon', 'cds'])
        >>> print(f"Parsed {len(genes)} genes")
    """
    builder = AnnotationHierarchyBuilder(child_features)
    return builder.build_from_file(annotation_file)
