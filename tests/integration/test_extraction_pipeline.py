"""
Integration tests for the extraction pipeline.

Tests the full pipeline from annotation to introns.
"""

from pathlib import Path
from statistics import median

import pytest

from intronIC.extraction.annotator import AnnotationHierarchyBuilder
from intronIC.extraction.filters import IntronFilter
from intronIC.extraction.intronator import IntronGenerator
from intronIC.extraction.sequences import SequenceExtractor

# Mark all tests in this module
pytestmark = [
    pytest.mark.integration,
    pytest.mark.extraction,
    pytest.mark.requires_chr19,
]

# Test data paths - use src layout path
TEST_DATA_DIR = (
    Path(__file__).parent.parent.parent / "src" / "intronIC" / "data" / "test_data"
)
CHR19_ANNOTATION = TEST_DATA_DIR / "Homo_sapiens.Chr19.Ensembl_91.gff3.gz"
CHR19_GENOME = TEST_DATA_DIR / "Homo_sapiens.Chr19.Ensembl_91.fa.gz"


@pytest.mark.skipif(
    not CHR19_ANNOTATION.exists(), reason="Chr19 test data not available"
)
class TestExtractionPipeline:
    """Test the full extraction pipeline."""

    def test_annotation_hierarchy_building(self):
        """Test building gene hierarchy from chr19 annotation."""
        builder = AnnotationHierarchyBuilder(["exon", "cds"])
        genes = builder.build_from_file(str(CHR19_ANNOTATION))

        # Should have genes
        assert len(genes) > 0
        print(f"Found {len(genes)} genes")

        # Check structure
        first_gene = genes[0]
        assert hasattr(first_gene, "children")
        assert hasattr(first_gene, "feature_id")

        # Should have transcripts as children
        transcript_count = sum(
            1
            for child_id in first_gene.children
            if child_id in builder.feature_index
            and hasattr(builder.feature_index[child_id], "children")
        )
        assert transcript_count > 0
        print(f"First gene has {transcript_count} transcripts")

    def test_intron_generation_from_genes(self):
        """Test generating introns from gene hierarchy."""
        # Build hierarchy
        builder = AnnotationHierarchyBuilder(["exon", "cds"])
        genes = builder.build_from_file(str(CHR19_ANNOTATION))

        # Generate introns
        generator = IntronGenerator()
        introns = list(generator.generate_from_genes(genes, builder.feature_index))

        # Should have many introns
        assert len(introns) > 1000  # Chr19 should have thousands
        print(f"Generated {len(introns)} introns")

        # Check intron structure
        first_intron = introns[0]
        assert hasattr(first_intron, "coordinates")
        assert hasattr(first_intron, "metadata")
        assert first_intron.coordinates.start < first_intron.coordinates.stop

    def test_sequence_extraction(self):
        """Test extracting sequences for introns."""
        # Build hierarchy
        builder = AnnotationHierarchyBuilder(["exon"])
        genes = builder.build_from_file(str(CHR19_ANNOTATION))

        # Generate introns (just first 10 for speed)
        generator = IntronGenerator()
        introns = list(generator.generate_from_genes(genes, builder.feature_index))[:10]

        # Extract sequences
        extractor = SequenceExtractor(str(CHR19_GENOME))
        introns_with_seqs = list(extractor.extract_sequences(introns))

        # Check sequences
        assert len(introns_with_seqs) == 10
        for intron in introns_with_seqs:
            assert intron.sequences is not None
            assert intron.sequences.seq
            assert len(intron.sequences.seq) > 0
            print(f"Intron {intron.intron_id}: {len(intron.sequences.seq)} bp")

    def test_full_pipeline_intron_count(self):
        """Test full pipeline and compare intron count to gold standard."""
        # Expected count from gold standard
        EXPECTED_INTRON_COUNT = 20252

        # Step 1: Build hierarchy
        print("Building annotation hierarchy...")
        builder = AnnotationHierarchyBuilder(["exon", "cds"])
        genes = builder.build_from_file(str(CHR19_ANNOTATION))
        print(f"  Found {len(genes)} genes")

        # Step 2: Generate introns
        print("Generating introns...")
        generator = IntronGenerator()
        introns = list(generator.generate_from_genes(genes, builder.feature_index))
        print(f"  Generated {len(introns)} introns")

        # Step 3: Filter introns (without sequences for speed)
        print("Filtering introns...")
        filter_obj = IntronFilter(
            min_length=30,
            allow_noncanonical=True,  # Include all for now
            allow_overlap=True,
            longest_only=False,
            include_duplicates=False,
        )

        # Just check omission without full filtering for speed
        omitted_count = 0
        duplicate_count = 0
        for intron in introns:
            filter_obj._check_omission(intron)
            filter_obj._tag_intron(intron)

            if intron.metadata.omitted:
                omitted_count += 1
            if intron.metadata.duplicate:
                duplicate_count += 1

        unique_introns = len(introns) - duplicate_count
        kept_introns = unique_introns - omitted_count

        print(f"  Total introns: {len(introns)}")
        print(f"  Duplicates: {duplicate_count}")
        print(f"  Unique introns: {unique_introns}")
        print(f"  Omitted: {omitted_count}")
        print(f"  Kept: {kept_introns}")

        # The count should be close to expected (within 10%)
        # We'll refine this as we verify against gold standard
        assert unique_introns > EXPECTED_INTRON_COUNT * 0.9
        assert unique_introns < EXPECTED_INTRON_COUNT * 1.1
        print(
            f"✓ Intron count within expected range (expected: {EXPECTED_INTRON_COUNT})"
        )

    def test_filtering_omission_codes(self):
        """Test that omission codes are correctly assigned."""
        # Build hierarchy
        builder = AnnotationHierarchyBuilder(["exon"])
        genes = builder.build_from_file(str(CHR19_ANNOTATION))

        # Generate introns
        generator = IntronGenerator()
        introns = list(generator.generate_from_genes(genes, builder.feature_index))[
            :100
        ]

        # Extract sequences for proper filtering
        extractor = SequenceExtractor(str(CHR19_GENOME))
        introns_with_seqs = list(extractor.extract_sequences(introns))

        # Filter
        filter_obj = IntronFilter(
            min_length=100,  # Higher threshold to catch some
            allow_noncanonical=False,
            allow_overlap=False,
            longest_only=True,
        )

        for intron in introns_with_seqs:
            filter_obj._check_omission(intron)

        # Check for various omission codes
        omission_codes = set()
        for intron in introns_with_seqs:
            if intron.metadata.omitted:
                omission_codes.add(intron.metadata.omitted)

        print(f"Found omission codes: {omission_codes}")
        # Should have at least some omissions
        assert len(omission_codes) > 0


@pytest.mark.skipif(
    not CHR19_ANNOTATION.exists(), reason="Chr19 test data not available"
)
def test_quick_extraction():
    """Quick test for development - just make sure nothing crashes."""
    builder = AnnotationHierarchyBuilder(["exon"])
    genes = builder.build_from_file(str(CHR19_ANNOTATION))

    generator = IntronGenerator()
    introns = list(generator.generate_from_genes(genes, builder.feature_index))

    print(f"✓ Generated {len(introns)} introns from chr19")
    assert len(introns) > 10000  # Should be thousands


@pytest.mark.skipif(
    not CHR19_ANNOTATION.exists(), reason="Chr19 test data not available"
)
def test_fractional_position_computed():
    """Regression guard for the frac_pos bug (all introns were 0.0) and for the CDS-relative
    definition: fractional_position is (CDS length upstream of the intron) / (total CDS
    length), so CDS introns lie in [0, 1], are non-degenerate, and are monotonic within a
    transcript, while exon-only (UTR) introns have no CDS position and are None.
    See intronator.generate_from_transcript."""
    from collections import defaultdict

    builder = AnnotationHierarchyBuilder(["exon", "cds"])
    genes = builder.build_from_file(str(CHR19_ANNOTATION))
    introns = list(IntronGenerator().generate_from_genes(genes, builder.feature_index))

    cds_fps = [i.metadata.fractional_position for i in introns
               if i.metadata.defined_by == "cds"]
    exon_fps = [i.metadata.fractional_position for i in introns
                if i.metadata.defined_by != "cds"]

    # CDS introns carry a CDS-relative position; exon-only (UTR) introns are None.
    assert all(f is not None for f in cds_fps), "CDS introns must have a frac_pos"
    assert all(f is None for f in exon_fps), "exon-only (UTR) introns must have frac_pos None"
    # Not the old degenerate all-zero output.
    nonzero = [f for f in cds_fps if f > 0.0]
    assert len(nonzero) > 0.9 * len(cds_fps), "frac_pos is degenerate (mostly 0.0)"
    assert len(set(cds_fps)) > 1000, "frac_pos should take many distinct values"
    assert all(0.0 <= f <= 1.0 for f in cds_fps), "frac_pos must lie in [0, 1]"
    # CDS-normalization removes the 3'UTR-driven 5' skew: the CDS-intron distribution
    # should be centered near the middle of the coding sequence, not bunched at 0.
    assert 0.4 < median(cds_fps) < 0.6, f"CDS frac_pos median off-center: {median(cds_fps):.3f}"

    by_tx = defaultdict(list)
    for i in introns:
        by_tx[i.metadata.parent].append(i)

    # Within any transcript, the CDS introns' frac_pos must be non-decreasing with ordinal
    # index (exon-only introns are None and are skipped).
    checked = 0
    for members in by_tx.values():
        cds_members = [m for m in members if m.metadata.defined_by == "cds"]
        if len(cds_members) < 3:
            continue
        ordered = sorted(cds_members, key=lambda x: x.metadata.index)
        seq = [x.metadata.fractional_position for x in ordered]
        assert seq == sorted(seq), f"frac_pos not monotonic in {ordered[0].metadata.parent}: {seq}"
        checked += 1
    assert checked > 100, "expected many multi-intron transcripts to check"

    # The chr19 set genuinely contains exon-only (UTR) introns, so the "None for UTR
    # introns" assertion above is exercised, not vacuous. frac_pos is CDS-relative and
    # therefore undefined for these; they are not placed on a shared transcript axis.
    assert len(exon_fps) > 0, "expected some exon-only (UTR) introns in the chr19 set"
