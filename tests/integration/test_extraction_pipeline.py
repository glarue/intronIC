"""
Integration tests for the extraction pipeline.

Tests the full pipeline from annotation to introns.
"""

from pathlib import Path

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
    assert len(introns) > 10000  # Should be thousands
