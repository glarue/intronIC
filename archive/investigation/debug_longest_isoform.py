#!/usr/bin/env python3
"""
Debug longest isoform identification.

Test the sorting logic with a sample of introns.
"""

# Simulate some introns with different parent lengths
class MockMetadata:
    def __init__(self, defined_by, parent_length, parent, grandparent, family_size, index):
        self.defined_by = defined_by
        self.parent_length = parent_length
        self.parent = parent
        self.grandparent = grandparent
        self.family_size = family_size
        self.index = index
        self.longest_isoform = False

class MockIntron:
    def __init__(self, name, metadata):
        self.name = name
        self.metadata = metadata

# Create test introns for gene1 with 3 transcripts of different lengths
introns = [
    # Gene1, Transcript A (length 1000), 2 introns
    MockIntron("Gene1_TxA_int1", MockMetadata("cds", 1000, "TxA", "Gene1", 2, 1)),
    MockIntron("Gene1_TxA_int2", MockMetadata("cds", 1000, "TxA", "Gene1", 2, 2)),

    # Gene1, Transcript B (length 1500 - LONGEST), 3 introns
    MockIntron("Gene1_TxB_int1", MockMetadata("cds", 1500, "TxB", "Gene1", 3, 1)),
    MockIntron("Gene1_TxB_int2", MockMetadata("cds", 1500, "TxB", "Gene1", 3, 2)),
    MockIntron("Gene1_TxB_int3", MockMetadata("cds", 1500, "TxB", "Gene1", 3, 3)),

    # Gene1, Transcript C (length 800), 1 intron
    MockIntron("Gene1_TxC_int1", MockMetadata("cds", 800, "TxC", "Gene1", 1, 1)),

    # Gene2, Transcript D (length 2000), 2 introns
    MockIntron("Gene2_TxD_int1", MockMetadata("exon", 2000, "TxD", "Gene2", 2, 1)),
    MockIntron("Gene2_TxD_int2", MockMetadata("exon", 2000, "TxD", "Gene2", 2, 2)),
]

print("Original order:")
for i in introns:
    print(f"  {i.name}: parent_length={i.metadata.parent_length}, parent={i.metadata.parent}")

# Sort using the refactored logic
sorted_introns = sorted(
    introns,
    key=lambda i: (
        i.metadata.defined_by or '',           # CDS before exon
        -(i.metadata.parent_length or 0),     # Descending by parent length
        i.metadata.parent or '',               # Transcript ID
        -(i.metadata.family_size or 0),       # Descending by family size
        i.metadata.index or 0                  # Intron index
    )
)

print("\nSorted order (by hierarchical_sort_attrs):")
for i in sorted_introns:
    print(f"  {i.name}: defined_by={i.metadata.defined_by}, parent_length={i.metadata.parent_length}, parent={i.metadata.parent}")

# Identify longest isoforms
longest_isoforms = {}
for intron in sorted_introns:
    grandparent = intron.metadata.grandparent
    parent = intron.metadata.parent

    if grandparent:
        if grandparent not in longest_isoforms:
            longest_isoforms[grandparent] = parent
            intron.metadata.longest_isoform = True
            print(f"\n  First intron for {grandparent}: {parent} (length={intron.metadata.parent_length}) marked as longest")
        elif longest_isoforms[grandparent] == parent:
            intron.metadata.longest_isoform = True
        else:
            intron.metadata.longest_isoform = False

print("\nLongest isoforms identified:")
for gene, transcript in longest_isoforms.items():
    print(f"  {gene}: {transcript}")

print("\nFinal intron tagging:")
for i in introns:
    status = "LONGEST" if i.metadata.longest_isoform else "not longest"
    print(f"  {i.name}: {status}")

print("\nExpected result:")
print("  Gene1: TxB should be longest (length 1500)")
print("  Gene2: TxD should be longest (length 2000)")
print("  All introns from TxB and TxD should be marked as longest_isoform=True")
