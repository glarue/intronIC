"""
Demonstration of the Tag → Attribute mapping system.

This shows how compact bracket tags in intron names map to verbose
attributes in the attributes column.

Author: intronIC refactoring project
Date: 2025-11-02
"""

from pathlib import Path
from file_io.writers import BEDWriter, MetaWriter, TAG_TO_ATTRIBUTE, generate_attributes
from core.intron import Intron, IntronMetadata
from utils.coordinates import GenomicCoordinate


def demonstrate_tag_mapping():
    """Demonstrate the tag to attribute mapping."""

    print("=" * 70)
    print("Tag → Attribute Mapping Dictionary")
    print("=" * 70)
    print("\nTAG_TO_ATTRIBUTE mapping:")
    for tag, attribute in TAG_TO_ATTRIBUTE.items():
        print(f"  [{tag}] → {attribute}")

    print("\n" + "=" * 70)
    print("Example Introns with Various Attributes")
    print("=" * 70)

    # Example 1: Noncanonical intron
    print("\n1. Noncanonical intron (AT-AC splice sites)")
    coord1 = GenomicCoordinate("chr1", 1001, 2000, '+', '1-based')
    meta1 = IntronMetadata(
        parent="TRANS001",
        index=1,
        family_size=5,
        noncanonical=True,
        longest_isoform=True
    )
    intron1 = Intron("intron1", coord1, metadata=meta1)
    print(f"   Tags in name: [n]")
    print(f"   Attributes:   {generate_attributes(intron1)}")

    # Example 2: Not longest isoform
    print("\n2. Intron from shorter isoform")
    meta2 = IntronMetadata(
        parent="TRANS002",
        index=2,
        family_size=3,
        longest_isoform=False
    )
    intron2 = Intron("intron2", coord1, metadata=meta2)
    print(f"   Tags in name: [i]")
    print(f"   Attributes:   {generate_attributes(intron2)}")

    # Example 3: Corrected boundaries
    print("\n3. Intron with corrected splice sites")
    meta3 = IntronMetadata(
        parent="TRANS003",
        index=1,
        family_size=4,
        corrected=True,
        longest_isoform=True
    )
    intron3 = Intron("intron3", coord1, metadata=meta3)
    print(f"   Tags in name: [c]")
    print(f"   Attributes:   {generate_attributes(intron3)}")

    # Example 4: Omitted (too short)
    print("\n4. Omitted intron (too short)")
    meta4 = IntronMetadata(
        parent="TRANS004",
        index=3,
        family_size=6,
        omitted='s',  # short
        longest_isoform=True
    )
    intron4 = Intron("intron4", coord1, metadata=meta4)
    print(f"   Tags in name: [o:s]")
    print(f"   Attributes:   {generate_attributes(intron4)}")

    # Example 5: Multiple attributes
    print("\n5. Intron with multiple attributes")
    meta5 = IntronMetadata(
        parent="TRANS005",
        index=2,
        family_size=8,
        noncanonical=True,
        longest_isoform=False,
        corrected=True,
        duplicate="intron_rep"
    )
    intron5 = Intron("intron5", coord1, metadata=meta5)
    print(f"   Tags in name: [n][i][c][d]")
    print(f"   Attributes:   {generate_attributes(intron5)}")

    # Example 6: Write to actual files
    print("\n" + "=" * 70)
    print("Example Output Files")
    print("=" * 70)

    # Create temp directory
    output_dir = Path("/tmp/intronic_demo")
    output_dir.mkdir(exist_ok=True)

    # Write BED file
    bed_file = output_dir / "demo.bed.iic"
    with BEDWriter(bed_file) as writer:
        writer.write_intron(intron1, species_name="demo")
        writer.write_intron(intron2, species_name="demo")
        writer.write_intron(intron5, species_name="demo")

    print(f"\nBED file written to: {bed_file}")
    print("Contents:")
    print(bed_file.read_text())

    # Write metadata file
    meta_file = output_dir / "demo.meta.iic"
    with MetaWriter(meta_file) as writer:
        writer.write_header()
        writer.write_intron(intron1, species_name="demo")
        writer.write_intron(intron2, species_name="demo")
        writer.write_intron(intron5, species_name="demo")

    print(f"Metadata file written to: {meta_file}")
    print("Contents (first 3 lines):")
    lines = meta_file.read_text().strip().split('\n')
    for line in lines[:3]:
        print(line)

    print("\n" + "=" * 70)
    print("Key Benefits of This System")
    print("=" * 70)
    print("""
1. **Portability**: Tags in name ensure compatibility with tools that
   only parse the name/ID column (bedtools, awk, etc.)

2. **Documentation**: Attributes column acts as a built-in legend,
   making it clear what [n][i][c] means

3. **Dual Access**: Both compact (tags) and verbose (attributes)
   available in the same file

4. **Maintainability**: Single TAG_TO_ATTRIBUTE dictionary defines
   the mapping for both writers

5. **Self-Explaining**: New users can understand tags by looking at
   the attributes column in the same row
""")


if __name__ == "__main__":
    demonstrate_tag_mapping()
