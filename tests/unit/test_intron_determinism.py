"""
Regression tests for deterministic intron extraction.

intronIC stores a transcript's exon IDs in a ``set`` (``Parent.children``), so
iterating them yields a string-hash-ordered sequence that varies with
``PYTHONHASHSEED`` from one process to the next. Intron generation sorts those
exons by genomic ``start`` — but that is NOT a total order when exons overlap
(an annotation artifact) and share a start coordinate. A stable sort then falls
back to the hash-dependent input order, so *which* overlapping-exon intron gets
skipped changed between runs. Extraction must instead be fully deterministic:
the exon sort breaks ties by ``_line_number`` (annotation/parse order), which is
unique and stable, so the emitted intron set is independent of input order.

These tests exercise the ordering directly with same-start exons, which is where
the non-determinism lived; clean (all-distinct-start) transcripts were never
affected because the tiebreak is never consulted.
"""

import pytest

from intronIC.core.models import Exon
from intronIC.extraction.intronator import IntronGenerator
from intronIC.utils.coordinates import GenomicCoordinate

pytestmark = pytest.mark.unit


def _exon(feature_id, start, stop, strand, line_number, parent="tx1"):
    return Exon(
        feature_id=feature_id,
        coordinates=GenomicCoordinate("chr1", start, stop, strand, "1-based"),
        parent_id=parent,
        attributes={"_line_number": line_number},
        phase=None,
        is_coding=False,
    )


def _perms(exons):
    """A few deterministic re-orderings of the same exon list."""
    return [
        list(exons),
        list(reversed(exons)),
        [exons[i] for i in (2, 0, 3, 1)] if len(exons) == 4 else list(reversed(exons)),
    ]


class TestExonSortDeterminism:
    def test_same_start_sort_is_input_order_invariant_plus_strand(self):
        """Overlapping exons sharing a start must sort by _line_number, not input order."""
        gen = IntronGenerator()
        # e2/e3 share start=100 (overlap); tie must resolve by _line_number (e2 before e3).
        e1 = _exon("e1", 10, 50, "+", line_number=1)
        e2 = _exon("e2", 100, 200, "+", line_number=2)
        e3 = _exon("e3", 100, 400, "+", line_number=3)
        e4 = _exon("e4", 500, 600, "+", line_number=4)
        base = [e1, e2, e3, e4]

        orders = [gen._sort_exons_by_coding_direction(p) for p in _perms(base)]
        id_orders = [[e.feature_id for e in o] for o in orders]
        assert all(o == id_orders[0] for o in id_orders), id_orders
        # ascending start, ties by ascending line_number
        assert id_orders[0] == ["e1", "e2", "e3", "e4"]

    def test_same_start_sort_is_input_order_invariant_minus_strand(self):
        """Negative strand: descending start, ties still ascending by _line_number."""
        gen = IntronGenerator()
        e1 = _exon("e1", 10, 50, "-", line_number=1)
        e2 = _exon("e2", 100, 200, "-", line_number=2)
        e3 = _exon("e3", 100, 400, "-", line_number=3)
        e4 = _exon("e4", 500, 600, "-", line_number=4)
        base = [e1, e2, e3, e4]

        orders = [gen._sort_exons_by_coding_direction(p) for p in _perms(base)]
        id_orders = [[e.feature_id for e in o] for o in orders]
        assert all(o == id_orders[0] for o in id_orders), id_orders
        # coding order = high coord first (e4, then the start=100 pair by line order, then e1)
        assert id_orders[0] == ["e4", "e2", "e3", "e1"]


class TestIntronGenerationDeterminism:
    def test_intron_set_invariant_to_exon_input_order(self):
        """The emitted intron coordinate set must not depend on input exon order."""
        gen = IntronGenerator()
        # Clean, non-overlapping exons -> 3 introns; order must never matter.
        exons = [
            _exon("e1", 10, 100, "+", 1),
            _exon("e2", 200, 300, "+", 2),
            _exon("e3", 400, 500, "+", 3),
            _exon("e4", 600, 700, "+", 4),
        ]
        sets = []
        for perm in _perms(exons):
            introns = list(gen.generate_from_exons(perm))
            sets.append(sorted((i.coordinates.start, i.coordinates.stop) for i in introns))
        assert all(s == sets[0] for s in sets), sets
        assert sets[0] == [(101, 199), (301, 399), (501, 599)]

    def test_overlapping_exon_intron_set_is_deterministic(self):
        """With a same-start overlap, the surviving intron set is fixed regardless of order."""
        gen = IntronGenerator()
        # e2/e3 overlap (share start=200). Whichever is chosen, the OUTPUT must be identical
        # across input orderings — that is the property that was violated by hash ordering.
        exons = [
            _exon("e1", 10, 100, "+", 1),
            _exon("e2", 200, 300, "+", 2),
            _exon("e3", 200, 500, "+", 3),
            _exon("e4", 600, 700, "+", 4),
        ]
        sets = []
        for perm in _perms(exons):
            introns = list(gen.generate_from_exons(perm))
            sets.append(sorted((i.coordinates.start, i.coordinates.stop) for i in introns))
        assert all(s == sets[0] for s in sets), sets
