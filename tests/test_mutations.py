"""Tests for mutation comparison."""

import pytest

from genomesketch.mutations import (
    compare,
    compare_aligned,
    compare_equal_length,
)


def test_equal_length_substitutions():
    result = compare_equal_length("ATGCGT", "ATGAGT")
    assert len(result.substitutions) == 1
    m = result.substitutions[0]
    assert m.position == 3
    assert m.ref == "C"
    assert m.alt == "A"


def test_equal_length_identity():
    result = compare_equal_length("AAAA", "AATA")
    assert result.identity == pytest.approx(75.0)


def test_equal_length_mismatch_raises():
    with pytest.raises(ValueError):
        compare_equal_length("ATG", "ATGC")


def test_aligned_detects_insertion():
    # alt has an extra base relative to ref -> insertion.
    result = compare_aligned("ATGCGT", "ATGACGT")
    kinds = {m.kind for m in result.mutations}
    assert "insertion" in kinds


def test_compare_dispatch_equal_length():
    result = compare("ATGC", "ATGA")
    assert result.mode == "equal-length"


def test_compare_dispatch_aligned_when_unequal():
    result = compare("ATGC", "ATGCA")
    assert result.mode == "aligned"


def test_perfect_identity():
    result = compare("ATGCGT", "ATGCGT")
    assert result.identity == pytest.approx(100.0)
    assert result.mutations == []
