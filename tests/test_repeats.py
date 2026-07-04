"""Tests for repeat detection."""

from genomesketch.repeats import (
    find_repeats,
    find_tandem_repeats,
    kmer_positions,
)


def test_kmer_positions_overlapping():
    positions = kmer_positions("ATGATG", 3)
    assert positions["ATG"] == [0, 3]
    assert positions["TGA"] == [1]


def test_find_repeats_basic():
    repeats = find_repeats("ATGATGATGCC", 3, min_count=3)
    top = repeats[0]
    assert top.motif == "ATG"
    assert top.copies == 3
    assert top.positions == [0, 3, 6]
    assert top.region == (0, 9)


def test_find_repeats_min_count_filter():
    repeats = find_repeats("ATGCC", 3, min_count=2)
    assert repeats == []


def test_tandem_repeats():
    tandems = find_tandem_repeats("ATGATGATGCC", 3, min_count=3)
    assert len(tandems) == 1
    t = tandems[0]
    assert t.motif == "ATG"
    assert t.copies == 3
    assert t.start == 0
    assert t.region == (0, 9)


def test_tandem_requires_adjacency():
    # ATG then a gap then ATG -> not a tandem run of 2.
    tandems = find_tandem_repeats("ATGCCCATG", 3, min_count=2)
    assert tandems == []
