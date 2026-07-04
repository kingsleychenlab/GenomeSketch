"""Tests for motif search."""

from genomesketch.motifs import count_motif, find_motif, search_motif


def test_overlapping_matches():
    assert find_motif("AAAAA", "AAA") == [0, 1, 2]


def test_no_match():
    assert find_motif("ATGCGT", "TTT") == []


def test_case_insensitive():
    assert find_motif("atgatg", "ATG") == [0, 3]


def test_search_motif_forward_only():
    matches = search_motif("ATGCATG", "ATG")
    assert [(m.start, m.strand) for m in matches] == [(0, "+"), (4, "+")]


def test_revcomp_search():
    # Motif ATG, reverse complement CAT. Sequence contains ATG at 0 and CAT at 4.
    matches = search_motif("ATGGCAT", "ATG", revcomp=True)
    strands = {(m.start, m.strand) for m in matches}
    assert (0, "+") in strands
    assert (4, "-") in strands


def test_count_motif():
    assert count_motif("AAAAA", "AA") == 4
