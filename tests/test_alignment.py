"""Tests for Needleman-Wunsch (and Smith-Waterman) alignment."""

import pytest

from genomesketch.alignment import needleman_wunsch, smith_waterman


def test_identical_sequences_score():
    res = needleman_wunsch("ATGC", "ATGC")
    assert res.score == 4
    assert res.identity == pytest.approx(100.0)
    assert res.gaps == 0
    assert res.aligned1 == "ATGC"
    assert res.aligned2 == "ATGC"


def test_single_gap_alignment():
    res = needleman_wunsch("ATGCGT", "ATGACGT")
    # 6 matches + 1 gap: score = 6*1 + 1*(-2) = 4
    assert res.score == 4
    assert len(res.aligned1) == len(res.aligned2)
    assert res.gaps == 1
    assert "-" in res.aligned1 or "-" in res.aligned2


def test_traceback_reconstructs_sequences():
    res = needleman_wunsch("ATGCGT", "ATGACGT")
    # Removing gaps must recover the originals.
    assert res.aligned1.replace("-", "") == "ATGCGT"
    assert res.aligned2.replace("-", "") == "ATGACGT"


def test_all_mismatch_score():
    res = needleman_wunsch("AAAA", "TTTT")
    # Best is 4 mismatches (score -4) vs gaps (worse).
    assert res.score == -4
    assert res.mismatches == 4


def test_custom_scoring():
    res = needleman_wunsch("ATGC", "ATGC", match=2)
    assert res.score == 8


def test_match_line_annotations():
    res = needleman_wunsch("ATGC", "ATGA")
    # 3 matches, 1 mismatch, no gaps -> '|||.'
    assert res.match == "|||."


def test_smith_waterman_local():
    # Local alignment should find the common core "GATTACA".
    res = smith_waterman("XXXGATTACAXXX", "YYGATTACAYY")
    assert "GATTACA" in res.aligned1
    assert res.score == 7  # 7 exact matches
    assert res.match_type == "local"
