"""Tests for nucleotide statistics."""

import math

from genomesketch.stats import (
    at_skew,
    compute_stats,
    gc_content,
    gc_skew,
    gc_windows,
    nucleotide_counts,
)


def test_counts():
    counts = nucleotide_counts("AATTCCGGN")
    assert counts == {"A": 2, "T": 2, "C": 2, "G": 2, "N": 1}


def test_gc_content_basic():
    # 2 G/C out of 4 known bases -> 0.5
    assert gc_content("ATGC") == 0.5


def test_gc_content_excludes_n():
    # GC over A+T+G+C, N excluded from denominator.
    assert gc_content("GCNN") == 1.0


def test_gc_content_zero_division_safe():
    assert gc_content("") == 0.0
    assert gc_content("NNNN") == 0.0


def test_gc_skew():
    # G=3, C=1 -> (3-1)/(3+1) = 0.5
    assert gc_skew("GGGC") == 0.5


def test_gc_skew_zero_division_safe():
    assert gc_skew("ATATAT") == 0.0


def test_at_skew():
    # A=1, T=3 -> (1-3)/(1+3) = -0.5
    assert at_skew("ATTT") == -0.5


def test_compute_stats_bundle():
    s = compute_stats("ATGCGTACCGTTAGCTAGCTAGGCTA")
    assert s.length == 26
    assert s.a + s.t + s.c + s.g + s.n == 26
    assert math.isclose(s.gc_content, gc_content("ATGCGTACCGTTAGCTAGCTAGGCTA"))


def test_gc_windows_shape():
    windows = gc_windows("ATGCATGCATGC", window=4, step=4)
    starts = [w[0] for w in windows]
    assert starts == [0, 4, 8]
    for _, gc in windows:
        assert 0.0 <= gc <= 1.0


def test_gc_windows_includes_partial_tail():
    windows = gc_windows("ATGCA", window=4, step=4)
    # 0..4 then partial tail starting at 4
    assert windows[0][0] == 0
    assert windows[-1][0] == 4
