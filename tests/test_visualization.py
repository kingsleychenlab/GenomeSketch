"""Smoke tests: plots are generated as real, non-empty image files."""

import os

from genomesketch.alignment import needleman_wunsch
from genomesketch.motifs import search_motif
from genomesketch.orfs import find_orfs
from genomesketch.visualization import (
    plot_alignment_mismatches,
    plot_gc_content,
    plot_motif_map,
    plot_nucleotide_frequency,
    plot_orf_map,
)

SEQ = "ATGCGTACCGTTAGCTAGCTAGGCTAGCTAGCATGGGCCCTAAGGCATGAAATTTGGGCCCTAA"


def _nonempty(path: str) -> bool:
    return os.path.exists(path) and os.path.getsize(path) > 0


def test_plot_gc(tmp_path):
    out = str(tmp_path / "gc.png")
    plot_gc_content(SEQ, out, window=10, step=2)
    assert _nonempty(out)


def test_plot_bases(tmp_path):
    out = str(tmp_path / "bases.png")
    plot_nucleotide_frequency(SEQ, out)
    assert _nonempty(out)


def test_plot_motif(tmp_path):
    out = str(tmp_path / "motif.png")
    matches = search_motif(SEQ, "ATG", revcomp=True)
    plot_motif_map(len(SEQ), matches, out)
    assert _nonempty(out)


def test_plot_orfs(tmp_path):
    out = str(tmp_path / "orfs.png")
    orfs = find_orfs(SEQ)
    plot_orf_map(len(SEQ), orfs, out)
    assert _nonempty(out)


def test_plot_alignment(tmp_path):
    out = str(tmp_path / "aln.png")
    res = needleman_wunsch("ATGCGT", "ATGACGT")
    plot_alignment_mismatches(res.match, out)
    assert _nonempty(out)
