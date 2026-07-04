"""Open reading frame (ORF) detection across all six reading frames.

Conventions
-----------
* Start codon: ``ATG``. Stop codons: ``TAA``, ``TAG``, ``TGA``.
* All 6 frames are scanned (``+1, +2, +3, -1, -2, -3``).
* Coordinates are **0-based, half-open** on the *forward* strand. ``start`` is
  the first base of the start codon, ``stop`` is one past the last base of the
  stop codon, so ``stop - start == length`` (which includes the stop codon).
* The reported protein is the translation of the ORF *excluding* the trailing
  stop codon.

For each frame an ORF runs from the first ``ATG`` to the next in-frame stop
codon; scanning then resumes after that stop. This yields the longest ORF that
terminates at each stop codon.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List

from genomesketch.transform import reverse_complement
from genomesketch.translate import STOP_CODONS, translate
from genomesketch.validation import normalize

START_CODON = "ATG"


@dataclass
class ORF:
    """A single open reading frame (forward-strand coordinates)."""

    frame: int  # +1/+2/+3 or -1/-2/-3
    start: int  # 0-based, forward strand (first base of start codon)
    stop: int  # 0-based half-open, forward strand (one past stop codon)
    strand: str  # "+" or "-"
    protein: str  # translation, excluding the stop codon

    @property
    def length(self) -> int:
        """ORF length in base pairs, including the stop codon."""
        return self.stop - self.start

    def describe(self, index: int = 1) -> str:
        """Human-readable multi-line block (1-based label)."""
        return (
            f"ORF {index}\n"
            f"Frame: {self.frame:+d}\n"
            f"Start: {self.start}\n"
            f"Stop: {self.stop}\n"
            f"Length: {self.length} bp\n"
            f"Protein: {self.protein}"
        )


def _scan_frame(strand_seq: str, offset: int) -> List[tuple[int, int, str]]:
    """Scan one frame of ``strand_seq``.

    Returns a list of ``(start, stop, protein)`` in *strand-local* 0-based
    half-open coordinates. ``start`` is the ATG position, ``stop`` is one past
    the stop codon.
    """

    results: List[tuple[int, int, str]] = []
    n = len(strand_seq)
    start: int | None = None
    i = offset
    while i + 3 <= n:
        codon = strand_seq[i : i + 3]
        if start is None:
            if codon == START_CODON:
                start = i
        else:
            if codon in STOP_CODONS:
                stop = i + 3
                protein = translate(strand_seq[start:stop], frame=1, to_stop=True)
                results.append((start, stop, protein))
                start = None
        i += 3
    return results


def find_orfs(
    sequence: str,
    *,
    min_length: int = 0,
    longest_only: bool = False,
) -> List[ORF]:
    """Find ORFs in all six reading frames.

    Parameters
    ----------
    min_length:
        Minimum ORF length in base pairs (inclusive). ORFs shorter than this
        are discarded.
    longest_only:
        When ``True`` only the single longest ORF (across all frames) is
        returned; ties break by earliest forward-strand start.
    """

    seq = normalize(sequence)
    n = len(seq)
    rc = reverse_complement(seq)
    orfs: List[ORF] = []

    for offset in range(3):
        # Forward strand: frame = offset + 1.
        for start, stop, protein in _scan_frame(seq, offset):
            orf = ORF(
                frame=offset + 1,
                start=start,
                stop=stop,
                strand="+",
                protein=protein,
            )
            if orf.length >= min_length:
                orfs.append(orf)

        # Reverse strand: frame = -(offset + 1). Map rc coords back to forward.
        for rc_start, rc_stop, protein in _scan_frame(rc, offset):
            fwd_start = n - rc_stop
            fwd_stop = n - rc_start
            orf = ORF(
                frame=-(offset + 1),
                start=fwd_start,
                stop=fwd_stop,
                strand="-",
                protein=protein,
            )
            if orf.length >= min_length:
                orfs.append(orf)

    orfs.sort(key=lambda o: (o.start, o.frame))

    if longest_only:
        if not orfs:
            return []
        best = max(orfs, key=lambda o: (o.length, -o.start))
        return [best]

    return orfs
