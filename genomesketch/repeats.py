"""Repeat detection: k-mer repeats and tandem repeats.

Two complementary views are offered:

* :func:`find_repeats` groups identical k-mers that occur at least
  ``min_count`` times anywhere in the sequence, reporting their positions and
  the spanned region.
* :func:`find_tandem_repeats` detects *tandem* (immediately adjacent) repeats
  of a k-mer, e.g. ``ATGATGATG``.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List

from genomesketch.validation import normalize


@dataclass
class Repeat:
    """A k-mer that recurs in the sequence."""

    motif: str
    positions: List[int]  # 0-based start positions (sorted)
    copies: int = field(init=False)

    def __post_init__(self) -> None:
        self.copies = len(self.positions)

    @property
    def region(self) -> tuple[int, int]:
        """Half-open region ``(start, end)`` spanned by all copies."""
        start = self.positions[0]
        end = self.positions[-1] + len(self.motif)
        return (start, end)

    def describe(self) -> str:
        start, end = self.region
        return (
            f"Repeat: {self.motif}\n"
            f"Region: {start}-{end}\n"
            f"Copies: {self.copies}"
        )


@dataclass
class TandemRepeat:
    """A k-mer repeated back-to-back."""

    motif: str
    start: int  # 0-based start of the tandem array
    copies: int

    @property
    def region(self) -> tuple[int, int]:
        end = self.start + self.copies * len(self.motif)
        return (self.start, end)


def kmer_positions(sequence: str, k: int) -> Dict[str, List[int]]:
    """Map every length-``k`` substring to its (overlapping) start positions."""
    if k <= 0:
        raise ValueError("k must be a positive integer")
    seq = normalize(sequence)
    positions: Dict[str, List[int]] = defaultdict(list)
    for i in range(len(seq) - k + 1):
        positions[seq[i : i + k]].append(i)
    return dict(positions)


def find_repeats(
    sequence: str,
    k: int,
    min_count: int = 2,
) -> List[Repeat]:
    """Find k-mers occurring at least ``min_count`` times.

    Results are sorted by descending copy count, then by first position.
    """
    if min_count < 1:
        raise ValueError("min_count must be >= 1")
    positions = kmer_positions(sequence, k)
    repeats = [
        Repeat(motif=motif, positions=sorted(pos))
        for motif, pos in positions.items()
        if len(pos) >= min_count
    ]
    repeats.sort(key=lambda r: (-r.copies, r.positions[0]))
    return repeats


def find_tandem_repeats(
    sequence: str,
    k: int,
    min_count: int = 2,
) -> List[TandemRepeat]:
    """Detect tandem (adjacent) repeats of length-``k`` motifs.

    A tandem array is a maximal run of the same k-mer repeated with a stride of
    ``k`` (``ATGATGATG`` -> motif ``ATG``, 3 copies). Only arrays with at least
    ``min_count`` copies are returned.
    """
    if min_count < 1:
        raise ValueError("min_count must be >= 1")
    seq = normalize(sequence)
    n = len(seq)
    tandems: List[TandemRepeat] = []
    i = 0
    while i + k <= n:
        motif = seq[i : i + k]
        copies = 1
        j = i + k
        while j + k <= n and seq[j : j + k] == motif:
            copies += 1
            j += k
        if copies >= min_count:
            tandems.append(TandemRepeat(motif=motif, start=i, copies=copies))
            i = j  # skip past this tandem array
        else:
            i += 1
    return tandems
