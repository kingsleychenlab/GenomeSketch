"""Exact motif search implemented from scratch.

Features:

* exact matching
* overlapping matches (e.g. ``AAA`` in ``AAAAA`` -> positions 0, 1, 2)
* optional reverse-complement matching
* case-insensitive search

All positions are 0-based indices into the normalised (upper-case) sequence.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List

from genomesketch.transform import reverse_complement
from genomesketch.validation import normalize


@dataclass
class MotifMatch:
    """A single motif hit."""

    start: int  # 0-based start position
    end: int  # 0-based half-open end position
    strand: str  # "+" for forward, "-" for reverse-complement
    motif: str  # the motif as searched on this strand

    @property
    def length(self) -> int:
        return self.end - self.start


def find_motif(sequence: str, motif: str) -> List[int]:
    """Return all 0-based start positions of ``motif`` (overlapping).

    Example
    -------
    >>> find_motif("AAAAA", "AAA")
    [0, 1, 2]
    """

    seq = normalize(sequence)
    pat = normalize(motif)
    if not pat:
        raise ValueError("motif must be a non-empty string")

    positions: List[int] = []
    start = 0
    limit = len(seq) - len(pat)
    while start <= limit:
        idx = seq.find(pat, start)
        if idx == -1:
            break
        positions.append(idx)
        start = idx + 1  # +1 gives overlapping matches
    return positions


def search_motif(
    sequence: str,
    motif: str,
    *,
    revcomp: bool = False,
) -> List[MotifMatch]:
    """Search for ``motif`` returning structured :class:`MotifMatch` objects.

    When ``revcomp`` is ``True`` the reverse complement of the motif is also
    searched against the forward strand; those hits are reported with a ``-``
    strand and sorted alongside forward hits by start position.
    """

    pat = normalize(motif)
    matches: List[MotifMatch] = []
    for start in find_motif(sequence, pat):
        matches.append(
            MotifMatch(start=start, end=start + len(pat), strand="+", motif=pat)
        )

    if revcomp:
        rc = reverse_complement(pat)
        for start in find_motif(sequence, rc):
            matches.append(
                MotifMatch(
                    start=start, end=start + len(rc), strand="-", motif=rc
                )
            )

    matches.sort(key=lambda m: (m.start, m.strand))
    return matches


def count_motif(sequence: str, motif: str, *, revcomp: bool = False) -> int:
    """Count motif occurrences (overlapping, optionally both strands)."""
    return len(search_motif(sequence, motif, revcomp=revcomp))
