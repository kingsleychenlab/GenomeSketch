"""Nucleotide composition statistics.

All formulas follow standard definitions and guard against division by zero:

* GC content = (G + C) / (A + T + G + C)
* AT content = (A + T) / (A + T + G + C)
* GC skew    = (G - C) / (G + C)
* AT skew    = (A - T) / (A + T)
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from typing import Dict

from genomesketch.validation import normalize


@dataclass
class SequenceStats:
    """Composition statistics for a single sequence."""

    length: int
    counts: Dict[str, int]
    gc_content: float
    at_content: float
    gc_skew: float
    at_skew: float
    other: int = field(default=0)

    @property
    def a(self) -> int:
        return self.counts.get("A", 0)

    @property
    def t(self) -> int:
        return self.counts.get("T", 0)

    @property
    def c(self) -> int:
        return self.counts.get("C", 0)

    @property
    def g(self) -> int:
        return self.counts.get("G", 0)

    @property
    def n(self) -> int:
        return self.counts.get("N", 0)

    def as_dict(self) -> Dict[str, float]:
        """Flat dictionary suitable for tabular export."""
        return {
            "length": self.length,
            "A": self.a,
            "T": self.t,
            "C": self.c,
            "G": self.g,
            "N": self.n,
            "other": self.other,
            "gc_content": self.gc_content,
            "at_content": self.at_content,
            "gc_skew": self.gc_skew,
            "at_skew": self.at_skew,
        }


def _safe_div(numerator: float, denominator: float) -> float:
    """Return ``numerator / denominator`` or ``0.0`` when the divisor is 0."""
    if denominator == 0:
        return 0.0
    return numerator / denominator


def nucleotide_counts(sequence: str) -> Dict[str, int]:
    """Count A/T/C/G/N (and anything else) in a normalised sequence."""
    normalised = normalize(sequence)
    counter = Counter(normalised)
    counts = {base: counter.get(base, 0) for base in ("A", "T", "C", "G", "N")}
    return counts


def gc_content(sequence: str) -> float:
    """Fraction of G+C over A+T+G+C (N is excluded from the denominator)."""
    counts = nucleotide_counts(sequence)
    g, c, a, t = counts["G"], counts["C"], counts["A"], counts["T"]
    return _safe_div(g + c, a + t + g + c)


def at_content(sequence: str) -> float:
    """Fraction of A+T over A+T+G+C."""
    counts = nucleotide_counts(sequence)
    g, c, a, t = counts["G"], counts["C"], counts["A"], counts["T"]
    return _safe_div(a + t, a + t + g + c)


def gc_skew(sequence: str) -> float:
    """(G - C) / (G + C)."""
    counts = nucleotide_counts(sequence)
    g, c = counts["G"], counts["C"]
    return _safe_div(g - c, g + c)


def at_skew(sequence: str) -> float:
    """(A - T) / (A + T)."""
    counts = nucleotide_counts(sequence)
    a, t = counts["A"], counts["T"]
    return _safe_div(a - t, a + t)


def compute_stats(sequence: str) -> SequenceStats:
    """Compute the full :class:`SequenceStats` bundle for a sequence."""
    normalised = normalize(sequence)
    counts = nucleotide_counts(normalised)
    known = counts["A"] + counts["T"] + counts["C"] + counts["G"] + counts["N"]
    other = len(normalised) - known
    return SequenceStats(
        length=len(normalised),
        counts=counts,
        gc_content=gc_content(normalised),
        at_content=at_content(normalised),
        gc_skew=gc_skew(normalised),
        at_skew=at_skew(normalised),
        other=other,
    )


def gc_windows(
    sequence: str, window: int = 100, step: int = 10
) -> list[tuple[int, float]]:
    """Sliding-window GC content.

    Returns a list of ``(window_start, gc_fraction)`` tuples using 0-based
    start positions. The final partial window is included if it contains at
    least one base.
    """

    if window <= 0:
        raise ValueError("window must be a positive integer")
    if step <= 0:
        raise ValueError("step must be a positive integer")

    normalised = normalize(sequence)
    results: list[tuple[int, float]] = []
    n = len(normalised)
    if n == 0:
        return results

    start = 0
    while start < n:
        segment = normalised[start : start + window]
        results.append((start, gc_content(segment)))
        if start + window >= n:
            break
        start += step
    return results
