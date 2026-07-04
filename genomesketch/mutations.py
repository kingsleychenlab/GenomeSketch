"""Mutation comparison between two sequences.

Two modes are provided:

* ``compare_equal_length`` -- fast substitution-only comparison for sequences
  of the same length.
* ``compare_aligned`` -- aligns the two sequences with Needleman-Wunsch first
  and then reports substitutions, insertions and deletions.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List

from genomesketch.alignment import GAP_CHAR, needleman_wunsch
from genomesketch.validation import normalize


@dataclass
class Mutation:
    """A single difference between two sequences.

    ``kind`` is one of ``"substitution"``, ``"insertion"`` or ``"deletion"``.
    ``position`` is a 0-based coordinate on the reference sequence.
    """

    position: int
    ref: str
    alt: str
    kind: str

    def __str__(self) -> str:
        if self.kind == "substitution":
            return f"Substitution at position {self.position}: {self.ref} -> {self.alt}"
        if self.kind == "insertion":
            return f"Insertion at position {self.position}: +{self.alt}"
        return f"Deletion at position {self.position}: -{self.ref}"


@dataclass
class ComparisonResult:
    """Outcome of comparing two sequences."""

    mutations: List[Mutation] = field(default_factory=list)
    identity: float = 0.0  # percent identity (0-100)
    aligned_ref: str = ""
    aligned_alt: str = ""
    length_ref: int = 0
    length_alt: int = 0
    mode: str = "equal-length"

    @property
    def substitutions(self) -> List[Mutation]:
        return [m for m in self.mutations if m.kind == "substitution"]

    @property
    def insertions(self) -> List[Mutation]:
        return [m for m in self.mutations if m.kind == "insertion"]

    @property
    def deletions(self) -> List[Mutation]:
        return [m for m in self.mutations if m.kind == "deletion"]


def compare_equal_length(ref: str, alt: str) -> ComparisonResult:
    """Compare two equal-length sequences by direct substitution scan."""
    a = normalize(ref)
    b = normalize(alt)
    if len(a) != len(b):
        raise ValueError(
            f"Sequences differ in length ({len(a)} vs {len(b)}). "
            "Use compare_aligned() for indel-aware comparison."
        )

    mutations: List[Mutation] = []
    matches = 0
    for i, (x, y) in enumerate(zip(a, b)):
        if x == y:
            matches += 1
        else:
            mutations.append(Mutation(position=i, ref=x, alt=y, kind="substitution"))

    identity = (matches / len(a) * 100.0) if a else 0.0
    return ComparisonResult(
        mutations=mutations,
        identity=identity,
        aligned_ref=a,
        aligned_alt=b,
        length_ref=len(a),
        length_alt=len(b),
        mode="equal-length",
    )


def compare_aligned(
    ref: str,
    alt: str,
    *,
    match: int = 1,
    mismatch: int = -1,
    gap: int = -2,
) -> ComparisonResult:
    """Align with Needleman-Wunsch then report substitutions and indels.

    Positions are reported against the reference coordinate system (0-based),
    counting only reference bases consumed so far.
    """

    result = needleman_wunsch(ref, alt, match=match, mismatch=mismatch, gap=gap)
    a, b = result.aligned1, result.aligned2
    mutations: List[Mutation] = []
    ref_pos = 0
    matches = 0
    for x, y in zip(a, b):
        if x == GAP_CHAR:
            # Insertion relative to reference (base present only in alt).
            mutations.append(
                Mutation(position=ref_pos, ref="-", alt=y, kind="insertion")
            )
        elif y == GAP_CHAR:
            mutations.append(
                Mutation(position=ref_pos, ref=x, alt="-", kind="deletion")
            )
            ref_pos += 1
        elif x == y:
            matches += 1
            ref_pos += 1
        else:
            mutations.append(
                Mutation(position=ref_pos, ref=x, alt=y, kind="substitution")
            )
            ref_pos += 1

    columns = len(a)
    identity = (matches / columns * 100.0) if columns else 0.0
    return ComparisonResult(
        mutations=mutations,
        identity=identity,
        aligned_ref=a,
        aligned_alt=b,
        length_ref=len(normalize(ref)),
        length_alt=len(normalize(alt)),
        mode="aligned",
    )


def compare(ref: str, alt: str, *, force_align: bool = False) -> ComparisonResult:
    """Compare two sequences, choosing the appropriate strategy.

    Equal-length sequences use the fast substitution scan unless
    ``force_align`` is set; otherwise Needleman-Wunsch alignment is used.
    """
    a = normalize(ref)
    b = normalize(alt)
    if len(a) == len(b) and not force_align:
        return compare_equal_length(a, b)
    return compare_aligned(a, b)
