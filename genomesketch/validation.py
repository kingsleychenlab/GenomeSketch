"""DNA sequence validation and normalisation.

The MVP alphabet is ``A, T, C, G, N`` (case-insensitive). ``N`` denotes an
unknown base. The validator is structured so that the full IUPAC ambiguity
alphabet can be enabled by passing ``iupac=True``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Set

# MVP alphabet.
DNA_ALPHABET: Set[str] = {"A", "T", "C", "G", "N"}

# Full IUPAC nucleotide ambiguity codes (https://www.bioinformatics.org/sms/iupac.html).
IUPAC_ALPHABET: Set[str] = {
    "A", "C", "G", "T", "U", "N",
    "R", "Y", "S", "W", "K", "M",
    "B", "D", "H", "V",
    "-",
}


@dataclass
class InvalidBase:
    """One invalid character found during validation."""

    position: int  # 0-based index into the (normalised) sequence
    base: str

    def __str__(self) -> str:  # pragma: no cover - trivial
        return f"{self.base!r} at position {self.position}"


class ValidationError(ValueError):
    """Raised when a sequence contains characters outside the alphabet."""

    def __init__(self, invalid: List[InvalidBase], sequence_id: str | None = None):
        self.invalid = invalid
        self.sequence_id = sequence_id
        preview = ", ".join(str(b) for b in invalid[:5])
        if len(invalid) > 5:
            preview += f", ... ({len(invalid)} total)"
        loc = f" in sequence {sequence_id!r}" if sequence_id else ""
        super().__init__(f"Invalid nucleotide(s){loc}: {preview}")


def normalize(sequence: str) -> str:
    """Return an upper-cased copy of ``sequence`` with whitespace removed."""
    return "".join(sequence.split()).upper()


def find_invalid(sequence: str, *, iupac: bool = False) -> List[InvalidBase]:
    """Return a list of invalid bases (with 0-based positions).

    The sequence is normalised (upper-cased, whitespace stripped) before the
    check, so positions refer to the normalised string.
    """

    alphabet = IUPAC_ALPHABET if iupac else DNA_ALPHABET
    normalised = normalize(sequence)
    invalid: List[InvalidBase] = []
    for index, base in enumerate(normalised):
        if base not in alphabet:
            invalid.append(InvalidBase(position=index, base=base))
    return invalid


def is_valid(sequence: str, *, iupac: bool = False) -> bool:
    """Return ``True`` when every character is in the chosen alphabet."""
    return not find_invalid(sequence, iupac=iupac)


def validate(
    sequence: str,
    *,
    iupac: bool = False,
    sequence_id: str | None = None,
) -> str:
    """Validate and normalise ``sequence``.

    Returns the normalised (upper-case) sequence when valid, otherwise raises
    :class:`ValidationError` describing every offending character.
    """

    invalid = find_invalid(sequence, iupac=iupac)
    if invalid:
        raise ValidationError(invalid, sequence_id=sequence_id)
    return normalize(sequence)
