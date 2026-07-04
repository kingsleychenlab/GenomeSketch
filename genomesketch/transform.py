"""Sequence transforms: complement, reverse complement, transcription."""

from __future__ import annotations

from genomesketch.validation import normalize

# Complement table for the MVP alphabet plus common IUPAC ambiguity codes so
# that reverse complement still behaves for ambiguous input.
_COMPLEMENT = {
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C",
    "N": "N",
    "U": "A",
    "R": "Y",
    "Y": "R",
    "S": "S",
    "W": "W",
    "K": "M",
    "M": "K",
    "B": "V",
    "V": "B",
    "D": "H",
    "H": "D",
    "-": "-",
}


def complement(sequence: str) -> str:
    """Return the base-by-base complement (5'->3' order preserved).

    Unknown characters are passed through unchanged after normalisation.
    """
    normalised = normalize(sequence)
    return "".join(_COMPLEMENT.get(base, base) for base in normalised)


def reverse_complement(sequence: str) -> str:
    """Return the reverse complement.

    Example
    -------
    >>> reverse_complement("ATGCCG")
    'CGGCAT'
    """
    return complement(sequence)[::-1]


def reverse(sequence: str) -> str:
    """Return the reversed (not complemented) sequence."""
    return normalize(sequence)[::-1]


def transcribe(sequence: str) -> str:
    """Transcribe DNA to RNA (T -> U).

    Example
    -------
    >>> transcribe("ATGCGT")
    'AUGCGU'
    """
    normalised = normalize(sequence)
    return normalised.replace("T", "U")


def back_transcribe(sequence: str) -> str:
    """Reverse transcription: RNA -> DNA (U -> T)."""
    normalised = normalize(sequence)
    return normalised.replace("U", "T")
