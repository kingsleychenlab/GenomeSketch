"""Tests for validation and normalisation."""

import pytest

from genomesketch.validation import (
    ValidationError,
    find_invalid,
    is_valid,
    normalize,
    validate,
)


def test_normalize_uppercases_and_strips():
    assert normalize("  atg c\n g ") == "ATGCG"


def test_valid_mvp_alphabet():
    assert is_valid("ATCGNatcgn")


def test_invalid_characters_reported_with_position():
    invalid = find_invalid("ATGXCGZ")
    assert [(b.position, b.base) for b in invalid] == [(3, "X"), (6, "Z")]


def test_validate_returns_normalised():
    assert validate("atgc") == "ATGC"


def test_validate_raises_on_invalid():
    with pytest.raises(ValidationError) as exc:
        validate("ATGQ", sequence_id="seq1")
    assert "seq1" in str(exc.value)
    assert exc.value.invalid[0].base == "Q"


def test_iupac_mode_allows_ambiguity_codes():
    assert not is_valid("ATGR")  # R not in MVP alphabet
    assert is_valid("ATGR", iupac=True)


def test_case_insensitive_validation():
    assert is_valid("atcgn")
    assert validate("atcgn") == "ATCGN"
