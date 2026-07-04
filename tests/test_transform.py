"""Tests for reverse complement and transcription."""

from genomesketch.transform import (
    back_transcribe,
    complement,
    reverse_complement,
    transcribe,
)


def test_reverse_complement():
    assert reverse_complement("ATGCCG") == "CGGCAT"


def test_reverse_complement_n():
    assert reverse_complement("ATGN") == "NCAT"


def test_complement():
    assert complement("ATGC") == "TACG"


def test_reverse_complement_case_insensitive():
    assert reverse_complement("atgccg") == "CGGCAT"


def test_transcribe():
    assert transcribe("ATGCGT") == "AUGCGU"


def test_back_transcribe():
    assert back_transcribe("AUGCGU") == "ATGCGT"


def test_transcribe_roundtrip():
    seq = "ATGCGTACCGTTAG"
    assert back_transcribe(transcribe(seq)) == seq
