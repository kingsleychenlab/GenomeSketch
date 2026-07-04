"""Tests for translation and reading frames."""

from genomesketch.translate import (
    six_frame_translation,
    translate,
    translate_codon,
)


def test_translate_basic():
    assert translate("ATGGAATTTTAA") == "MEF*"


def test_translate_ignores_incomplete_codon():
    # 13 bases -> 4 full codons + 1 leftover base ignored.
    assert translate("ATGGAATTTTAAA") == "MEF*"


def test_translate_rna_input():
    assert translate("AUGGAAUUUUAA") == "MEF*"


def test_translate_n_gives_x():
    assert translate("ATGNNNTTT") == "MXF"


def test_translate_frame_2():
    # Shift by one base.
    seq = "AATGGAATTTTAA"
    assert translate(seq, frame=2) == "MEF*"


def test_translate_frame_3():
    seq = "AAATGGAATTTTAA"
    assert translate(seq, frame=3) == "MEF*"


def test_translate_reverse_frame():
    from genomesketch.transform import reverse_complement

    seq = "ATGGAATTTTAA"
    rc = reverse_complement(seq)
    # Forward frame 1 of rc equals reverse frame -1 of seq.
    assert translate(seq, frame=-1) == translate(rc, frame=1)


def test_to_stop_truncates():
    assert translate("ATGGAATAAGGG", to_stop=True) == "ME"


def test_six_frames_keys():
    frames = six_frame_translation("ATGGAATTTTAA")
    assert set(frames) == {1, 2, 3, -1, -2, -3}


def test_translate_codon():
    assert translate_codon("ATG") == "M"
    assert translate_codon("TAA") == "*"
    assert translate_codon("NNN") == "X"
