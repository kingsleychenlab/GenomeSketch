"""Tests for FASTA / plain-text parsing."""

import pytest

from genomesketch.io import (
    FastaParseError,
    SequenceRecord,
    parse_fasta,
    read_sequences,
    write_fasta,
)


def test_parse_single_record():
    text = ">sample_1\nATGCGTACCGTTAGCTAGCTAGGCTA\n"
    records = parse_fasta(text)
    assert len(records) == 1
    assert records[0] == SequenceRecord(
        id="sample_1", description="", sequence="ATGCGTACCGTTAGCTAGCTAGGCTA"
    )


def test_parse_header_description():
    text = ">seq1 this is a description\nATGC\n"
    (record,) = parse_fasta(text)
    assert record.id == "seq1"
    assert record.description == "this is a description"
    assert record.header == "seq1 this is a description"


def test_parse_multiline_sequence():
    text = ">seq\nATGC\nGGTT\nAACC\n"
    (record,) = parse_fasta(text)
    assert record.sequence == "ATGCGGTTAACC"


def test_parse_blank_lines_and_whitespace():
    text = "\n>seq\n  ATGC  \n\n  GGTT\n\n"
    (record,) = parse_fasta(text)
    assert record.sequence == "ATGCGGTT"


def test_parse_multiple_sequences():
    text = ">a\nATG\n>b\nCCC\n>c desc\nGGG\n"
    records = parse_fasta(text)
    assert [r.id for r in records] == ["a", "b", "c"]
    assert [r.sequence for r in records] == ["ATG", "CCC", "GGG"]
    assert records[2].description == "desc"


def test_parse_plain_text_no_header():
    text = "ATGCGT\nAACCGG\n"
    (record,) = parse_fasta(text)
    assert record.id == "sequence_1"
    assert record.sequence == "ATGCGTAACCGG"


def test_read_sequences_from_file(tmp_path):
    path = tmp_path / "seq.fasta"
    path.write_text(">x\nATGC\n>y\nGGGG\n")
    records = read_sequences(path)
    assert len(records) == 2
    assert records[0].id == "x"


def test_read_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        read_sequences(tmp_path / "nope.fasta")


def test_read_empty_file_raises(tmp_path):
    path = tmp_path / "empty.fasta"
    path.write_text("\n\n")
    with pytest.raises(FastaParseError):
        read_sequences(path)


def test_write_fasta_roundtrip():
    records = [
        SequenceRecord(id="a", description="d", sequence="ATGCATGC"),
        SequenceRecord(id="b", sequence="GGGG"),
    ]
    text = write_fasta(records, width=4)
    reparsed = parse_fasta(text)
    assert reparsed[0].id == "a"
    assert reparsed[0].sequence == "ATGCATGC"
    assert reparsed[1].sequence == "GGGG"


def test_example_multi_sequence_file():
    records = read_sequences("examples/multi_sequence.fasta")
    assert [r.id for r in records] == ["seq_A", "seq_B", "seq_C"]
