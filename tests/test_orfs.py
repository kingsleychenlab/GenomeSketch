"""Tests for ORF detection."""

from genomesketch.orfs import find_orfs


def test_simple_orf():
    orfs = find_orfs("AAAATGAAATAACCC")
    # Expect an ORF ATG..TAA on the forward strand.
    forward = [o for o in orfs if o.strand == "+"]
    assert forward
    orf = forward[0]
    assert orf.start == 3
    assert orf.stop == 12  # half-open, includes stop codon
    assert orf.length == 9
    assert orf.protein == "MK"
    assert orf.frame == 1


def test_min_length_filter():
    seq = "AAAATGAAATAACCC"
    all_orfs = find_orfs(seq)
    filtered = find_orfs(seq, min_length=100)
    assert all_orfs
    assert filtered == []


def test_longest_only_returns_single():
    seq = "ATGAAATAA" + "GGG" + "ATGAAAAAAAAAAAATAA"
    orfs = find_orfs(seq, longest_only=True)
    assert len(orfs) == 1


def test_reverse_strand_orf():
    from genomesketch.transform import reverse_complement

    # Build a sequence whose reverse complement contains a clean ORF.
    fwd_orf = "ATGAAATAA"
    seq = reverse_complement(fwd_orf)  # so an ORF appears on the - strand
    orfs = find_orfs(seq)
    minus = [o for o in orfs if o.strand == "-"]
    assert minus
    assert minus[0].protein == "MK"


def test_orf_protein_excludes_stop():
    orfs = find_orfs("ATGGAATTTTAA")
    forward = [o for o in orfs if o.strand == "+"]
    assert forward[0].protein == "MEF"
