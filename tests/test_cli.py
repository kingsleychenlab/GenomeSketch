"""CLI smoke tests exercising the argparse entry point."""

import pytest

from genomesketch.cli import main

SAMPLE = "examples/sample.fasta"
REF = "examples/reference.fasta"
MUT = "examples/mutated.fasta"


def run(argv, capsys):
    code = main(argv)
    out = capsys.readouterr()
    return code, out.out, out.err


def test_stats(capsys):
    code, out, _ = run(["stats", SAMPLE], capsys)
    assert code == 0
    assert "Length:" in out
    assert "GC content:" in out


def test_motif(capsys):
    code, out, _ = run(["motif", SAMPLE, "--motif", "ATG"], capsys)
    assert code == 0
    assert "Motif: ATG" in out


def test_motif_revcomp(capsys):
    code, out, _ = run(["motif", SAMPLE, "--motif", "ATG", "--revcomp"], capsys)
    assert code == 0
    assert "reverse-complement" in out


def test_revcomp(capsys):
    code, out, _ = run(["revcomp", SAMPLE], capsys)
    assert code == 0
    assert out.startswith(">")


def test_transcribe(capsys):
    code, out, _ = run(["transcribe", SAMPLE], capsys)
    assert code == 0
    assert "U" in out


def test_translate(capsys):
    code, out, _ = run(["translate", SAMPLE, "--frame", "1"], capsys)
    assert code == 0


def test_orfs(capsys):
    code, out, _ = run(["orfs", SAMPLE], capsys)
    assert code == 0
    assert "ORFs found" in out


def test_compare(capsys):
    code, out, _ = run(["compare", REF, MUT], capsys)
    assert code == 0
    assert "Percent identity:" in out


def test_align(capsys):
    code, out, _ = run(["align", REF, MUT], capsys)
    assert code == 0
    assert "Score:" in out


def test_repeats(capsys):
    code, out, _ = run(["repeats", SAMPLE, "--k", "3", "--min-count", "2"], capsys)
    assert code == 0
    assert "Repeats" in out


def test_plot_gc(tmp_path, capsys):
    out_path = str(tmp_path / "gc.png")
    code, out, _ = run(
        ["plot-gc", SAMPLE, "--window", "10", "--step", "2", "-o", out_path],
        capsys,
    )
    assert code == 0
    import os

    assert os.path.getsize(out_path) > 0


def test_plot_bases(tmp_path, capsys):
    out_path = str(tmp_path / "bases.png")
    code, _, _ = run(["plot-bases", SAMPLE, "-o", out_path], capsys)
    assert code == 0
    import os

    assert os.path.getsize(out_path) > 0


def test_missing_file_errors(capsys):
    with pytest.raises(SystemExit) as exc:
        main(["stats", "does_not_exist.fasta"])
    assert exc.value.code == 1


def test_invalid_nucleotide_errors(tmp_path, capsys):
    bad = tmp_path / "bad.fasta"
    bad.write_text(">bad\nATGQ\n")
    with pytest.raises(SystemExit):
        main(["stats", str(bad)])
