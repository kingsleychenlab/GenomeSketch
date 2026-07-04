"""Command line interface for GenomeSketch (argparse-based, no dependencies).

Run ``python -m genomesketch.cli --help`` for the full command list.
"""

from __future__ import annotations

import argparse
import sys
from typing import List, Sequence

from genomesketch import __version__
from genomesketch.io import FastaParseError, SequenceRecord, read_sequences
from genomesketch.validation import ValidationError, validate


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def _load(path: str, *, iupac: bool = False) -> List[SequenceRecord]:
    """Load and validate every sequence in ``path``.

    Exits with a helpful message on I/O or validation errors.
    """
    try:
        records = read_sequences(path)
    except FileNotFoundError:
        _die(f"File not found: {path}")
    except FastaParseError as exc:
        _die(str(exc))

    for record in records:
        try:
            record.sequence = validate(
                record.sequence, iupac=iupac, sequence_id=record.id
            )
        except ValidationError as exc:
            _die(
                f"{exc}\n"
                "Tip: MVP supports A, T, C, G, N. Re-run with --iupac to allow "
                "full IUPAC ambiguity codes."
            )
    return records


def _load_one(path: str, *, iupac: bool = False) -> SequenceRecord:
    """Load a single sequence; if the file has several, use the first."""
    records = _load(path, iupac=iupac)
    if len(records) > 1:
        print(
            f"[note] {path} contains {len(records)} sequences; using the first "
            f"({records[0].id}).",
            file=sys.stderr,
        )
    return records[0]


def _die(message: str) -> "None":
    print(f"error: {message}", file=sys.stderr)
    raise SystemExit(1)


def _fmt_pct(value: float) -> str:
    return f"{value * 100:.2f}%"


# --------------------------------------------------------------------------- #
# Command implementations
# --------------------------------------------------------------------------- #
def cmd_stats(args: argparse.Namespace) -> None:
    from genomesketch.stats import compute_stats

    records = _load(args.file, iupac=args.iupac)
    for i, record in enumerate(records):
        if i:
            print()
        s = compute_stats(record.sequence)
        print(f"Sequence: {record.id}")
        print(f"Length: {s.length} bp")
        print(f"A: {s.a}")
        print(f"T: {s.t}")
        print(f"C: {s.c}")
        print(f"G: {s.g}")
        print(f"N: {s.n}")
        print(f"GC content: {_fmt_pct(s.gc_content)}")
        print(f"AT content: {_fmt_pct(s.at_content)}")
        print(f"GC skew: {s.gc_skew:.3f}")
        print(f"AT skew: {s.at_skew:.3f}")


def cmd_motif(args: argparse.Namespace) -> None:
    from genomesketch.motifs import search_motif

    records = _load(args.file, iupac=args.iupac)
    for i, record in enumerate(records):
        if i:
            print()
        matches = search_motif(record.sequence, args.motif, revcomp=args.revcomp)
        print(f"Sequence: {record.id}")
        print(f"Motif: {args.motif.upper()}  (matches: {len(matches)})")
        if args.revcomp:
            print("Strand search: forward + reverse-complement")
        for m in matches:
            print(f"  pos {m.start}-{m.end}  strand {m.strand}  {m.motif}")


def cmd_revcomp(args: argparse.Namespace) -> None:
    from genomesketch.transform import reverse_complement

    records = _load(args.file, iupac=args.iupac)
    for record in records:
        print(f">{record.id} reverse_complement")
        print(reverse_complement(record.sequence))


def cmd_transcribe(args: argparse.Namespace) -> None:
    from genomesketch.transform import transcribe

    records = _load(args.file, iupac=args.iupac)
    for record in records:
        print(f">{record.id} transcribed")
        print(transcribe(record.sequence))


def cmd_translate(args: argparse.Namespace) -> None:
    from genomesketch.translate import translate

    records = _load(args.file, iupac=args.iupac)
    for record in records:
        protein = translate(record.sequence, frame=args.frame)
        print(f">{record.id} frame={args.frame:+d}")
        print(protein)


def cmd_orfs(args: argparse.Namespace) -> None:
    from genomesketch.orfs import find_orfs

    records = _load(args.file, iupac=args.iupac)
    for i, record in enumerate(records):
        if i:
            print()
        orfs = find_orfs(
            record.sequence,
            min_length=args.min_length,
            longest_only=args.longest,
        )
        print(f"Sequence: {record.id}  (ORFs found: {len(orfs)})")
        for idx, orf in enumerate(orfs, start=1):
            print()
            print(orf.describe(idx))


def cmd_compare(args: argparse.Namespace) -> None:
    from genomesketch.mutations import compare

    ref = _load_one(args.reference, iupac=args.iupac)
    alt = _load_one(args.mutated, iupac=args.iupac)
    result = compare(ref.sequence, alt.sequence, force_align=args.align)

    print(f"Reference: {ref.id} ({result.length_ref} bp)")
    print(f"Mutated:   {alt.id} ({result.length_alt} bp)")
    print(f"Mode: {result.mode}")
    print(f"Percent identity: {result.identity:.1f}%")
    print(f"Substitutions: {len(result.substitutions)}", end="")
    if result.mode == "aligned":
        print(
            f"  Insertions: {len(result.insertions)}  "
            f"Deletions: {len(result.deletions)}"
        )
    else:
        print()
    for mutation in result.mutations:
        print(f"  {mutation}")


def cmd_align(args: argparse.Namespace) -> None:
    from genomesketch.alignment import needleman_wunsch, smith_waterman

    ref = _load_one(args.reference, iupac=args.iupac)
    alt = _load_one(args.mutated, iupac=args.iupac)
    if args.local:
        result = smith_waterman(
            ref.sequence, alt.sequence,
            match=args.match, mismatch=args.mismatch, gap=args.gap,
        )
        kind = "Smith-Waterman (local)"
    else:
        result = needleman_wunsch(
            ref.sequence, alt.sequence,
            match=args.match, mismatch=args.mismatch, gap=args.gap,
        )
        kind = "Needleman-Wunsch (global)"

    print(f"Algorithm: {kind}")
    print(f"Score: {result.score}")
    print(f"Identity: {result.identity:.1f}%")
    print(
        f"Matches: {result.matches}  Mismatches: {result.mismatches}  "
        f"Gaps: {result.gaps}"
    )
    print()
    print(result.format_pretty(ref.id, alt.id))


def cmd_repeats(args: argparse.Namespace) -> None:
    from genomesketch.repeats import find_repeats, find_tandem_repeats

    records = _load(args.file, iupac=args.iupac)
    for i, record in enumerate(records):
        if i:
            print()
        print(f"Sequence: {record.id}")
        if args.tandem:
            tandems = find_tandem_repeats(record.sequence, args.k, args.min_count)
            print(f"Tandem repeats (k={args.k}, min-count={args.min_count}): "
                  f"{len(tandems)}")
            for t in tandems:
                start, end = t.region
                print(f"  {t.motif} x{t.copies}  region {start}-{end}")
        else:
            repeats = find_repeats(record.sequence, args.k, args.min_count)
            print(f"Repeats (k={args.k}, min-count={args.min_count}): "
                  f"{len(repeats)}")
            for r in repeats:
                start, end = r.region
                positions = ", ".join(str(p) for p in r.positions)
                print(
                    f"  {r.motif}  copies {r.copies}  region {start}-{end}  "
                    f"positions [{positions}]"
                )


def cmd_plot_gc(args: argparse.Namespace) -> None:
    from genomesketch.visualization import plot_gc_content

    record = _load_one(args.file, iupac=args.iupac)
    out = plot_gc_content(
        record.sequence, args.output, window=args.window, step=args.step
    )
    print(f"Wrote GC-content plot to {out}")


def cmd_plot_bases(args: argparse.Namespace) -> None:
    from genomesketch.visualization import plot_nucleotide_frequency

    record = _load_one(args.file, iupac=args.iupac)
    out = plot_nucleotide_frequency(record.sequence, args.output)
    print(f"Wrote nucleotide-frequency plot to {out}")


def cmd_plot_orfs(args: argparse.Namespace) -> None:
    from genomesketch.orfs import find_orfs
    from genomesketch.visualization import plot_orf_map

    record = _load_one(args.file, iupac=args.iupac)
    orfs = find_orfs(record.sequence, min_length=args.min_length)
    out = plot_orf_map(len(record.sequence), orfs, args.output)
    print(f"Wrote ORF map ({len(orfs)} ORFs) to {out}")


def cmd_plot_motif(args: argparse.Namespace) -> None:
    from genomesketch.motifs import search_motif
    from genomesketch.visualization import plot_motif_map

    record = _load_one(args.file, iupac=args.iupac)
    matches = search_motif(record.sequence, args.motif, revcomp=args.revcomp)
    out = plot_motif_map(len(record.sequence), matches, args.output)
    print(f"Wrote motif map ({len(matches)} hits) to {out}")


def cmd_kmers(args: argparse.Namespace) -> None:
    from genomesketch.repeats import kmer_positions

    records = _load(args.file, iupac=args.iupac)
    for i, record in enumerate(records):
        if i:
            print()
        positions = kmer_positions(record.sequence, args.k)
        counts = sorted(
            ((kmer, len(pos)) for kmer, pos in positions.items()),
            key=lambda kv: (-kv[1], kv[0]),
        )
        print(f"Sequence: {record.id}  (unique {args.k}-mers: {len(counts)})")
        for kmer, count in counts[: args.top]:
            print(f"  {kmer}\t{count}")


# --------------------------------------------------------------------------- #
# Argument parser
# --------------------------------------------------------------------------- #
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="genomesketch",
        description="Offline DNA sequence analysis toolkit.",
    )
    parser.add_argument("--version", action="version", version=f"GenomeSketch {__version__}")
    sub = parser.add_subparsers(dest="command", metavar="<command>")
    sub.required = True

    def add_iupac(p: argparse.ArgumentParser) -> None:
        p.add_argument(
            "--iupac",
            action="store_true",
            help="Allow full IUPAC ambiguity codes (default: A/T/C/G/N only).",
        )

    # stats
    p = sub.add_parser("stats", help="Nucleotide composition statistics.")
    p.add_argument("file")
    add_iupac(p)
    p.set_defaults(func=cmd_stats)

    # motif
    p = sub.add_parser("motif", help="Exact motif search (overlapping).")
    p.add_argument("file")
    p.add_argument("--motif", required=True, help="Motif to search for.")
    p.add_argument("--revcomp", action="store_true", help="Also search the reverse complement.")
    add_iupac(p)
    p.set_defaults(func=cmd_motif)

    # revcomp
    p = sub.add_parser("revcomp", help="Reverse complement.")
    p.add_argument("file")
    add_iupac(p)
    p.set_defaults(func=cmd_revcomp)

    # transcribe
    p = sub.add_parser("transcribe", help="Transcribe DNA to RNA.")
    p.add_argument("file")
    add_iupac(p)
    p.set_defaults(func=cmd_transcribe)

    # translate
    p = sub.add_parser("translate", help="Translate DNA/RNA to protein.")
    p.add_argument("file")
    p.add_argument(
        "--frame", type=int, default=1, choices=[1, 2, 3, -1, -2, -3],
        help="Reading frame (default: 1).",
    )
    add_iupac(p)
    p.set_defaults(func=cmd_translate)

    # orfs
    p = sub.add_parser("orfs", help="Find open reading frames (6 frames).")
    p.add_argument("file")
    p.add_argument("--min-length", type=int, default=0, dest="min_length",
                   help="Minimum ORF length in bp (default: 0).")
    p.add_argument("--longest", action="store_true", help="Report only the longest ORF.")
    add_iupac(p)
    p.set_defaults(func=cmd_orfs)

    # compare
    p = sub.add_parser("compare", help="Compare two sequences for mutations.")
    p.add_argument("reference")
    p.add_argument("mutated")
    p.add_argument("--align", action="store_true",
                   help="Force alignment-based (indel-aware) comparison.")
    add_iupac(p)
    p.set_defaults(func=cmd_compare)

    # align
    p = sub.add_parser("align", help="Pairwise alignment (Needleman-Wunsch).")
    p.add_argument("reference")
    p.add_argument("mutated")
    p.add_argument("--local", action="store_true", help="Smith-Waterman local alignment.")
    p.add_argument("--match", type=int, default=1, help="Match score (default: 1).")
    p.add_argument("--mismatch", type=int, default=-1, help="Mismatch score (default: -1).")
    p.add_argument("--gap", type=int, default=-2, help="Gap penalty (default: -2).")
    add_iupac(p)
    p.set_defaults(func=cmd_align)

    # repeats
    p = sub.add_parser("repeats", help="Detect repeated k-mers / tandem repeats.")
    p.add_argument("file")
    p.add_argument("--k", type=int, default=3, help="k-mer size (default: 3).")
    p.add_argument("--min-count", type=int, default=2, dest="min_count",
                   help="Minimum copy count (default: 2).")
    p.add_argument("--tandem", action="store_true", help="Only tandem (adjacent) repeats.")
    add_iupac(p)
    p.set_defaults(func=cmd_repeats)

    # kmers (stretch)
    p = sub.add_parser("kmers", help="k-mer frequency table.")
    p.add_argument("file")
    p.add_argument("--k", type=int, default=3, help="k-mer size (default: 3).")
    p.add_argument("--top", type=int, default=20, help="Show the top N k-mers (default: 20).")
    add_iupac(p)
    p.set_defaults(func=cmd_kmers)

    # plot-gc
    p = sub.add_parser("plot-gc", help="Sliding-window GC-content plot.")
    p.add_argument("file")
    p.add_argument("--window", type=int, default=100, help="Window size (default: 100).")
    p.add_argument("--step", type=int, default=10, help="Step size (default: 10).")
    p.add_argument("-o", "--output", default="gc_plot.png", help="Output image path.")
    add_iupac(p)
    p.set_defaults(func=cmd_plot_gc)

    # plot-bases
    p = sub.add_parser("plot-bases", help="Nucleotide-frequency bar chart.")
    p.add_argument("file")
    p.add_argument("-o", "--output", default="bases_plot.png", help="Output image path.")
    add_iupac(p)
    p.set_defaults(func=cmd_plot_bases)

    # plot-orfs
    p = sub.add_parser("plot-orfs", help="ORF map figure.")
    p.add_argument("file")
    p.add_argument("--min-length", type=int, default=0, dest="min_length",
                   help="Minimum ORF length in bp (default: 0).")
    p.add_argument("-o", "--output", default="orf_map.png", help="Output image path.")
    add_iupac(p)
    p.set_defaults(func=cmd_plot_orfs)

    # plot-motif
    p = sub.add_parser("plot-motif", help="Motif position map figure.")
    p.add_argument("file")
    p.add_argument("--motif", required=True, help="Motif to map.")
    p.add_argument("--revcomp", action="store_true", help="Include reverse-complement hits.")
    p.add_argument("-o", "--output", default="motif_map.png", help="Output image path.")
    add_iupac(p)
    p.set_defaults(func=cmd_plot_motif)

    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        args.func(args)
    except ValidationError as exc:  # pragma: no cover - defensive
        _die(str(exc))
    except BrokenPipeError:  # pragma: no cover - piping into head, etc.
        return 0
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
