"""Matplotlib visualisations saved to local image files.

Every function renders a real plot from the supplied data and writes it to
``output`` (PNG or any Matplotlib-supported format). A non-interactive backend
is forced so plots can be generated headlessly and offline.
"""

from __future__ import annotations

from typing import List, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")  # headless, offline-safe backend

import matplotlib.pyplot as plt  # noqa: E402  (must come after use())

from genomesketch.motifs import MotifMatch
from genomesketch.orfs import ORF
from genomesketch.stats import gc_windows, nucleotide_counts


def plot_gc_content(
    sequence: str,
    output: str,
    *,
    window: int = 100,
    step: int = 10,
    title: str = "GC content (sliding window)",
) -> str:
    """Sliding-window GC-content line plot. Returns ``output``."""
    windows = gc_windows(sequence, window=window, step=step)
    xs = [start for start, _ in windows]
    ys = [gc * 100.0 for _, gc in windows]

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(xs, ys, color="#2a6f97", linewidth=1.5)
    ax.axhline(50, color="grey", linestyle="--", linewidth=0.8, alpha=0.7)
    ax.set_xlabel("Position (bp)")
    ax.set_ylabel("GC content (%)")
    ax.set_title(title)
    ax.set_ylim(0, 100)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(output, dpi=120)
    plt.close(fig)
    return output


def plot_nucleotide_frequency(
    sequence: str,
    output: str,
    *,
    title: str = "Nucleotide frequency",
) -> str:
    """Bar chart of A/T/C/G/N counts. Returns ``output``."""
    counts = nucleotide_counts(sequence)
    bases = ["A", "T", "C", "G", "N"]
    values = [counts[b] for b in bases]
    colors = ["#e07a5f", "#3d405b", "#81b29a", "#f2cc8f", "#999999"]

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(bases, values, color=colors)
    for i, v in enumerate(values):
        ax.text(i, v, str(v), ha="center", va="bottom")
    ax.set_xlabel("Nucleotide")
    ax.set_ylabel("Count")
    ax.set_title(title)
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(output, dpi=120)
    plt.close(fig)
    return output


def plot_motif_map(
    sequence_length: int,
    matches: Sequence[MotifMatch],
    output: str,
    *,
    title: str = "Motif positions",
) -> str:
    """Position map of motif hits along the sequence. Returns ``output``."""
    fig, ax = plt.subplots(figsize=(10, 3))
    ax.hlines(0, 0, sequence_length, color="black", linewidth=1)
    for match in matches:
        color = "#2a9d8f" if match.strand == "+" else "#e76f51"
        y = 0.15 if match.strand == "+" else -0.15
        ax.plot(
            [match.start, match.end],
            [y, y],
            color=color,
            linewidth=6,
            solid_capstyle="butt",
        )
    ax.set_ylim(-1, 1)
    ax.set_xlim(0, max(sequence_length, 1))
    ax.set_yticks([0.15, -0.15])
    ax.set_yticklabels(["+ strand", "- strand"])
    ax.set_xlabel("Position (bp)")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(output, dpi=120)
    plt.close(fig)
    return output


def plot_orf_map(
    sequence_length: int,
    orfs: Sequence[ORF],
    output: str,
    *,
    title: str = "ORF map",
) -> str:
    """Draw ORFs on a per-frame track. Returns ``output``."""
    fig, ax = plt.subplots(figsize=(10, 5))
    frame_order = [3, 2, 1, -1, -2, -3]
    y_for = {frame: idx for idx, frame in enumerate(frame_order)}

    ax.hlines(
        list(y_for.values()),
        0,
        sequence_length,
        color="#dddddd",
        linewidth=1,
    )
    for orf in orfs:
        y = y_for.get(orf.frame, 0)
        color = "#264653" if orf.strand == "+" else "#e76f51"
        ax.plot(
            [orf.start, orf.stop],
            [y, y],
            color=color,
            linewidth=8,
            solid_capstyle="butt",
        )
    ax.set_yticks(list(y_for.values()))
    ax.set_yticklabels([f"frame {f:+d}" for f in frame_order])
    ax.set_xlim(0, max(sequence_length, 1))
    ax.set_xlabel("Position (bp)")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(output, dpi=120)
    plt.close(fig)
    return output


def plot_alignment_mismatches(
    match_line: str,
    output: str,
    *,
    title: str = "Alignment mismatch map",
) -> str:
    """Bar strip marking matches / mismatches / gaps across the alignment."""
    xs = list(range(len(match_line)))
    categories: List[Tuple[str, str]] = []
    for symbol in match_line:
        if symbol == "|":
            categories.append(("match", "#2a9d8f"))
        elif symbol == ".":
            categories.append(("mismatch", "#e63946"))
        else:
            categories.append(("gap", "#adb5bd"))

    fig, ax = plt.subplots(figsize=(10, 2.5))
    for x, (_, color) in zip(xs, categories):
        ax.vlines(x, 0, 1, color=color, linewidth=1.5)
    ax.set_xlim(0, max(len(match_line), 1))
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_xlabel("Alignment column")
    ax.set_title(title)

    # Manual legend.
    from matplotlib.patches import Patch

    legend = [
        Patch(color="#2a9d8f", label="match"),
        Patch(color="#e63946", label="mismatch"),
        Patch(color="#adb5bd", label="gap"),
    ]
    ax.legend(handles=legend, loc="upper right", ncol=3, fontsize=8)
    fig.tight_layout()
    fig.savefig(output, dpi=120)
    plt.close(fig)
    return output
