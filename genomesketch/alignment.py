"""Pairwise sequence alignment via dynamic programming (from scratch).

Implements Needleman-Wunsch global alignment and, as a stretch feature,
Smith-Waterman local alignment. Default scoring:

* match:    +1
* mismatch: -1
* gap:      -2

The alignment builds an explicit scoring matrix and a traceback matrix; no
external alignment library is used.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Tuple

from genomesketch.validation import normalize

# Traceback direction markers.
_DIAG = 1  # match / mismatch (came from [i-1][j-1])
_UP = 2  # gap in sequence 2 (came from [i-1][j])
_LEFT = 3  # gap in sequence 1 (came from [i][j-1])
_STOP = 0  # boundary / local-alignment terminator

GAP_CHAR = "-"


@dataclass
class AlignmentResult:
    """Result of a pairwise alignment."""

    aligned1: str
    aligned2: str
    score: int
    match: str  # symbol per column, "match", "gap", "identity" info encoded
    match_type: str  # "global" or "local"
    identity: float  # percent identity over aligned columns (0-100)
    matches: int
    mismatches: int
    gaps: int

    def format_pretty(
        self,
        name1: str = "seq1",
        name2: str = "seq2",
        width: int = 60,
    ) -> str:
        """Return an EMBOSS-Needle-style block for terminal display."""
        lines: List[str] = []
        label_w = max(len(name1), len(name2), 4)
        for i in range(0, len(self.aligned1), width):
            a = self.aligned1[i : i + width]
            b = self.aligned2[i : i + width]
            m = self.match[i : i + width]
            lines.append(f"{name1:<{label_w}}  {a}")
            lines.append(f"{'':<{label_w}}  {m}")
            lines.append(f"{name2:<{label_w}}  {b}")
            lines.append("")
        return "\n".join(lines).rstrip("\n")


def _match_line(a: str, b: str) -> Tuple[str, int, int, int]:
    """Build the middle annotation line and count matches/mismatches/gaps."""
    symbols: List[str] = []
    matches = mismatches = gaps = 0
    for x, y in zip(a, b):
        if x == GAP_CHAR or y == GAP_CHAR:
            symbols.append(" ")
            gaps += 1
        elif x == y:
            symbols.append("|")
            matches += 1
        else:
            symbols.append(".")
            mismatches += 1
    return "".join(symbols), matches, mismatches, gaps


def needleman_wunsch(
    seq1: str,
    seq2: str,
    *,
    match: int = 1,
    mismatch: int = -1,
    gap: int = -2,
) -> AlignmentResult:
    """Global alignment (Needleman-Wunsch).

    Example
    -------
    >>> res = needleman_wunsch("ATGCGT", "ATGACGT")
    >>> res.score
    4
    """

    a = normalize(seq1)
    b = normalize(seq2)
    n, m = len(a), len(b)

    # Scoring matrix (n+1) x (m+1) and traceback matrix.
    score = [[0] * (m + 1) for _ in range(n + 1)]
    trace = [[_STOP] * (m + 1) for _ in range(n + 1)]

    for i in range(1, n + 1):
        score[i][0] = i * gap
        trace[i][0] = _UP
    for j in range(1, m + 1):
        score[0][j] = j * gap
        trace[0][j] = _LEFT

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = score[i - 1][j - 1] + (match if a[i - 1] == b[j - 1] else mismatch)
            up = score[i - 1][j] + gap
            left = score[i][j - 1] + gap
            best = max(diag, up, left)
            score[i][j] = best
            if best == diag:
                trace[i][j] = _DIAG
            elif best == up:
                trace[i][j] = _UP
            else:
                trace[i][j] = _LEFT

    aligned1, aligned2 = _traceback(a, b, trace, i=n, j=m, stop_at_zero=False, score=score)
    final_score = score[n][m]
    return _finish(aligned1, aligned2, final_score, "global")


def smith_waterman(
    seq1: str,
    seq2: str,
    *,
    match: int = 1,
    mismatch: int = -1,
    gap: int = -2,
) -> AlignmentResult:
    """Local alignment (Smith-Waterman). Stretch feature."""

    a = normalize(seq1)
    b = normalize(seq2)
    n, m = len(a), len(b)

    score = [[0] * (m + 1) for _ in range(n + 1)]
    trace = [[_STOP] * (m + 1) for _ in range(n + 1)]

    best_score = 0
    best_pos = (0, 0)
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = score[i - 1][j - 1] + (match if a[i - 1] == b[j - 1] else mismatch)
            up = score[i - 1][j] + gap
            left = score[i][j - 1] + gap
            best = max(0, diag, up, left)
            score[i][j] = best
            if best == 0:
                trace[i][j] = _STOP
            elif best == diag:
                trace[i][j] = _DIAG
            elif best == up:
                trace[i][j] = _UP
            else:
                trace[i][j] = _LEFT
            if best > best_score:
                best_score = best
                best_pos = (i, j)

    aligned1, aligned2 = _traceback(
        a, b, trace, i=best_pos[0], j=best_pos[1], stop_at_zero=True, score=score
    )
    return _finish(aligned1, aligned2, best_score, "local")


def _traceback(
    a: str,
    b: str,
    trace: List[List[int]],
    *,
    i: int,
    j: int,
    stop_at_zero: bool,
    score: List[List[int]],
) -> Tuple[str, str]:
    """Walk the traceback matrix from ``(i, j)`` back to the origin/zero."""
    aligned1: List[str] = []
    aligned2: List[str] = []
    while i > 0 or j > 0:
        direction = trace[i][j]
        if stop_at_zero and (direction == _STOP or score[i][j] == 0):
            break
        if direction == _DIAG:
            aligned1.append(a[i - 1])
            aligned2.append(b[j - 1])
            i -= 1
            j -= 1
        elif direction == _UP:
            aligned1.append(a[i - 1])
            aligned2.append(GAP_CHAR)
            i -= 1
        elif direction == _LEFT:
            aligned1.append(GAP_CHAR)
            aligned2.append(b[j - 1])
            j -= 1
        else:  # _STOP at a boundary in global mode
            if i > 0:
                aligned1.append(a[i - 1])
                aligned2.append(GAP_CHAR)
                i -= 1
            elif j > 0:
                aligned1.append(GAP_CHAR)
                aligned2.append(b[j - 1])
                j -= 1
            else:
                break
    aligned1.reverse()
    aligned2.reverse()
    return "".join(aligned1), "".join(aligned2)


def _finish(
    aligned1: str, aligned2: str, score: int, kind: str
) -> AlignmentResult:
    match_line, matches, mismatches, gaps = _match_line(aligned1, aligned2)
    columns = len(aligned1)
    identity = (matches / columns * 100.0) if columns else 0.0
    return AlignmentResult(
        aligned1=aligned1,
        aligned2=aligned2,
        score=score,
        match=match_line,
        match_type=kind,
        identity=identity,
        matches=matches,
        mismatches=mismatches,
        gaps=gaps,
    )
