"""FASTA / plain-text DNA parsing implemented from scratch.

The parser intentionally avoids Biopython so that the behaviour is fully
transparent and offline. It supports:

* ``.fa``, ``.fasta``, ``.fna`` and plain-text DNA files
* multiple sequences per file
* header lines that begin with ``>``
* sequences that span multiple lines
* blank lines and surrounding whitespace (ignored)

A plain-text file with no ``>`` header is treated as a single anonymous
sequence.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Iterable, Iterator, List


@dataclass
class SequenceRecord:
    """A single parsed sequence.

    Attributes
    ----------
    id:
        The identifier taken from the FASTA header (the first whitespace
        delimited token after ``>``).
    description:
        The remainder of the header line after the id. Empty string when the
        header carries no description.
    sequence:
        The concatenated sequence with all internal whitespace removed. The raw
        characters are preserved as-is (case is *not* changed here; use
        :mod:`genomesketch.validation` to normalise).
    """

    id: str
    description: str = ""
    sequence: str = ""
    _line: int = field(default=0, repr=False, compare=False)

    def __len__(self) -> int:
        return len(self.sequence)

    @property
    def header(self) -> str:
        """Reconstruct the ``id description`` header (without the ``>``)."""
        if self.description:
            return f"{self.id} {self.description}"
        return self.id


class FastaParseError(ValueError):
    """Raised when a FASTA/DNA file cannot be parsed."""


def _split_header(header: str) -> tuple[str, str]:
    """Split a header line (already stripped of ``>``) into id + description."""
    header = header.strip()
    if not header:
        return "", ""
    parts = header.split(None, 1)
    seq_id = parts[0]
    description = parts[1].strip() if len(parts) > 1 else ""
    return seq_id, description


def parse_fasta(text: str) -> List[SequenceRecord]:
    """Parse FASTA/plain-text ``text`` into a list of :class:`SequenceRecord`.

    Blank lines and whitespace are ignored. A body without any ``>`` header is
    returned as a single record with an auto-generated id of ``sequence_1``.
    """

    records: List[SequenceRecord] = []
    current_id: str | None = None
    current_desc = ""
    current_line = 0
    chunks: List[str] = []
    anonymous_chunks: List[str] = []
    anonymous_line = 0

    def flush() -> None:
        nonlocal current_id, current_desc, chunks
        if current_id is not None:
            records.append(
                SequenceRecord(
                    id=current_id,
                    description=current_desc,
                    sequence="".join(chunks),
                    _line=current_line,
                )
            )
        current_id = None
        current_desc = ""
        chunks = []

    for lineno, raw_line in enumerate(text.splitlines(), start=1):
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            flush()
            current_id, current_desc = _split_header(line[1:])
            current_line = lineno
            if not current_id:
                current_id = f"sequence_{len(records) + 1}"
        elif line.startswith(";"):
            # Legacy FASTA comment lines start with ';' -- ignore them.
            continue
        else:
            # Remove any internal whitespace (rare but permitted).
            cleaned = "".join(line.split())
            if current_id is None:
                if not anonymous_chunks:
                    anonymous_line = lineno
                anonymous_chunks.append(cleaned)
            else:
                chunks.append(cleaned)

    flush()

    if anonymous_chunks:
        # A plain-text DNA body with no header. Insert it as the first record.
        records.insert(
            0,
            SequenceRecord(
                id="sequence_1",
                description="",
                sequence="".join(anonymous_chunks),
                _line=anonymous_line,
            ),
        )
        # Re-number any auto-generated ids so they stay unique/monotonic.
        _renumber_anonymous(records)

    return records


def _renumber_anonymous(records: List[SequenceRecord]) -> None:
    counter = 1
    for record in records:
        if record.id.startswith("sequence_"):
            record.id = f"sequence_{counter}"
            counter += 1


def read_sequences(path: str | os.PathLike[str]) -> List[SequenceRecord]:
    """Read and parse a FASTA/plain-text DNA file from disk.

    Raises
    ------
    FileNotFoundError:
        If ``path`` does not exist.
    FastaParseError:
        If the file is empty or contains no sequence data.
    """

    path = os.fspath(path)
    if not os.path.exists(path):
        raise FileNotFoundError(f"No such file: {path!r}")
    if os.path.isdir(path):
        raise FastaParseError(f"Expected a file but got a directory: {path!r}")

    with open(path, "r", encoding="utf-8") as handle:
        text = handle.read()

    records = parse_fasta(text)
    if not records:
        raise FastaParseError(
            f"No sequences found in {path!r}. Is this a FASTA/DNA file?"
        )
    for record in records:
        if not record.sequence:
            raise FastaParseError(
                f"Sequence {record.id!r} in {path!r} is empty."
            )
    return records


def write_fasta(
    records: Iterable[SequenceRecord],
    path: str | os.PathLike[str] | None = None,
    width: int = 70,
) -> str:
    """Serialise ``records`` as FASTA text, optionally writing to ``path``.

    Returns the FASTA text regardless of whether ``path`` is given.
    """

    lines: List[str] = []
    for record in records:
        header = record.header if record.header else record.id
        lines.append(f">{header}")
        seq = record.sequence
        if width and width > 0:
            for i in range(0, len(seq), width):
                lines.append(seq[i : i + width])
        else:
            lines.append(seq)
    text = "\n".join(lines) + "\n"
    if path is not None:
        with open(os.fspath(path), "w", encoding="utf-8") as handle:
            handle.write(text)
    return text


def iter_sequences(path: str | os.PathLike[str]) -> Iterator[SequenceRecord]:
    """Yield records one at a time (thin wrapper over :func:`read_sequences`)."""
    yield from read_sequences(path)
