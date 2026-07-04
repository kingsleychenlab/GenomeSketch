"""GenomeSketch: a fully offline DNA sequence analysis toolkit.

GenomeSketch parses local FASTA/plain-text DNA files and runs a collection of
bioinformatics algorithms implemented from scratch: nucleotide statistics,
motif search, reverse complement, transcription, translation, ORF finding,
mutation comparison, global alignment, repeat detection, and visualizations.

No network access, no external bioinformatics dependencies for core algorithms.
"""

from genomesketch.io import SequenceRecord, parse_fasta, read_sequences

__all__ = [
    "SequenceRecord",
    "parse_fasta",
    "read_sequences",
    "__version__",
]

__version__ = "0.1.0"
