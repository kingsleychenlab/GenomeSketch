"""Translation of DNA/RNA into protein using the standard genetic code.

The codon table is built manually (NCBI translation table 1, the standard
code). Codons containing ``N`` or any unknown base translate to ``X``; stop
codons translate to ``*``.
"""

from __future__ import annotations

from typing import Dict

from genomesketch.transform import reverse_complement
from genomesketch.validation import normalize

# NCBI standard genetic code (transl_table=1), written against DNA codons.
STANDARD_CODON_TABLE: Dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

START_CODONS = {"ATG"}
STOP_CODONS = {"TAA", "TAG", "TGA"}


def translate_codon(codon: str) -> str:
    """Translate a single 3-letter codon.

    ``U`` is accepted (RNA) and treated as ``T``. Any codon containing an
    unknown base (including ``N``) returns ``X``.
    """
    codon = normalize(codon).replace("U", "T")
    if len(codon) != 3:
        raise ValueError(f"Codon must be exactly 3 bases, got {codon!r}")
    return STANDARD_CODON_TABLE.get(codon, "X")


def translate(
    sequence: str,
    frame: int = 1,
    *,
    to_stop: bool = False,
) -> str:
    """Translate ``sequence`` in the given reading ``frame``.

    Frames ``1, 2, 3`` read the forward strand starting at offset ``0, 1, 2``.
    Frames ``-1, -2, -3`` read the reverse complement with the same offsets.
    Any incomplete trailing codon is ignored. When ``to_stop`` is ``True`` the
    protein is truncated at (and excludes) the first stop codon.
    """

    if frame not in (1, 2, 3, -1, -2, -3):
        raise ValueError("frame must be one of 1, 2, 3, -1, -2, -3")

    seq = normalize(sequence).replace("U", "T")
    if frame < 0:
        seq = reverse_complement(seq)
    offset = abs(frame) - 1
    seq = seq[offset:]

    residues = []
    for i in range(0, len(seq) - len(seq) % 3, 3):
        codon = seq[i : i + 3]
        aa = STANDARD_CODON_TABLE.get(codon, "X")
        if to_stop and aa == "*":
            break
        residues.append(aa)
    return "".join(residues)


def six_frame_translation(sequence: str) -> Dict[int, str]:
    """Return a mapping ``{frame: protein}`` for all six reading frames."""
    return {frame: translate(sequence, frame) for frame in (1, 2, 3, -1, -2, -3)}
