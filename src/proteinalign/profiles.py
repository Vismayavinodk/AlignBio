"""Profile frequency estimation and pseudocount regularization for multiple sequence alignments."""

from typing import Dict, List, Optional, Sequence
from .scoring import AMINO_ACIDS

# Standard proteome-wide background amino acid frequencies (Robinson & Robinson, 1991)
ROBINSON_FREQUENCIES: Dict[str, float] = {
    "A": 0.078, "R": 0.051, "N": 0.045, "D": 0.054, "C": 0.019,
    "Q": 0.043, "E": 0.063, "G": 0.074, "H": 0.022, "I": 0.052,
    "L": 0.091, "K": 0.057, "M": 0.023, "F": 0.040, "P": 0.052,
    "S": 0.071, "T": 0.058, "W": 0.013, "Y": 0.032, "V": 0.064,
}

UNIFORM_FREQUENCIES: Dict[str, float] = {aa: 1.0 / len(AMINO_ACIDS) for aa in AMINO_ACIDS}


def extract_alignment_columns(sequences: Sequence[str]) -> List[List[str]]:
    """Transposes aligned sequences into list of column residue lists."""
    if not sequences:
        return []
    col_len = len(sequences[0])
    for s in sequences:
        if len(s) != col_len:
            raise ValueError("All sequences in alignment must have equal length.")

    return [[seq[c] for seq in sequences] for c in range(col_len)]


def compute_column_counts(columns: List[List[str]]) -> List[Dict[str, int]]:
    """Calculates observed residue counts for each column, ignoring gap characters."""
    counts = []
    for col in columns:
        col_count = {aa: 0 for aa in AMINO_ACIDS}
        for char in col:
            u = char.upper()
            if u in col_count:
                col_count[u] += 1
        counts.append(col_count)
    return counts


def compute_profile_frequencies(
    counts: List[Dict[str, int]],
    pseudocount_weight: float = 1.0,
    background: Optional[Dict[str, float]] = None,
) -> List[Dict[str, float]]:
    """Regularizes column counts into frequency vectors using background pseudocounts.

    f_i(a) = (c_i(a) + B * p_a) / (N_i + B)
    """
    if background is None:
        background = ROBINSON_FREQUENCIES

    frequencies = []
    for col_dict in counts:
        n_obs = sum(col_dict.values())
        denominator = n_obs + pseudocount_weight
        col_freq = {}
        for aa in AMINO_ACIDS:
            numerator = col_dict[aa] + pseudocount_weight * background.get(aa, 0.05)
            col_freq[aa] = numerator / denominator if denominator > 0 else background.get(aa, 0.05)
        frequencies.append(col_freq)
    return frequencies
