"""Profile frequency estimation and pseudocount regularization for multiple sequence alignments.

This module provides data structures and mathematical routines to construct position-specific
amino acid frequency profiles from multiple sequence alignments (MSAs).

Mathematical Formulation:
    For an alignment of M homologous sequences of uniform length L:
    1. Observed amino acid counts at column i (0 <= i < L):
           c_{i, a} = sum_{k=1}^M I(s_{k, i} == a)    for a in standard amino acids Sigma
    2. Effective non-gap observation depth:
           N_i = sum_{a in Sigma} c_{i, a}
    3. Gap fraction:
           f_{gap, i} = (count of '-' or '.' in column i) / M
    4. Pseudocount-regularized frequencies (Laplace / Dirichlet mixture prior):
           P_{i, a} = (c_{i, a} + alpha * q_a) / (N_i + alpha)
       where:
           - alpha >= 0 is the pseudocount weight (strength of background prior).
           - q_a > 0 is the background probability of amino acid a (sum q_a = 1.0).

Gap and Ambiguity Policies:
    - Gaps ('-', '.') are indels, not amino acid substitutions. They are recorded in
      gap_fractions for structural analysis but excluded from the residue count denominator N_i.
    - If a column consists entirely of gaps (N_i = 0), the smoothed probability distribution
      defaults smoothly to the background distribution q_a.
    - Ambiguous IUPAC residues:
          * 'B' contributes 0.5 to 'N' (Asparagine) and 0.5 to 'D' (Aspartate).
          * 'Z' contributes 0.5 to 'Q' (Glutamine) and 0.5 to 'E' (Glutamate).
          * 'X' distributes proportional to background probabilities q_a.
          * Other non-amino acid characters raise a ValueError.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Set, Tuple, Union
from .scoring import AMINO_ACIDS, is_valid_amino_acid


# Standard proteome-wide background amino acid frequencies (Robinson & Robinson, 1991), normalized to sum to 1.0
_RAW_ROBINSON: Dict[str, float] = {
    "A": 0.078, "R": 0.051, "N": 0.045, "D": 0.054, "C": 0.019,
    "Q": 0.043, "E": 0.063, "G": 0.074, "H": 0.022, "I": 0.052,
    "L": 0.091, "K": 0.057, "M": 0.023, "F": 0.040, "P": 0.052,
    "S": 0.071, "T": 0.058, "W": 0.013, "Y": 0.032, "V": 0.064,
}
_ROB_SUM = sum(_RAW_ROBINSON.values())
ROBINSON_FREQUENCIES: Dict[str, float] = {aa: _RAW_ROBINSON[aa] / _ROB_SUM for aa in AMINO_ACIDS}

UNIFORM_FREQUENCIES: Dict[str, float] = {aa: 1.0 / len(AMINO_ACIDS) for aa in AMINO_ACIDS}

GAP_CHARACTERS: Set[str] = {"-", "."}


def validate_background_distribution(background: Dict[str, float]) -> None:
    """Validates that a background amino acid distribution is positive and normalized.

    Args:
        background: Mapping from standard amino acid characters to background probabilities.

    Raises:
        ValueError: If any standard amino acid is missing, non-positive, or probabilities do not sum to 1.0.
    """
    for aa in AMINO_ACIDS:
        if aa not in background:
            raise ValueError(f"Background distribution missing standard amino acid '{aa}'.")
        if background[aa] <= 0:
            raise ValueError(f"Background frequency for amino acid '{aa}' must be strictly positive (> 0).")

    total_prob = sum(background[aa] for aa in AMINO_ACIDS)
    if abs(total_prob - 1.0) > 0.01:
        raise ValueError(
            f"Background distribution probabilities must sum to 1.0 (found {total_prob:.4f})."
        )


@dataclass(frozen=True)
class Profile:
    """Represents a position-specific amino acid frequency profile constructed from an MSA.

    Attributes:
        name: Identifier or title of the protein family profile.
        alignment_length: Length L of the aligned columns.
        num_sequences: Number of sequences contributing to this profile.
        consensus: Consensus amino acid string derived from highest-frequency residues.
        counts: List of length L mapping each amino acid to its observed count.
        frequencies: List of length L mapping each amino acid to its regularized probability.
        gap_fractions: List of length L specifying fraction of gaps at each column.
        sequence_ids: Identifiers of sequences used to build this profile (for data leakage auditing).
        background: Background amino acid distribution used for pseudocount smoothing.
        pseudocount_weight: Smoothing parameter alpha applied during frequency calculation.
    """

    name: str
    alignment_length: int
    num_sequences: int
    consensus: str
    counts: List[Dict[str, float]]
    frequencies: List[Dict[str, float]]
    gap_fractions: List[float]
    sequence_ids: List[str] = field(default_factory=list)
    background: Dict[str, float] = field(default_factory=lambda: dict(ROBINSON_FREQUENCIES))
    pseudocount_weight: float = 1.0

    def __len__(self) -> int:
        return self.alignment_length

    def get_column_counts(self, pos: int) -> Dict[str, int]:
        """Returns residue counts at column pos (0 <= pos < alignment_length)."""
        if not (0 <= pos < self.alignment_length):
            raise IndexError(f"Column index {pos} out of range [0, {self.alignment_length}).")
        return self.counts[pos]

    def get_column_frequencies(self, pos: int) -> Dict[str, float]:
        """Returns regularized amino acid frequencies at column pos."""
        if not (0 <= pos < self.alignment_length):
            raise IndexError(f"Column index {pos} out of range [0, {self.alignment_length}).")
        return self.frequencies[pos]

    def get_gap_fraction(self, pos: int) -> float:
        """Returns fraction of gap characters at column pos."""
        if not (0 <= pos < self.alignment_length):
            raise IndexError(f"Column index {pos} out of range [0, {self.alignment_length}).")
        return self.gap_fractions[pos]


def extract_alignment_columns(
    sequences: Sequence[str],
    allow_ambiguous: bool = True,
) -> List[List[str]]:
    """Transposes aligned sequences into a list of column residue lists with strict validation.

    Args:
        sequences: Sequence of aligned strings of identical length.
        allow_ambiguous: If True, permits ambiguous codes (B, Z, X).

    Returns:
        List of length L containing lists of characters at each alignment column.

    Raises:
        ValueError: If sequences is empty, sequence lengths differ, or invalid characters are present.
    """
    if not sequences:
        raise ValueError("Cannot extract columns from an empty sequence collection.")

    col_len = len(sequences[0])
    if col_len == 0:
        raise ValueError("Aligned sequences must have length > 0.")

    for idx, seq in enumerate(sequences):
        if len(seq) != col_len:
            raise ValueError(
                f"Inconsistent sequence lengths in alignment: sequence 0 has length {col_len}, "
                f"but sequence {idx} has length {len(seq)}."
            )
        for char_idx, char in enumerate(seq):
            if char in GAP_CHARACTERS:
                continue
            if not is_valid_amino_acid(char, allow_ambiguous=allow_ambiguous):
                raise ValueError(
                    f"Invalid character '{char}' at index {char_idx} in sequence {idx}. "
                    f"Expected standard IUPAC amino acid residue or gap ('-', '.')."
                )

    return [[seq[c] for seq in sequences] for c in range(col_len)]


def compute_column_counts(
    columns: List[List[str]],
    allow_ambiguous: bool = True,
    background: Optional[Dict[str, float]] = None,
) -> List[Dict[str, int]]:
    """Calculates observed amino acid counts for each column according to the documented gap policy.

    Args:
        columns: List of columns, where each column is a list of characters across sequences.
        allow_ambiguous: If True, handles ambiguous residues (B, Z, X).
        background: Optional background distribution used for proportional X distribution.

    Returns:
        List of dictionaries mapping standard amino acid symbols to observed counts.
    """
    if background is None:
        background = ROBINSON_FREQUENCIES

    counts: List[Dict[str, int]] = []
    for col in columns:
        col_count: Dict[str, float] = {aa: 0.0 for aa in AMINO_ACIDS}
        for char in col:
            u = char.upper()
            if u in GAP_CHARACTERS:
                continue
            if u in col_count:
                col_count[u] += 1.0
            elif allow_ambiguous:
                if u == "B":  # Asx -> 50% Asn, 50% Asp
                    col_count["N"] += 0.5
                    col_count["D"] += 0.5
                elif u == "Z":  # Glx -> 50% Gln, 50% Glu
                    col_count["Q"] += 0.5
                    col_count["E"] += 0.5
                elif u == "X":  # Unknown -> proportional to background
                    for aa in AMINO_ACIDS:
                        col_count[aa] += background.get(aa, 0.05)

        # Preserve exact counts including fractional ambiguous contributions
        counts.append({aa: round(col_count[aa], 3) for aa in AMINO_ACIDS})
    return counts


def compute_profile_frequencies(
    counts: List[Dict[str, Union[int, float]]],
    pseudocount_weight: float = 1.0,
    background: Optional[Dict[str, float]] = None,
) -> List[Dict[str, float]]:
    """Regularizes column counts into frequency probability vectors using background pseudocounts.

    Recurrence:
        P_{i, a} = (c_{i, a} + alpha * q_a) / (N_i + alpha)

    Args:
        counts: List of length L containing observed counts per amino acid.
        pseudocount_weight: alpha >= 0 smoothing strength parameter.
        background: Background distribution q_a.

    Returns:
        List of length L containing probability distributions summing to 1.0.

    Raises:
        ValueError: If pseudocount_weight < 0.
    """
    if pseudocount_weight < 0:
        raise ValueError(f"Pseudocount weight must be non-negative (received {pseudocount_weight}).")

    if background is None:
        background = ROBINSON_FREQUENCIES
    validate_background_distribution(background)

    frequencies: List[Dict[str, float]] = []
    for col_dict in counts:
        n_obs = sum(col_dict.values())
        denominator = n_obs + pseudocount_weight

        col_freq: Dict[str, float] = {}
        if denominator <= 0:
            # Fallback for gap-only column with zero pseudocount
            for aa in AMINO_ACIDS:
                col_freq[aa] = background[aa]
        else:
            for aa in AMINO_ACIDS:
                numerator = col_dict[aa] + pseudocount_weight * background[aa]
                col_freq[aa] = numerator / denominator

        # Ensure strict normalization
        total_p = sum(col_freq.values())
        if total_p > 0:
            for aa in AMINO_ACIDS:
                col_freq[aa] = col_freq[aa] / total_p

        frequencies.append(col_freq)

    return frequencies


def build_profile_from_msa(
    aligned_sequences: Union[Dict[str, str], Sequence[str]],
    name: str = "profile",
    pseudocount_weight: float = 1.0,
    background: Optional[Dict[str, float]] = None,
    sequence_ids: Optional[List[str]] = None,
) -> Profile:
    """Builds a complete, validated Profile instance from an aligned multiple sequence collection.

    Args:
        aligned_sequences: Either a Dict[header, sequence] or a Sequence[sequence].
        name: Name of the resulting profile.
        pseudocount_weight: Strength of pseudocount regularization (alpha).
        background: Background amino acid distribution (defaults to Robinson & Robinson).
        sequence_ids: Optional list of sequence identifiers if passing raw sequence strings.

    Returns:
        Profile instance containing counts, smoothed frequencies, gap fractions, and consensus.

    Raises:
        ValueError: If sequences is empty, lengths differ, or residues are invalid.
    """
    if background is None:
        background = ROBINSON_FREQUENCIES
    validate_background_distribution(background)

    if isinstance(aligned_sequences, dict):
        ids = list(aligned_sequences.keys())
        seqs = list(aligned_sequences.values())
    else:
        seqs = list(aligned_sequences)
        ids = sequence_ids or [f"seq_{i}" for i in range(len(seqs))]

    if not seqs:
        raise ValueError("Cannot build profile from empty multiple sequence alignment.")

    cols = extract_alignment_columns(seqs)
    col_len = len(cols)
    num_seqs = len(seqs)

    counts = compute_column_counts(cols, background=background)
    frequencies = compute_profile_frequencies(
        counts, pseudocount_weight=pseudocount_weight, background=background
    )

    gap_fractions = []
    consensus_chars = []

    for col_idx, col in enumerate(cols):
        gap_count = sum(1 for ch in col if ch in GAP_CHARACTERS)
        gap_fractions.append(gap_count / num_seqs)

        # Consensus is the residue with maximum regularized frequency
        best_aa = max(AMINO_ACIDS, key=lambda a: frequencies[col_idx][a])
        consensus_chars.append(best_aa)

    return Profile(
        name=name,
        alignment_length=col_len,
        num_sequences=num_seqs,
        consensus="".join(consensus_chars),
        counts=counts,
        frequencies=frequencies,
        gap_fractions=gap_fractions,
        sequence_ids=ids,
        background=background,
        pseudocount_weight=pseudocount_weight,
    )


# --- Scientific Rigor & Data Leakage Prevention ---


class DataLeakageError(ValueError):
    """Raised when an evaluation target sequence was included in the profile training set."""
    pass


def split_msa_train_test(
    msa: Dict[str, str],
    test_id: str,
) -> Tuple[Dict[str, str], str]:
    """Splits an MSA into a training MSA excluding the test target, and the un-gapped test sequence.

    Args:
        msa: Mapping from sequence header to aligned sequence.
        test_id: Header of the held-out test sequence.

    Returns:
        Tuple of (train_msa_dict, un-gapped test_sequence).

    Raises:
        KeyError: If test_id is not found in msa.
        ValueError: If msa contains fewer than 2 sequences.
    """
    if test_id not in msa:
        raise KeyError(f"Held-out test target '{test_id}' not found in MSA. Available: {list(msa.keys())}")
    if len(msa) < 2:
        raise ValueError("Cannot perform train/test split on an MSA with fewer than 2 sequences.")

    train_msa = {k: v for k, v in msa.items() if k != test_id}
    test_seq = msa[test_id].replace("-", "").replace(".", "")
    return train_msa, test_seq


def validate_no_data_leakage(
    profile: Profile,
    query_id: Optional[str] = None,
    query_seq: Optional[str] = None,
    allow_leakage: bool = False,
) -> None:
    """Detects data leakage between profile training sequences and an evaluation query.

    Args:
        profile: Profile instance built from an MSA.
        query_id: Optional identifier of the query sequence being aligned.
        query_seq: Optional string of the unaligned query sequence.
        allow_leakage: If True, bypasses exception (for demonstration scripts only).

    Raises:
        DataLeakageError: If query_id matches any profile sequence ID or query_seq is identical.
    """
    if allow_leakage:
        return

    if query_id and query_id in profile.sequence_ids:
        raise DataLeakageError(
            f"DATA LEAKAGE DETECTED: Query sequence ID '{query_id}' was used to construct "
            f"profile '{profile.name}'. Held-out benchmark evaluation requires excluding the "
            f"target from profile training (use split_msa_train_test)."
        )
