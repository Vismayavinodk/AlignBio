from typing import Dict, List, Optional, Set, Tuple
from .scoring import BLOSUM62, SubstitutionMatrix


def compute_pid(aligned1: str, aligned2: str, mode: str = "aligned_columns") -> float:
    """Computes percent sequence identity between two aligned sequences.

    Args:
        aligned1: First aligned sequence string (with '-').
        aligned2: Second aligned sequence string (with '-').
        mode: 'aligned_columns' (total columns) or 'non_gaps' (only columns without gaps).

    Returns:
        Percentage of identical residues (0.0 to 100.0).
    """
    if len(aligned1) != len(aligned2) or len(aligned1) == 0:
        return 0.0

    matches = 0
    valid_cols = 0

    for a, b in zip(aligned1, aligned2):
        if mode == "non_gaps":
            if a != "-" and b != "-":
                valid_cols += 1
                if a.upper() == b.upper():
                    matches += 1
        else:  # aligned_columns
            valid_cols += 1
            if a.upper() == b.upper() and a != "-":
                matches += 1

    return (matches / valid_cols) * 100.0 if valid_cols > 0 else 0.0


def extract_aligned_pairs(seq1: str, seq2: str, uppercase_only: bool = True) -> Set[Tuple[int, int]]:
    """Extracts the set of aligned residue index pairs between two sequences.

    Indexes returned are (orig_idx1, orig_idx2) based on 0-indexed positions in ungapped sequences.
    If uppercase_only is True, only columns where both residues are uppercase (core conserved) are kept.
    """
    pairs = set()
    idx1 = 0
    idx2 = 0

    for c1, c2 in zip(seq1, seq2):
        is_res1 = c1 != "-"
        is_res2 = c2 != "-"

        if is_res1 and is_res2:
            if not uppercase_only or (c1.isupper() and c2.isupper()):
                pairs.add((idx1, idx2))

        if is_res1:
            idx1 += 1
        if is_res2:
            idx2 += 1

    return pairs


def compute_q_score(test_alignment: Dict[str, str], ref_alignment: Dict[str, str], uppercase_only: bool = True) -> float:
    """Computes benchmark Q-score (developer score) comparing a test alignment against a reference alignment.

    Q = (number of correctly aligned residue pairs) / (total aligned pairs in reference)
    """
    common_keys = sorted(set(test_alignment.keys()) & set(ref_alignment.keys()))
    if len(common_keys) < 2:
        return 0.0

    total_ref_pairs = 0
    correct_test_pairs = 0

    for i in range(len(common_keys)):
        for j in range(i + 1, len(common_keys)):
            k1, k2 = common_keys[i], common_keys[j]

            ref_pairs = extract_aligned_pairs(ref_alignment[k1], ref_alignment[k2], uppercase_only=uppercase_only)
            test_pairs = extract_aligned_pairs(test_alignment[k1], test_alignment[k2], uppercase_only=False)

            total_ref_pairs += len(ref_pairs)
            correct_test_pairs += len(ref_pairs & test_pairs)

    return (correct_test_pairs / total_ref_pairs) if total_ref_pairs > 0 else 0.0


def compute_sp_score(alignment: Dict[str, str], matrix: Optional[SubstitutionMatrix] = None, gap_penalty: float = 8.0) -> float:
    """Computes Sum-of-Pairs (SP) score across all sequence pairs in an MSA."""
    if matrix is None:
        matrix = BLOSUM62

    keys = list(alignment.keys())
    if len(keys) < 2:
        return 0.0

    total_score = 0.0
    num_cols = len(alignment[keys[0]])

    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            s1 = alignment[keys[i]]
            s2 = alignment[keys[j]]
            for c1, c2 in zip(s1, s2):
                if c1 == "-" and c2 == "-":
                    continue
                elif c1 == "-" or c2 == "-":
                    total_score -= gap_penalty
                else:
                    total_score += matrix.score(c1, c2)

    return total_score
