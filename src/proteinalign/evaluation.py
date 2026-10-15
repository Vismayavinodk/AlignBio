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


def compute_modeler_score(test_alignment: Dict[str, str], ref_alignment: Dict[str, str], uppercase_only: bool = True) -> float:
    """Computes benchmark Modeler score (precision): fraction of test aligned residue pairs present in the reference.

    M = (number of correctly aligned residue pairs) / (total aligned pairs in test)
    """
    common_keys = sorted(set(test_alignment.keys()) & set(ref_alignment.keys()))
    if len(common_keys) < 2:
        return 0.0

    total_test_pairs = 0
    correct_test_pairs = 0

    for i in range(len(common_keys)):
        for j in range(i + 1, len(common_keys)):
            k1, k2 = common_keys[i], common_keys[j]
            ref_pairs = extract_aligned_pairs(ref_alignment[k1], ref_alignment[k2], uppercase_only=uppercase_only)
            test_pairs = extract_aligned_pairs(test_alignment[k1], test_alignment[k2], uppercase_only=False)

            total_test_pairs += len(test_pairs)
            correct_test_pairs += len(ref_pairs & test_pairs)

    return (correct_test_pairs / total_test_pairs) if total_test_pairs > 0 else 0.0


def compute_f1_score(test_alignment: Dict[str, str], ref_alignment: Dict[str, str], uppercase_only: bool = True) -> float:
    """Computes harmonic mean between Q-score (sensitivity/recall) and Modeler score (precision)."""
    q = compute_q_score(test_alignment, ref_alignment, uppercase_only=uppercase_only)
    m = compute_modeler_score(test_alignment, ref_alignment, uppercase_only=uppercase_only)
    if q + m <= 0.0:
        return 0.0
    return 2.0 * (q * m) / (q + m)


def compute_gap_statistics(aligned_seq: str) -> Dict[str, float]:
    """Computes total gap characters, gap opening count, and mean indel run length."""
    import re
    gap_runs = re.findall(r"-+", aligned_seq)
    total_gaps = sum(len(g) for g in gap_runs)
    num_runs = len(gap_runs)
    mean_length = (total_gaps / num_runs) if num_runs > 0 else 0.0
    max_length = max((len(g) for g in gap_runs), default=0)
    return {
        "gap_characters": float(total_gaps),
        "gap_openings": float(num_runs),
        "mean_gap_length": round(mean_length, 3),
        "max_gap_length": float(max_length),
    }


def compute_alignment_gap_stats(aligned1: str, aligned2: str) -> Dict[str, float]:
    """Computes combined gap statistics for a pairwise alignment."""
    stats1 = compute_gap_statistics(aligned1)
    stats2 = compute_gap_statistics(aligned2)
    total_gaps = stats1["gap_characters"] + stats2["gap_characters"]
    total_openings = stats1["gap_openings"] + stats2["gap_openings"]
    mean_length = (total_gaps / total_openings) if total_openings > 0 else 0.0
    return {
        "total_gaps": total_gaps,
        "gap_openings": total_openings,
        "mean_gap_length": round(mean_length, 3),
        "max_gap_length": max(stats1["max_gap_length"], stats2["max_gap_length"]),
    }


def recompute_score(
    aligned1: str,
    aligned2: str,
    matrix: Optional[SubstitutionMatrix] = None,
    gap_open: float = 8.0,
    gap_extend: Optional[float] = None,
) -> float:
    """Recomputes alignment score directly from aligned sequence strings.

    Supports both linear (when gap_extend is None or equal to gap_open)
    and Gotoh affine gap penalties (when gap_extend is specified).
    """
    if matrix is None:
        matrix = BLOSUM62

    if len(aligned1) != len(aligned2):
        raise ValueError("Aligned sequences must have identical lengths.")

    is_affine = gap_extend is not None
    score = 0.0

    if not is_affine:
        for a, b in zip(aligned1, aligned2):
            if a == "-" and b == "-":
                continue
            elif a == "-" or b == "-":
                score -= gap_open
            else:
                score += matrix.score(a, b)
        return score

    for a, b in zip(aligned1, aligned2):
        if a != "-" and b != "-":
            score += matrix.score(a, b)

    import re
    for match in re.finditer(r"-+", aligned1):
        k = len(match.group())
        score -= (gap_open + (k - 1) * gap_extend)

    for match in re.finditer(r"-+", aligned2):
        k = len(match.group())
        score -= (gap_open + (k - 1) * gap_extend)

    return score


def compute_pairwise_metrics(
    test_aln1: str,
    test_aln2: str,
    ref_aln1: str,
    ref_aln2: str,
    uppercase_only: bool = True,
) -> Dict[str, float]:
    """Computes comprehensive alignment metrics for a single pair of sequences."""
    test_dict = {"seq1": test_aln1, "seq2": test_aln2}
    ref_dict = {"seq1": ref_aln1, "seq2": ref_aln2}

    q = compute_q_score(test_dict, ref_dict, uppercase_only=uppercase_only)
    modeler = compute_modeler_score(test_dict, ref_dict, uppercase_only=uppercase_only)
    f1 = compute_f1_score(test_dict, ref_dict, uppercase_only=uppercase_only)
    pid = compute_pid(test_aln1, test_aln2, mode="aligned_columns")
    gap_stats = compute_alignment_gap_stats(test_aln1, test_aln2)

    return {
        "q_score": round(q, 4),
        "modeler_score": round(modeler, 4),
        "f1_score": round(f1, 4),
        "pid": round(pid, 2),
        "total_gaps": gap_stats["total_gaps"],
        "gap_openings": gap_stats["gap_openings"],
        "mean_gap_length": gap_stats["mean_gap_length"],
        "alignment_length": float(len(test_aln1)),
    }


def compute_profile_placement_accuracy(
    aligned_consensus: str,
    aligned_query: str,
    ref_target_aln: str,
) -> Dict[str, float]:
    """Computes placement accuracy of held-out query residues against reference MSA columns.

    Args:
        aligned_consensus: Predicted aligned consensus string (from PSSM or matrix alignment).
        aligned_query: Predicted aligned query string.
        ref_target_aln: Ground-truth reference aligned string (from benchmark reference MSA).

    Returns:
        Dictionary with profile_q_score, total_core_residues, and correct_core_residues.
    """
    ref_cols = []
    is_core = []
    for col_idx, char in enumerate(ref_target_aln):
        if char != "-":
            ref_cols.append(col_idx)
            is_core.append(char.isupper())

    pred_cols = []
    curr_cons_idx = 0
    for a, b in zip(aligned_consensus, aligned_query):
        if b != "-":
            if a != "-":
                pred_cols.append(curr_cons_idx)
            else:
                pred_cols.append(-1)
        if a != "-":
            curr_cons_idx += 1

    if len(ref_cols) != len(pred_cols):
        raise ValueError(
            f"Length mismatch: {len(ref_cols)} ungapped residues in ref, {len(pred_cols)} in predicted alignment."
        )

    total_core = sum(1 for c in is_core if c)
    correct_core = 0
    for r_col, p_col, core in zip(ref_cols, pred_cols, is_core):
        if core and r_col == p_col:
            correct_core += 1

    q_score = (correct_core / total_core) if total_core > 0 else 0.0
    return {
        "profile_q_score": round(q_score, 4),
        "total_core_residues": float(total_core),
        "correct_core_residues": float(correct_core),
    }


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
