"""Gotoh 3-state dynamic programming for affine gap penalties.

This module implements the Gotoh (1982) three-state dynamic programming algorithm
for optimal global (Needleman-Wunsch) and local (Smith-Waterman) pairwise sequence
alignment with affine gap penalties.

Gap Cost Convention:
    Under the affine model, a contiguous gap of length k >= 1 incurs a penalty:
        G(k) = gap_open + (k - 1) * gap_extend
    where gap_open > 0 and gap_extend > 0.
    - A gap of length 1 costs exactly gap_open.
    - A gap of length 2 costs gap_open + gap_extend.
    - Each additional gap character costs gap_extend.

Mathematical States:
    - M(i, j): Optimal alignment of prefixes X[1..i] and Y[1..j] ending with a match/mismatch.
    - Ix(i, j): Optimal alignment of prefixes X[1..i] and Y[1..j] ending with a gap in Y (deletion of x_i).
    - Iy(i, j): Optimal alignment of prefixes X[1..i] and Y[1..j] ending with a gap in X (insertion of y_j).

Complexity:
    - Time Complexity: O(m * n) arithmetic operations across 3 coupled DP layers.
    - Space Complexity: O(m * n) storage for score matrices and traceback pointers.
"""

from typing import Optional
from .models import AlignmentResult
from .scoring import BLOSUM62, SubstitutionMatrix, validate_sequence


NEG_INF = -1e9


def needleman_wunsch_affine(
    seq1: str,
    seq2: str,
    matrix: Optional[SubstitutionMatrix] = None,
    gap_open: float = 10.0,
    gap_extend: float = 1.0,
    allow_ambiguous: bool = True,
) -> AlignmentResult:
    """Computes optimal global pairwise sequence alignment with affine gap penalties using Gotoh's 3-state DP algorithm.

    Args:
        seq1: Query protein sequence.
        seq2: Target protein sequence.
        matrix: Substitution matrix (defaults to BLOSUM62).
        gap_open: Penalty magnitude for opening a gap (d > 0).
        gap_extend: Penalty magnitude for extending a gap (e > 0).
        allow_ambiguous: If True, permits ambiguous IUPAC amino acid codes (X, B, Z).

    Returns:
        AlignmentResult with the optimal global affine alignment.

    Raises:
        ValueError: If input sequences contain invalid characters, gap penalties are negative,
                    or gap_open < gap_extend.
    """
    validate_sequence(seq1, allow_ambiguous=allow_ambiguous)
    validate_sequence(seq2, allow_ambiguous=allow_ambiguous)

    if gap_open < 0 or gap_extend < 0:
        raise ValueError("Gap penalties must be non-negative.")
    if gap_open < gap_extend:
        raise ValueError(
            f"Gap opening penalty ({gap_open}) should typically be >= gap extension penalty ({gap_extend})."
        )

    if matrix is None:
        matrix = BLOSUM62

    m, n = len(seq1), len(seq2)

    # Empty cases
    if m == 0 and n == 0:
        return AlignmentResult(
            seq1_aligned="",
            seq2_aligned="",
            score=0.0,
            start1=0,
            end1=0,
            start2=0,
            end2=0,
            alignment_type="affine_global",
            gap_open=gap_open,
            gap_extend=gap_extend,
            matrix_name=matrix.name,
            seq1=seq1,
            seq2=seq2,
        )
    if m == 0:
        cost = -(gap_open + (n - 1) * gap_extend)
        return AlignmentResult(
            seq1_aligned="-" * n,
            seq2_aligned=seq2,
            score=cost,
            start1=0,
            end1=0,
            start2=0,
            end2=n,
            alignment_type="affine_global",
            gap_open=gap_open,
            gap_extend=gap_extend,
            matrix_name=matrix.name,
            seq1=seq1,
            seq2=seq2,
        )
    if n == 0:
        cost = -(gap_open + (m - 1) * gap_extend)
        return AlignmentResult(
            seq1_aligned=seq1,
            seq2_aligned="-" * m,
            score=cost,
            start1=0,
            end1=m,
            start2=0,
            end2=0,
            alignment_type="affine_global",
            gap_open=gap_open,
            gap_extend=gap_extend,
            matrix_name=matrix.name,
            seq1=seq1,
            seq2=seq2,
        )

    M = [[NEG_INF] * (n + 1) for _ in range(m + 1)]
    Ix = [[NEG_INF] * (n + 1) for _ in range(m + 1)]
    Iy = [[NEG_INF] * (n + 1) for _ in range(m + 1)]

    M[0][0] = 0.0

    for i in range(1, m + 1):
        Ix[i][0] = -(gap_open + (i - 1) * gap_extend)
    for j in range(1, n + 1):
        Iy[0][j] = -(gap_open + (j - 1) * gap_extend)

    for i in range(1, m + 1):
        char1 = seq1[i - 1]
        for j in range(1, n + 1):
            char2 = seq2[j - 1]
            sub_score = matrix.score(char1, char2)

            prev_max = max(M[i - 1][j - 1], Ix[i - 1][j - 1], Iy[i - 1][j - 1])
            M[i][j] = sub_score + prev_max
            Ix[i][j] = max(M[i - 1][j] - gap_open, Ix[i - 1][j] - gap_extend)
            Iy[i][j] = max(M[i][j - 1] - gap_open, Iy[i][j - 1] - gap_extend)

    best_score = max(M[m][n], Ix[m][n], Iy[m][n])
    if abs(best_score - M[m][n]) < 1e-9:
        curr_state = 0
    elif abs(best_score - Ix[m][n]) < 1e-9:
        curr_state = 1
    else:
        curr_state = 2

    aligned1, aligned2 = [], []
    curr_i, curr_j = m, n

    while curr_i > 0 or curr_j > 0:
        if curr_i == 0:
            aligned1.append("-")
            aligned2.append(seq2[curr_j - 1])
            curr_j -= 1
            continue
        if curr_j == 0:
            aligned1.append(seq1[curr_i - 1])
            aligned2.append("-")
            curr_i -= 1
            continue

        if curr_state == 0:
            char1 = seq1[curr_i - 1]
            char2 = seq2[curr_j - 1]
            aligned1.append(char1)
            aligned2.append(char2)
            prev_val = M[curr_i][curr_j] - matrix.score(char1, char2)

            if abs(prev_val - M[curr_i - 1][curr_j - 1]) < 1e-9:
                curr_state = 0
            elif abs(prev_val - Ix[curr_i - 1][curr_j - 1]) < 1e-9:
                curr_state = 1
            else:
                curr_state = 2

            curr_i -= 1
            curr_j -= 1

        elif curr_state == 1:
            aligned1.append(seq1[curr_i - 1])
            aligned2.append("-")
            if abs(Ix[curr_i][curr_j] - (M[curr_i - 1][curr_j] - gap_open)) < 1e-9:
                curr_state = 0
            else:
                curr_state = 1
            curr_i -= 1

        elif curr_state == 2:
            aligned1.append("-")
            aligned2.append(seq2[curr_j - 1])
            if abs(Iy[curr_i][curr_j] - (M[curr_i][curr_j - 1] - gap_open)) < 1e-9:
                curr_state = 0
            else:
                curr_state = 2
            curr_j -= 1

    aligned1.reverse()
    aligned2.reverse()

    return AlignmentResult(
        seq1_aligned="".join(aligned1),
        seq2_aligned="".join(aligned2),
        score=best_score,
        start1=0,
        end1=m,
        start2=0,
        end2=n,
        alignment_type="affine_global",
        gap_open=gap_open,
        gap_extend=gap_extend,
        matrix_name=matrix.name,
        seq1=seq1,
        seq2=seq2,
    )


def smith_waterman_affine(
    seq1: str,
    seq2: str,
    matrix: Optional[SubstitutionMatrix] = None,
    gap_open: float = 10.0,
    gap_extend: float = 1.0,
    allow_ambiguous: bool = True,
) -> AlignmentResult:
    """Computes optimal local pairwise sequence alignment with affine gap penalties using 3-state DP.

    Args:
        seq1: Query protein sequence.
        seq2: Target protein sequence.
        matrix: Substitution matrix (defaults to BLOSUM62).
        gap_open: Penalty magnitude for opening a gap (d > 0).
        gap_extend: Penalty magnitude for extending a gap (e > 0).
        allow_ambiguous: If True, permits ambiguous IUPAC amino acid codes (X, B, Z).

    Returns:
        AlignmentResult with the optimal local affine alignment and extracted coordinates.

    Raises:
        ValueError: If input sequences contain invalid characters, gap penalties are negative,
                    or gap_open < gap_extend.
    """
    validate_sequence(seq1, allow_ambiguous=allow_ambiguous)
    validate_sequence(seq2, allow_ambiguous=allow_ambiguous)

    if gap_open < 0 or gap_extend < 0:
        raise ValueError("Gap penalties must be non-negative.")
    if gap_open < gap_extend:
        raise ValueError(
            f"Gap opening penalty ({gap_open}) should typically be >= gap extension penalty ({gap_extend})."
        )

    if matrix is None:
        matrix = BLOSUM62

    m, n = len(seq1), len(seq2)
    if m == 0 or n == 0:
        return AlignmentResult(
            seq1_aligned="",
            seq2_aligned="",
            score=0.0,
            start1=0,
            end1=0,
            start2=0,
            end2=0,
            alignment_type="affine_local",
            gap_open=gap_open,
            gap_extend=gap_extend,
            matrix_name=matrix.name,
            seq1=seq1,
            seq2=seq2,
        )

    M = [[0.0] * (n + 1) for _ in range(m + 1)]
    Ix = [[0.0] * (n + 1) for _ in range(m + 1)]
    Iy = [[0.0] * (n + 1) for _ in range(m + 1)]

    max_score = 0.0
    max_i, max_j = 0, 0

    for i in range(1, m + 1):
        char1 = seq1[i - 1]
        for j in range(1, n + 1):
            char2 = seq2[j - 1]
            sub_score = matrix.score(char1, char2)

            prev_max = max(0.0, M[i - 1][j - 1], Ix[i - 1][j - 1], Iy[i - 1][j - 1])
            m_val = max(0.0, sub_score + prev_max)
            M[i][j] = m_val

            ix_val = max(0.0, M[i - 1][j] - gap_open, Ix[i - 1][j] - gap_extend)
            Ix[i][j] = ix_val

            iy_val = max(0.0, M[i][j - 1] - gap_open, Iy[i][j - 1] - gap_extend)
            Iy[i][j] = iy_val

            if m_val > max_score:
                max_score = m_val
                max_i, max_j = i, j

    if max_score <= 0.0:
        return AlignmentResult(
            seq1_aligned="",
            seq2_aligned="",
            score=0.0,
            start1=0,
            end1=0,
            start2=0,
            end2=0,
            alignment_type="affine_local",
            gap_open=gap_open,
            gap_extend=gap_extend,
            matrix_name=matrix.name,
            seq1=seq1,
            seq2=seq2,
        )

    aligned1, aligned2 = [], []
    curr_i, curr_j = max_i, max_j
    curr_state = 0

    while curr_i > 0 and curr_j > 0:
        if curr_state == 0:
            if M[curr_i][curr_j] <= 0.0:
                break
            char1 = seq1[curr_i - 1]
            char2 = seq2[curr_j - 1]
            prev_val = M[curr_i][curr_j] - matrix.score(char1, char2)

            aligned1.append(char1)
            aligned2.append(char2)

            if prev_val <= 0.0:
                curr_i -= 1
                curr_j -= 1
                break

            if abs(prev_val - M[curr_i - 1][curr_j - 1]) < 1e-9:
                curr_state = 0
            elif abs(prev_val - Ix[curr_i - 1][curr_j - 1]) < 1e-9:
                curr_state = 1
            else:
                curr_state = 2

            curr_i -= 1
            curr_j -= 1

        elif curr_state == 1:
            if Ix[curr_i][curr_j] <= 0.0:
                break
            aligned1.append(seq1[curr_i - 1])
            aligned2.append("-")
            if abs(Ix[curr_i][curr_j] - (M[curr_i - 1][curr_j] - gap_open)) < 1e-9:
                curr_state = 0
            else:
                curr_state = 1
            curr_i -= 1

        elif curr_state == 2:
            if Iy[curr_i][curr_j] <= 0.0:
                break
            aligned1.append("-")
            aligned2.append(seq2[curr_j - 1])
            if abs(Iy[curr_i][curr_j] - (M[curr_i][curr_j - 1] - gap_open)) < 1e-9:
                curr_state = 0
            else:
                curr_state = 2
            curr_j -= 1

    aligned1.reverse()
    aligned2.reverse()

    return AlignmentResult(
        seq1_aligned="".join(aligned1),
        seq2_aligned="".join(aligned2),
        score=max_score,
        start1=curr_i,
        end1=max_i,
        start2=curr_j,
        end2=max_j,
        alignment_type="affine_local",
        gap_open=gap_open,
        gap_extend=gap_extend,
        matrix_name=matrix.name,
        seq1=seq1,
        seq2=seq2,
    )
