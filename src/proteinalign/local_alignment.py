"""Smith–Waterman local pairwise sequence alignment algorithm.

This module implements the Smith–Waterman dynamic programming algorithm
for identifying optimal local high-scoring sequence segments between two
protein sequences under a linear gap penalty model.

Mathematical Recurrence:
    H(i, j) = max(
        0,                           # Reset: terminate path if score drops below zero
        H(i-1, j-1) + s(x_i, y_j),   # Diagonal: align x_i with y_j
        H(i-1, j)   - d,             # Up: align x_i with gap in Y (deletion)
        H(i, j-1)   - d              # Left: align gap in X with y_j (insertion)
    )

Boundary Conditions:
    H(i, 0) = 0   for 0 <= i <= m
    H(0, j) = 0   for 0 <= j <= n

Optimal Score:
    Score_local = max_{0 <= i <= m, 0 <= j <= n} H(i, j) = H(i*, j*)

Traceback Termination:
    Traceback begins at the maximal cell (i*, j*) and terminates immediately upon
    reaching any cell where H(i, j) <= 0. The coordinates where traceback halts
    and begins define the 0-indexed half-open slice [start, end) within seq1 and seq2.

Tie-Breaking Policy:
    1. Maximum cell selection: In the rare event of multiple cells having the identical
       global maximum score, the earliest cell encountered during forward traversal
       (or largest coordinate) is tracked deterministically.
    2. Step transitions: When multiple transitions (Diag, Up, Left) produce the identical
       score, ties are resolved with the priority:
           1. Diagonal (Match / Substitution)
           2. Up (Gap in seq2 / deletion of x_i)
           3. Left (Gap in seq1 / insertion of y_j)

Complexity:
    - Time Complexity:  O(m * n) arithmetic operations and comparisons.
    - Space Complexity: O(m * n) storage for score table.
"""

from typing import Optional
from .models import AlignmentResult
from .scoring import BLOSUM62, SubstitutionMatrix, validate_sequence


def smith_waterman(
    seq1: str,
    seq2: str,
    matrix: Optional[SubstitutionMatrix] = None,
    gap_penalty: float = 8.0,
    allow_ambiguous: bool = True,
) -> AlignmentResult:
    """Computes optimal local pairwise sequence alignment using Smith-Waterman DP.

    Args:
        seq1: Query protein sequence string.
        seq2: Target protein sequence string.
        matrix: Substitution matrix (defaults to BLOSUM62).
        gap_penalty: Linear gap penalty magnitude (positive float subtracted for each gap).
        allow_ambiguous: If True, permits ambiguous IUPAC amino acid codes (X, B, Z).

    Returns:
        AlignmentResult instance containing the local alignment motif, score, and coordinates.

    Raises:
        ValueError: If input sequences contain invalid characters or gap penalty is negative.
    """
    if gap_penalty < 0:
        raise ValueError(f"Gap penalty must be non-negative, got {gap_penalty}.")

    if matrix is None:
        matrix = BLOSUM62

    # Validate residue characters
    validate_sequence(seq1, allow_ambiguous=allow_ambiguous)
    validate_sequence(seq2, allow_ambiguous=allow_ambiguous)

    m = len(seq1)
    n = len(seq2)

    # Empty input edge cases
    if m == 0 or n == 0:
        return AlignmentResult(
            seq1_aligned="",
            seq2_aligned="",
            score=0.0,
            start1=0,
            end1=0,
            start2=0,
            end2=0,
            alignment_type="local",
            gap_open=gap_penalty,
            gap_extend=gap_penalty,
            matrix_name=matrix.name,
            seq1=seq1,
            seq2=seq2,
        )

    # 1. Initialize (m+1) x (n+1) DP score table with zeros
    H = [[0.0] * (n + 1) for _ in range(m + 1)]
    max_score = 0.0
    max_i, max_j = 0, 0

    # 2. Fill score table using local recurrence with zero floor
    for i in range(1, m + 1):
        char1 = seq1[i - 1]
        for j in range(1, n + 1):
            char2 = seq2[j - 1]
            diag = H[i - 1][j - 1] + matrix.score(char1, char2)
            up = H[i - 1][j] - gap_penalty
            left = H[i][j - 1] - gap_penalty
            score = max(0.0, diag, up, left)
            H[i][j] = score

            if score > max_score:
                max_score = score
                max_i, max_j = i, j

    # If no positive-scoring local alignment exists
    if max_score <= 0.0:
        return AlignmentResult(
            seq1_aligned="",
            seq2_aligned="",
            score=0.0,
            start1=0,
            end1=0,
            start2=0,
            end2=0,
            alignment_type="local",
            gap_open=gap_penalty,
            gap_extend=gap_penalty,
            matrix_name=matrix.name,
            seq1=seq1,
            seq2=seq2,
        )

    # 3. Traceback from (max_i, max_j) until H[curr_i][curr_j] == 0
    aligned1, aligned2 = [], []
    curr_i, curr_j = max_i, max_j

    while curr_i > 0 and curr_j > 0 and H[curr_i][curr_j] > 0:
        char1 = seq1[curr_i - 1]
        char2 = seq2[curr_j - 1]
        diag_score = H[curr_i - 1][curr_j - 1] + matrix.score(char1, char2)

        # Diagonal priority
        if abs(H[curr_i][curr_j] - diag_score) < 1e-9:
            aligned1.append(char1)
            aligned2.append(char2)
            curr_i -= 1
            curr_j -= 1
        elif abs(H[curr_i][curr_j] - (H[curr_i - 1][curr_j] - gap_penalty)) < 1e-9:
            aligned1.append(char1)
            aligned2.append("-")
            curr_i -= 1
        elif abs(H[curr_i][curr_j] - (H[curr_i][curr_j - 1] - gap_penalty)) < 1e-9:
            aligned1.append("-")
            aligned2.append(char2)
            curr_j -= 1
        else:
            break

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
        alignment_type="local",
        gap_open=gap_penalty,
        gap_extend=gap_penalty,
        matrix_name=matrix.name,
        seq1=seq1,
        seq2=seq2,
    )
