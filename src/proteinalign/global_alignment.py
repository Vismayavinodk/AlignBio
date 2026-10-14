"""Needleman–Wunsch global pairwise sequence alignment algorithm.

This module implements the classic Needleman–Wunsch dynamic programming algorithm
for finding the optimal end-to-end alignment between two protein sequences under a
linear gap penalty model.

Mathematical Recurrence:
    F(i, j) = max(
        F(i-1, j-1) + s(x_i, y_j),   # Diagonal: align x_i with y_j
        F(i-1, j)   - d,             # Up: align x_i with gap in Y (deletion)
        F(i, j-1)   - d              # Left: align gap in X with y_j (insertion)
    )

Boundary Conditions:
    F(0, 0) = 0
    F(i, 0) = -i * d   for 1 <= i <= m
    F(0, j) = -j * d   for 1 <= j <= n

Tie-Breaking Policy:
    When two or three transitions produce the identical optimal score, ties are
    resolved with the following deterministic priority:
        1. Diagonal (Match / Substitution)
        2. Up (Gap in seq2 / deletion of x_i)
        3. Left (Gap in seq1 / insertion of y_j)
    This standard policy prefers residue pairing over indel introduction and ensures
    reproducibility across runs.

Complexity:
    - Time Complexity:  O(m * n) arithmetic operations and comparisons.
    - Space Complexity: O(m * n) storage for score table and traceback paths.
"""

from typing import Optional
from .models import AlignmentResult
from .scoring import BLOSUM62, SubstitutionMatrix, validate_sequence


def needleman_wunsch(
    seq1: str,
    seq2: str,
    matrix: Optional[SubstitutionMatrix] = None,
    gap_penalty: float = 8.0,
    allow_ambiguous: bool = True,
) -> AlignmentResult:
    """Computes the optimal global pairwise sequence alignment using Needleman-Wunsch DP.

    Args:
        seq1: Query protein sequence string.
        seq2: Target protein sequence string.
        matrix: Substitution matrix (defaults to BLOSUM62).
        gap_penalty: Linear gap penalty magnitude (positive float subtracted for each gap).
        allow_ambiguous: If True, permits ambiguous IUPAC amino acid codes (X, B, Z).

    Returns:
        AlignmentResult instance containing aligned sequences, DP score, coordinates, and metadata.

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

    # Empty sequence boundary cases
    if m == 0 and n == 0:
        return AlignmentResult(
            seq1_aligned="",
            seq2_aligned="",
            score=0.0,
            start1=0,
            end1=0,
            start2=0,
            end2=0,
            alignment_type="global",
            gap_open=gap_penalty,
            gap_extend=gap_penalty,
            matrix_name=matrix.name,
            seq1=seq1,
            seq2=seq2,
        )
    if m == 0:
        return AlignmentResult(
            seq1_aligned="-" * n,
            seq2_aligned=seq2,
            score=-n * gap_penalty,
            start1=0,
            end1=0,
            start2=0,
            end2=n,
            alignment_type="global",
            gap_open=gap_penalty,
            gap_extend=gap_penalty,
            matrix_name=matrix.name,
            seq1=seq1,
            seq2=seq2,
        )
    if n == 0:
        return AlignmentResult(
            seq1_aligned=seq1,
            seq2_aligned="-" * m,
            score=-m * gap_penalty,
            start1=0,
            end1=m,
            start2=0,
            end2=0,
            alignment_type="global",
            gap_open=gap_penalty,
            gap_extend=gap_penalty,
            matrix_name=matrix.name,
            seq1=seq1,
            seq2=seq2,
        )

    # 1. Initialize (m+1) x (n+1) DP score table
    F = [[0.0] * (n + 1) for _ in range(m + 1)]

    # Boundary conditions
    for i in range(1, m + 1):
        F[i][0] = -i * gap_penalty
    for j in range(1, n + 1):
        F[0][j] = -j * gap_penalty

    # 2. Fill matrix using the Needleman-Wunsch recurrence
    for i in range(1, m + 1):
        char1 = seq1[i - 1]
        for j in range(1, n + 1):
            char2 = seq2[j - 1]
            diag = F[i - 1][j - 1] + matrix.score(char1, char2)
            up = F[i - 1][j] - gap_penalty
            left = F[i][j - 1] - gap_penalty
            F[i][j] = max(diag, up, left)

    # 3. Traceback from (m, n) to (0, 0)
    aligned1, aligned2 = [], []
    curr_i, curr_j = m, n

    while curr_i > 0 or curr_j > 0:
        # Check Diagonal (Match/Substitution) with highest priority
        if curr_i > 0 and curr_j > 0:
            char1 = seq1[curr_i - 1]
            char2 = seq2[curr_j - 1]
            diag_score = F[curr_i - 1][curr_j - 1] + matrix.score(char1, char2)
            if abs(F[curr_i][curr_j] - diag_score) < 1e-9:
                aligned1.append(char1)
                aligned2.append(char2)
                curr_i -= 1
                curr_j -= 1
                continue

        # Check Up (Gap in seq2 / deletion of seq1 residue)
        if curr_i > 0:
            up_score = F[curr_i - 1][curr_j] - gap_penalty
            if abs(F[curr_i][curr_j] - up_score) < 1e-9:
                aligned1.append(seq1[curr_i - 1])
                aligned2.append("-")
                curr_i -= 1
                continue

        # Check Left (Gap in seq1 / insertion of seq2 residue)
        aligned1.append("-")
        aligned2.append(seq2[curr_j - 1])
        curr_j -= 1

    aligned1.reverse()
    aligned2.reverse()

    return AlignmentResult(
        seq1_aligned="".join(aligned1),
        seq2_aligned="".join(aligned2),
        score=F[m][n],
        start1=0,
        end1=m,
        start2=0,
        end2=n,
        alignment_type="global",
        gap_open=gap_penalty,
        gap_extend=gap_penalty,
        matrix_name=matrix.name,
        seq1=seq1,
        seq2=seq2,
    )
