"""Position-Specific Scoring Matrix (PSSM) construction and profile-sequence alignment."""

import math
from typing import Dict, List, Optional, Sequence
from .models import AlignmentResult
from .scoring import AMINO_ACIDS
from .profiles import (
    ROBINSON_FREQUENCIES,
    extract_alignment_columns,
    compute_column_counts,
    compute_profile_frequencies,
)


class PSSM:
    """Position-Specific Scoring Matrix representing a sequence profile."""

    def __init__(
        self,
        matrix: List[Dict[str, float]],
        name: str = "pssm",
        consensus: str = "",
        background: Optional[Dict[str, float]] = None,
    ) -> None:
        self.matrix = matrix
        self.name = name
        self.consensus = consensus
        self.background = background or ROBINSON_FREQUENCIES

    def __len__(self) -> int:
        return len(self.matrix)

    def score(self, pos: int, residue: str) -> float:
        """Returns the log-odds score for residue at 0-indexed profile column pos."""
        if 0 <= pos < len(self.matrix):
            return self.matrix[pos].get(residue.upper(), -4.0)
        return -4.0

    @classmethod
    def from_msa(
        cls,
        aligned_sequences: Sequence[str],
        name: str = "profile",
        pseudocount_weight: float = 1.0,
        background: Optional[Dict[str, float]] = None,
        scale: float = 2.0,
    ) -> "PSSM":
        """Constructs a PSSM from aligned multiple sequences with pseudocount regularization.

        Score P(i, a) = round(scale * log2(f_i(a) / p_a))
        """
        if background is None:
            background = ROBINSON_FREQUENCIES

        cols = extract_alignment_columns(aligned_sequences)
        counts = compute_column_counts(cols)
        frequencies = compute_profile_frequencies(counts, pseudocount_weight=pseudocount_weight, background=background)

        matrix = []
        consensus_chars = []

        for f_dict in frequencies:
            col_scores = {}
            best_aa = "A"
            best_freq = -1.0

            for aa in AMINO_ACIDS:
                prob = f_dict[aa]
                bg = background.get(aa, 0.05)
                # Compute log-odds ratio in half-bits / scaled units
                ratio = max(prob / bg, 1e-6)
                log_odds = scale * math.log2(ratio)
                col_scores[aa] = round(log_odds, 2)

                if prob > best_freq:
                    best_freq = prob
                    best_aa = aa

            matrix.append(col_scores)
            consensus_chars.append(best_aa)

        return cls(
            matrix=matrix,
            name=name,
            consensus="".join(consensus_chars),
            background=background,
        )

    def align_sequence_global(self, query: str, gap_penalty: float = 8.0) -> AlignmentResult:
        """Aligns a query sequence to this PSSM profile using global Needleman-Wunsch."""
        m = len(self)
        n = len(query)

        F = [[0.0] * (n + 1) for _ in range(m + 1)]
        for i in range(1, m + 1):
            F[i][0] = -i * gap_penalty
        for j in range(1, n + 1):
            F[0][j] = -j * gap_penalty

        for i in range(1, m + 1):
            for j in range(1, n + 1):
                char_q = query[j - 1]
                diag = F[i - 1][j - 1] + self.score(i - 1, char_q)
                up = F[i - 1][j] - gap_penalty
                left = F[i][j - 1] - gap_penalty
                F[i][j] = max(diag, up, left)

        # Traceback
        prof_aln, seq_aln = [], []
        curr_i, curr_j = m, n

        while curr_i > 0 or curr_j > 0:
            if curr_i > 0 and curr_j > 0:
                char_q = query[curr_j - 1]
                diag_score = F[curr_i - 1][curr_j - 1] + self.score(curr_i - 1, char_q)
                if abs(F[curr_i][curr_j] - diag_score) < 1e-9:
                    prof_aln.append(self.consensus[curr_i - 1] if self.consensus else f"P{curr_i}")
                    seq_aln.append(char_q)
                    curr_i -= 1
                    curr_j -= 1
                    continue

            if curr_i > 0:
                up_score = F[curr_i - 1][curr_j] - gap_penalty
                if abs(F[curr_i][curr_j] - up_score) < 1e-9:
                    prof_aln.append(self.consensus[curr_i - 1] if self.consensus else f"P{curr_i}")
                    seq_aln.append("-")
                    curr_i -= 1
                    continue

            prof_aln.append("-")
            seq_aln.append(query[curr_j - 1])
            curr_j -= 1

        prof_aln.reverse()
        seq_aln.reverse()

        return AlignmentResult(
            seq1_aligned="".join(prof_aln),
            seq2_aligned="".join(seq_aln),
            score=F[m][n],
            start1=0,
            end1=m,
            start2=0,
            end2=n,
            alignment_type="pssm_global",
            gap_open=gap_penalty,
            gap_extend=gap_penalty,
            matrix_name=self.name,
        )
