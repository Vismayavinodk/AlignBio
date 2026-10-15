"""Position-Specific Scoring Matrix (PSSM) construction and profile-sequence alignment.

This module implements:
1. PSSM generation from position-wise frequency profiles or aligned multiple sequence alignments.
2. Position-specific log-odds scoring in configurable units (e.g. half-bits, third-bits, nats).
3. Dynamic programming profile-to-sequence alignment for global (Needleman-Wunsch) and
   local (Smith-Waterman) modes with linear and Gotoh (1982) affine gap penalties.
4. Independent score recomputation verification directly from reconstructed alignments.
5. Research-grade data leakage auditing against profile training sequence IDs.

Mathematical Recurrence:
    For profile column i (0 <= i < L) and amino acid a:
        P_{i, a} = (c_{i, a} + alpha * q_a) / (N_i + alpha)
        S_{i, a} = lambda * log_b ( P_{i, a} / q_a )
    where:
        - alpha is the pseudocount strength.
        - q_a is the background amino acid probability.
        - b is the logarithm base (default 2.0).
        - lambda is the scaling factor (default 2.0, giving half-bits).
"""

import math
from typing import Dict, List, Optional, Sequence, Union
from .models import AlignmentResult
from .scoring import AMINO_ACIDS, validate_sequence
from .profiles import (
    ROBINSON_FREQUENCIES,
    Profile,
    build_profile_from_msa,
    validate_background_distribution,
    validate_no_data_leakage,
)


NEG_INF = -1e9


class PSSM:
    """Position-Specific Scoring Matrix representing an evolutionary profile."""

    def __init__(
        self,
        matrix: List[Dict[str, float]],
        name: str = "pssm",
        consensus: str = "",
        background: Optional[Dict[str, float]] = None,
        scale: float = 2.0,
        log_base: float = 2.0,
        profile: Optional[Profile] = None,
        sequence_ids: Optional[List[str]] = None,
    ) -> None:
        """Initializes a PSSM instance.

        Args:
            matrix: List of L dictionaries mapping amino acid characters to log-odds scores.
            name: Identifier for this scoring profile.
            consensus: Consensus sequence string across profile columns.
            background: Background amino acid distribution.
            scale: Scaling factor lambda (multiplier for log-odds scores).
            log_base: Base of the logarithm used (e.g. 2.0 or math.e).
            profile: Optional underlying Profile instance.
            sequence_ids: List of training sequence IDs used to build this matrix.
        """
        self.matrix = matrix
        self.name = name
        self.consensus = consensus
        self.background = background or dict(ROBINSON_FREQUENCIES)
        self.scale = scale
        self.log_base = log_base
        self.profile = profile
        self.sequence_ids = sequence_ids or (profile.sequence_ids if profile else [])

    def __len__(self) -> int:
        return len(self.matrix)

    def score(self, pos: int, residue: str) -> float:
        """Returns the position-specific log-odds score for a residue at 0-indexed column pos.

        Args:
            pos: 0-indexed profile column (0 <= pos < len(self)).
            residue: Single-character amino acid code.

        Returns:
            Float log-odds score.

        Raises:
            IndexError: If pos is out of range [0, len(self)).
            ValueError: If residue is invalid.
        """
        if not (0 <= pos < len(self.matrix)):
            raise IndexError(f"Profile position {pos} out of range [0, {len(self.matrix)}).")

        u = residue.upper()
        col = self.matrix[pos]

        if u in col:
            return col[u]

        # Handle ambiguous residues according to documented expectation
        if u == "B":  # Asx: 50% Asn + 50% Asp
            return 0.5 * col.get("N", 0.0) + 0.5 * col.get("D", 0.0)
        elif u == "Z":  # Glx: 50% Gln + 50% Glu
            return 0.5 * col.get("Q", 0.0) + 0.5 * col.get("E", 0.0)
        elif u == "X":  # Unknown: Expected log-odds under background distribution
            return sum(self.background.get(aa, 0.05) * col.get(aa, 0.0) for aa in AMINO_ACIDS)

        raise ValueError(f"Unknown amino acid character '{residue}' at position {pos}.")

    def expected_score(self, pos: int) -> float:
        """Computes the expected score at column pos under the background distribution."""
        if not (0 <= pos < len(self.matrix)):
            raise IndexError(f"Profile position {pos} out of range [0, {len(self.matrix)}).")
        return sum(self.background.get(aa, 0.05) * self.matrix[pos].get(aa, 0.0) for aa in AMINO_ACIDS)

    @classmethod
    def from_profile(
        cls,
        profile: Profile,
        scale: float = 2.0,
        log_base: float = 2.0,
    ) -> "PSSM":
        """Constructs a PSSM directly from an existing Profile instance.

        Score:
            S_{i, a} = lambda * log_b ( P_{i, a} / q_a )

        Args:
            profile: Profile instance with regularized frequencies and background distribution.
            scale: Multiplier lambda (e.g., 2.0 for half-bits with log_base 2.0).
            log_base: Logarithm base (default 2.0).

        Returns:
            PSSM instance.

        Raises:
            ValueError: If scale <= 0 or log_base <= 1.0.
        """
        if scale <= 0:
            raise ValueError(f"PSSM scale parameter lambda must be positive (received {scale}).")
        if log_base <= 1.0:
            raise ValueError(f"Logarithm base must be > 1.0 (received {log_base}).")

        validate_background_distribution(profile.background)

        matrix: List[Dict[str, float]] = []
        for col_idx in range(profile.alignment_length):
            f_dict = profile.get_column_frequencies(col_idx)
            col_scores: Dict[str, float] = {}

            for aa in AMINO_ACIDS:
                p_ia = f_dict[aa]
                q_a = profile.background[aa]
                # Ratio of observed target frequency to background probability
                ratio = max(p_ia / q_a, 1e-9)
                if log_base == 2.0:
                    val = scale * math.log2(ratio)
                elif abs(log_base - math.e) < 1e-9:
                    val = scale * math.log(ratio)
                else:
                    val = scale * (math.log(ratio) / math.log(log_base))

                col_scores[aa] = round(val, 3)

            matrix.append(col_scores)

        return cls(
            matrix=matrix,
            name=f"pssm_{profile.name}",
            consensus=profile.consensus,
            background=profile.background,
            scale=scale,
            log_base=log_base,
            profile=profile,
            sequence_ids=profile.sequence_ids,
        )

    @classmethod
    def from_msa(
        cls,
        aligned_sequences: Union[Dict[str, str], Sequence[str]],
        name: str = "profile",
        pseudocount_weight: float = 1.0,
        background: Optional[Dict[str, float]] = None,
        scale: float = 2.0,
        log_base: float = 2.0,
        sequence_ids: Optional[List[str]] = None,
    ) -> "PSSM":
        """Constructs a PSSM directly from multiple aligned sequences with pseudocounts.

        Args:
            aligned_sequences: Dict[header, seq] or list of equal-length aligned strings.
            name: Profile name.
            pseudocount_weight: Strength of background prior (alpha).
            background: Background amino acid distribution.
            scale: Score multiplier lambda.
            log_base: Logarithm base.
            sequence_ids: Optional list of sequence headers.

        Returns:
            PSSM instance.
        """
        prof = build_profile_from_msa(
            aligned_sequences,
            name=name,
            pseudocount_weight=pseudocount_weight,
            background=background,
            sequence_ids=sequence_ids,
        )
        return cls.from_profile(prof, scale=scale, log_base=log_base)

    # -------------------------------------------------------------
    # Profile-to-Sequence Dynamic Programming Alignment
    # -------------------------------------------------------------

    def align_global(
        self,
        query: str,
        gap_penalty: float = 8.0,
        gap_open: Optional[float] = None,
        gap_extend: Optional[float] = None,
        allow_ambiguous: bool = True,
        query_id: Optional[str] = None,
        allow_leakage: bool = False,
    ) -> AlignmentResult:
        """Aligns a query sequence globally against this PSSM profile using dynamic programming.

        Supports both linear and Gotoh 3-state affine gap penalties.

        Args:
            query: Target protein sequence string.
            gap_penalty: Linear gap penalty magnitude (used when gap_open is None).
            gap_open: Affine gap-opening penalty magnitude (d > 0).
            gap_extend: Affine gap-extension penalty magnitude (e > 0).
            allow_ambiguous: If True, permits ambiguous IUPAC amino acid codes.
            query_id: Optional ID of the query sequence (checked for data leakage).
            allow_leakage: If True, bypasses data leakage validation (demonstration only).

        Returns:
            AlignmentResult with profile consensus as seq1 and query as seq2.
        """
        validate_sequence(query, allow_ambiguous=allow_ambiguous)
        if self.profile:
            validate_no_data_leakage(self.profile, query_id=query_id, allow_leakage=allow_leakage)

        is_affine = gap_open is not None and gap_extend is not None
        if is_affine:
            if gap_open < 0 or gap_extend < 0:
                raise ValueError("Gap penalties must be non-negative.")
            if gap_open < gap_extend:
                raise ValueError(f"gap_open ({gap_open}) must be >= gap_extend ({gap_extend}).")
            g_open, g_ext = float(gap_open), float(gap_extend)
        else:
            if gap_penalty < 0:
                raise ValueError("Gap penalty must be non-negative.")
            g_open, g_ext = float(gap_penalty), float(gap_penalty)

        m = len(self)
        n = len(query)

        # Empty cases
        if m == 0 and n == 0:
            return AlignmentResult("", "", 0.0, 0, 0, 0, 0, "pssm_global", g_open, g_ext, self.name, "", "")
        if m == 0:
            cost = -(g_open + (n - 1) * g_ext) if is_affine else -n * g_open
            return AlignmentResult("-" * n, query, cost, 0, 0, 0, n, "pssm_global", g_open, g_ext, self.name, "", query)
        if n == 0:
            cost = -(g_open + (m - 1) * g_ext) if is_affine else -m * g_open
            return AlignmentResult(self.consensus, "-" * m, cost, 0, m, 0, 0, "pssm_global", g_open, g_ext, self.name, self.consensus, "")

        if is_affine:
            return self._align_global_affine(query, g_open, g_ext)
        return self._align_global_linear(query, g_open)

    def _align_global_linear(self, query: str, d: float) -> AlignmentResult:
        m, n = len(self), len(query)
        F = [[0.0] * (n + 1) for _ in range(m + 1)]

        for i in range(1, m + 1):
            F[i][0] = -i * d
        for j in range(1, n + 1):
            F[0][j] = -j * d

        for i in range(1, m + 1):
            for j in range(1, n + 1):
                char_q = query[j - 1]
                diag = F[i - 1][j - 1] + self.score(i - 1, char_q)
                up = F[i - 1][j] - d
                left = F[i][j - 1] - d
                F[i][j] = max(diag, up, left)

        # Traceback with deterministic priority: Diag > Up > Left
        prof_aln, seq_aln = [], []
        curr_i, curr_j = m, n

        while curr_i > 0 or curr_j > 0:
            if curr_i == 0:
                prof_aln.append("-")
                seq_aln.append(query[curr_j - 1])
                curr_j -= 1
                continue
            if curr_j == 0:
                prof_aln.append(self.consensus[curr_i - 1])
                seq_aln.append("-")
                curr_i -= 1
                continue

            char_q = query[curr_j - 1]
            diag_score = F[curr_i - 1][curr_j - 1] + self.score(curr_i - 1, char_q)
            if abs(F[curr_i][curr_j] - diag_score) < 1e-9:
                prof_aln.append(self.consensus[curr_i - 1])
                seq_aln.append(char_q)
                curr_i -= 1
                curr_j -= 1
                continue

            up_score = F[curr_i - 1][curr_j] - d
            if abs(F[curr_i][curr_j] - up_score) < 1e-9:
                prof_aln.append(self.consensus[curr_i - 1])
                seq_aln.append("-")
                curr_i -= 1
                continue

            prof_aln.append("-")
            seq_aln.append(char_q)
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
            gap_open=d,
            gap_extend=d,
            matrix_name=self.name,
            seq1=self.consensus,
            seq2=query,
        )

    def _align_global_affine(self, query: str, gap_open: float, gap_extend: float) -> AlignmentResult:
        m, n = len(self), len(query)
        M = [[NEG_INF] * (n + 1) for _ in range(m + 1)]
        Ix = [[NEG_INF] * (n + 1) for _ in range(m + 1)]
        Iy = [[NEG_INF] * (n + 1) for _ in range(m + 1)]

        M[0][0] = 0.0
        for i in range(1, m + 1):
            Ix[i][0] = -(gap_open + (i - 1) * gap_extend)
        for j in range(1, n + 1):
            Iy[0][j] = -(gap_open + (j - 1) * gap_extend)

        for i in range(1, m + 1):
            for j in range(1, n + 1):
                char_q = query[j - 1]
                sub_score = self.score(i - 1, char_q)

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

        prof_aln, seq_aln = [], []
        curr_i, curr_j = m, n

        while curr_i > 0 or curr_j > 0:
            if curr_i == 0:
                prof_aln.append("-")
                seq_aln.append(query[curr_j - 1])
                curr_j -= 1
                continue
            if curr_j == 0:
                prof_aln.append(self.consensus[curr_i - 1])
                seq_aln.append("-")
                curr_i -= 1
                continue

            if curr_state == 0:
                char_q = query[curr_j - 1]
                prof_aln.append(self.consensus[curr_i - 1])
                seq_aln.append(char_q)
                prev_val = M[curr_i][curr_j] - self.score(curr_i - 1, char_q)

                if abs(prev_val - M[curr_i - 1][curr_j - 1]) < 1e-9:
                    curr_state = 0
                elif abs(prev_val - Ix[curr_i - 1][curr_j - 1]) < 1e-9:
                    curr_state = 1
                else:
                    curr_state = 2

                curr_i -= 1
                curr_j -= 1

            elif curr_state == 1:
                prof_aln.append(self.consensus[curr_i - 1])
                seq_aln.append("-")
                if abs(Ix[curr_i][curr_j] - (M[curr_i - 1][curr_j] - gap_open)) < 1e-9:
                    curr_state = 0
                else:
                    curr_state = 1
                curr_i -= 1

            elif curr_state == 2:
                prof_aln.append("-")
                seq_aln.append(query[curr_j - 1])
                if abs(Iy[curr_i][curr_j] - (M[curr_i][curr_j - 1] - gap_open)) < 1e-9:
                    curr_state = 0
                else:
                    curr_state = 2
                curr_j -= 1

        prof_aln.reverse()
        seq_aln.reverse()

        return AlignmentResult(
            seq1_aligned="".join(prof_aln),
            seq2_aligned="".join(seq_aln),
            score=best_score,
            start1=0,
            end1=m,
            start2=0,
            end2=n,
            alignment_type="pssm_global_affine",
            gap_open=gap_open,
            gap_extend=gap_extend,
            matrix_name=self.name,
            seq1=self.consensus,
            seq2=query,
        )

    def align_local(
        self,
        query: str,
        gap_penalty: float = 8.0,
        gap_open: Optional[float] = None,
        gap_extend: Optional[float] = None,
        allow_ambiguous: bool = True,
        query_id: Optional[str] = None,
        allow_leakage: bool = False,
    ) -> AlignmentResult:
        """Aligns a query sequence locally against this PSSM profile using Smith-Waterman DP.

        Supports both linear and Gotoh 3-state affine gap penalties.

        Args:
            query: Target protein sequence string.
            gap_penalty: Linear gap penalty magnitude.
            gap_open: Affine gap-opening penalty magnitude.
            gap_extend: Affine gap-extension penalty magnitude.
            allow_ambiguous: If True, permits ambiguous IUPAC codes.
            query_id: Optional ID of the query sequence.
            allow_leakage: If True, bypasses data leakage validation.

        Returns:
            AlignmentResult with local alignment coordinates and scores.
        """
        validate_sequence(query, allow_ambiguous=allow_ambiguous)
        if self.profile:
            validate_no_data_leakage(self.profile, query_id=query_id, allow_leakage=allow_leakage)

        is_affine = gap_open is not None and gap_extend is not None
        if is_affine:
            if gap_open < 0 or gap_extend < 0:
                raise ValueError("Gap penalties must be non-negative.")
            if gap_open < gap_extend:
                raise ValueError(f"gap_open ({gap_open}) must be >= gap_extend ({gap_extend}).")
            g_open, g_ext = float(gap_open), float(gap_extend)
        else:
            if gap_penalty < 0:
                raise ValueError("Gap penalty must be non-negative.")
            g_open, g_ext = float(gap_penalty), float(gap_penalty)

        m = len(self)
        n = len(query)

        if m == 0 or n == 0:
            return AlignmentResult("", "", 0.0, 0, 0, 0, 0, "pssm_local", g_open, g_ext, self.name, self.consensus, query)

        if is_affine:
            return self._align_local_affine(query, g_open, g_ext)
        return self._align_local_linear(query, g_open)

    def _align_local_linear(self, query: str, d: float) -> AlignmentResult:
        m, n = len(self), len(query)
        H = [[0.0] * (n + 1) for _ in range(m + 1)]

        max_score = 0.0
        max_i, max_j = 0, 0

        for i in range(1, m + 1):
            for j in range(1, n + 1):
                char_q = query[j - 1]
                diag = H[i - 1][j - 1] + self.score(i - 1, char_q)
                up = H[i - 1][j] - d
                left = H[i][j - 1] - d
                val = max(0.0, diag, up, left)
                H[i][j] = val

                if val > max_score:
                    max_score = val
                    max_i, max_j = i, j

        if max_score <= 0.0:
            return AlignmentResult("", "", 0.0, 0, 0, 0, 0, "pssm_local", d, d, self.name, self.consensus, query)

        prof_aln, seq_aln = [], []
        curr_i, curr_j = max_i, max_j

        while curr_i > 0 and curr_j > 0 and H[curr_i][curr_j] > 0.0:
            char_q = query[curr_j - 1]
            diag_score = H[curr_i - 1][curr_j - 1] + self.score(curr_i - 1, char_q)

            if abs(H[curr_i][curr_j] - diag_score) < 1e-9:
                prof_aln.append(self.consensus[curr_i - 1])
                seq_aln.append(char_q)
                curr_i -= 1
                curr_j -= 1
            elif abs(H[curr_i][curr_j] - (H[curr_i - 1][curr_j] - d)) < 1e-9:
                prof_aln.append(self.consensus[curr_i - 1])
                seq_aln.append("-")
                curr_i -= 1
            elif abs(H[curr_i][curr_j] - (H[curr_i][curr_j - 1] - d)) < 1e-9:
                prof_aln.append("-")
                seq_aln.append(char_q)
                curr_j -= 1
            else:
                break

        prof_aln.reverse()
        seq_aln.reverse()

        return AlignmentResult(
            seq1_aligned="".join(prof_aln),
            seq2_aligned="".join(seq_aln),
            score=max_score,
            start1=curr_i,
            end1=max_i,
            start2=curr_j,
            end2=max_j,
            alignment_type="pssm_local",
            gap_open=d,
            gap_extend=d,
            matrix_name=self.name,
            seq1=self.consensus,
            seq2=query,
        )

    def _align_local_affine(self, query: str, gap_open: float, gap_extend: float) -> AlignmentResult:
        m, n = len(self), len(query)
        M = [[0.0] * (n + 1) for _ in range(m + 1)]
        Ix = [[0.0] * (n + 1) for _ in range(m + 1)]
        Iy = [[0.0] * (n + 1) for _ in range(m + 1)]

        max_score = 0.0
        max_i, max_j = 0, 0

        for i in range(1, m + 1):
            for j in range(1, n + 1):
                char_q = query[j - 1]
                sub_score = self.score(i - 1, char_q)

                prev_max = max(0.0, M[i - 1][j - 1], Ix[i - 1][j - 1], Iy[i - 1][j - 1])
                m_val = max(0.0, sub_score + prev_max)
                M[i][j] = m_val

                Ix[i][j] = max(0.0, M[i - 1][j] - gap_open, Ix[i - 1][j] - gap_extend)
                Iy[i][j] = max(0.0, M[i][j - 1] - gap_open, Iy[i][j - 1] - gap_extend)

                if m_val > max_score:
                    max_score = m_val
                    max_i, max_j = i, j

        if max_score <= 0.0:
            return AlignmentResult("", "", 0.0, 0, 0, 0, 0, "pssm_local_affine", gap_open, gap_extend, self.name, self.consensus, query)

        prof_aln, seq_aln = [], []
        curr_i, curr_j = max_i, max_j
        curr_state = 0

        while curr_i > 0 and curr_j > 0:
            if curr_state == 0:
                if M[curr_i][curr_j] <= 0.0:
                    break
                char_q = query[curr_j - 1]
                prev_val = M[curr_i][curr_j] - self.score(curr_i - 1, char_q)

                prof_aln.append(self.consensus[curr_i - 1])
                seq_aln.append(char_q)

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
                prof_aln.append(self.consensus[curr_i - 1])
                seq_aln.append("-")
                if abs(Ix[curr_i][curr_j] - (M[curr_i - 1][curr_j] - gap_open)) < 1e-9:
                    curr_state = 0
                else:
                    curr_state = 1
                curr_i -= 1

            elif curr_state == 2:
                if Iy[curr_i][curr_j] <= 0.0:
                    break
                prof_aln.append("-")
                seq_aln.append(query[curr_j - 1])
                if abs(Iy[curr_i][curr_j] - (M[curr_i][curr_j - 1] - gap_open)) < 1e-9:
                    curr_state = 0
                else:
                    curr_state = 2
                curr_j -= 1

        prof_aln.reverse()
        seq_aln.reverse()

        return AlignmentResult(
            seq1_aligned="".join(prof_aln),
            seq2_aligned="".join(seq_aln),
            score=max_score,
            start1=curr_i,
            end1=max_i,
            start2=curr_j,
            end2=max_j,
            alignment_type="pssm_local_affine",
            gap_open=gap_open,
            gap_extend=gap_extend,
            matrix_name=self.name,
            seq1=self.consensus,
            seq2=query,
        )

    # -------------------------------------------------------------
    # Independent Score Verification
    # -------------------------------------------------------------

    def recompute_score(self, result: AlignmentResult) -> float:
        """Independently recomputes the alignment score by summing column scores and gap penalties.

        Args:
            result: AlignmentResult produced by profile alignment.

        Returns:
            Recomputed score float.
        """
        p_aln = result.seq1_aligned
        q_aln = result.seq2_aligned
        if len(p_aln) != len(q_aln):
            raise ValueError("Aligned sequence strings must have identical lengths.")

        if not p_aln:
            return 0.0

        is_affine = result.gap_open != result.gap_extend
        g_open = result.gap_open
        g_ext = result.gap_extend

        total_score = 0.0
        prof_idx = result.start1
        in_p_gap = False
        in_q_gap = False

        for char_p, char_q in zip(p_aln, q_aln):
            if char_p != "-" and char_q != "-":
                total_score += self.score(prof_idx, char_q)
                prof_idx += 1
                in_p_gap = False
                in_q_gap = False

            elif char_p != "-" and char_q == "-":
                # Deletion of profile position in query (gap in query)
                prof_idx += 1
                in_p_gap = False
                if is_affine:
                    total_score -= (g_ext if in_q_gap else g_open)
                else:
                    total_score -= g_open
                in_q_gap = True

            elif char_p == "-" and char_q != "-":
                # Insertion in query (gap in profile)
                in_q_gap = False
                if is_affine:
                    total_score -= (g_ext if in_p_gap else g_open)
                else:
                    total_score -= g_open
                in_p_gap = True

        return round(total_score, 4)
