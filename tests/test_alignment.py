"""Unit and Property-Based Tests for Sequence Alignment Algorithms.

Tests manually verifiable protein sequences and structural invariants:
1. Identity Property: Self-alignment has 100% PID and equals sum of self-substitution scores.
2. Symmetry Property: align(A, B).score == align(B, A).score for symmetric matrices.
3. Degap Invariant: Removing '-' from aligned outputs must reconstruct original sequences.
4. Equal Length Invariant: Both aligned sequence strings must have identical length.
5. Manually Verified Score: Small toy alignments whose score can be calculated by hand.
6. Local Isolation: Smith-Waterman correctly extracts embedded motif ignoring unaligned flanks.
7. Affine Gap Penalty Invariant: Affine penalties favor contiguous blocks over scattered indels.
8. Benchmark Reference Integrity: Real SABRE file check and Q-score self-comparison.
"""

from pathlib import Path
import pytest
from proteinalign.models import AlignmentResult
from proteinalign.scoring import BLOSUM62, SubstitutionMatrix
from proteinalign.global_align import needleman_wunsch
from proteinalign.local_align import smith_waterman
from proteinalign.affine_align import needleman_wunsch_affine, smith_waterman_affine
from proteinalign.io import read_msa
from proteinalign.evaluation import compute_pid, compute_q_score


# --- Test Sequences ---
SEQ_HEAGAWGHEE = "HEAGAWGHEE"
# Hand calculation for HEAGAWGHEE self-scores under BLOSUM62:
# H: 8, E: 5, A: 4, G: 6, A: 4, W: 11, G: 6, H: 8, E: 5, E: 5 => Total = 62.0
EXPECTED_SELF_SCORE_HEAGAWGHEE = 62.0

SEQ_SHORT_1 = "HEAG"
SEQ_SHORT_2 = "HEAP"
# BLOSUM62: H-H (+8), E-E (+5), A-A (+4), G-P (-2) => Total = 15.0
EXPECTED_HEAG_HEAP_SCORE = 15.0


class TestNeedlemanWunsch:
    """Tests for Global Alignment."""

    def test_identity_self_alignment(self):
        result = needleman_wunsch(SEQ_HEAGAWGHEE, SEQ_HEAGAWGHEE, matrix=BLOSUM62, gap_penalty=8.0)
        assert result.seq1_aligned == SEQ_HEAGAWGHEE
        assert result.seq2_aligned == SEQ_HEAGAWGHEE
        assert result.score == EXPECTED_SELF_SCORE_HEAGAWGHEE
        assert result.identity == 100.0

    def test_symmetry(self):
        res_ab = needleman_wunsch("PAWHEAE", "HEAGAWGHEE", matrix=BLOSUM62, gap_penalty=8.0)
        res_ba = needleman_wunsch("HEAGAWGHEE", "PAWHEAE", matrix=BLOSUM62, gap_penalty=8.0)
        assert res_ab.score == res_ba.score

    def test_manually_verifiable_mismatch(self):
        result = needleman_wunsch(SEQ_SHORT_1, SEQ_SHORT_2, matrix=BLOSUM62, gap_penalty=8.0)
        assert result.score == EXPECTED_HEAG_HEAP_SCORE
        assert result.seq1_aligned == "HEAG"
        assert result.seq2_aligned == "HEAP"

    def test_single_deletion(self):
        # seq1 has extra 'G' at end
        result = needleman_wunsch("HEA", "HEAG", matrix=BLOSUM62, gap_penalty=8.0)
        # H-H (+8) + E-E (+5) + A-A (+4) - Gap (8) = 9.0
        assert result.score == 9.0
        assert result.seq1_aligned == "HEA-"
        assert result.seq2_aligned == "HEAG"

    def test_reconstruction_invariant(self):
        s1 = "MKFLILLFNILCLFPVLAAD"
        s2 = "MKFLILLVLFPAAD"
        result = needleman_wunsch(s1, s2, matrix=BLOSUM62, gap_penalty=8.0)
        assert result.seq1_aligned.replace("-", "") == s1
        assert result.seq2_aligned.replace("-", "") == s2
        assert len(result.seq1_aligned) == len(result.seq2_aligned)


class TestSmithWaterman:
    """Tests for Local Alignment."""

    def test_local_motif_isolation(self):
        # Motif 'PAWHEAE' embedded inside dissimilar flanking residues (D vs W: -4, R vs E: 0)
        s1 = "DDDDPAWHEAERRRR"
        s2 = "WWWWPAWHEAEWWWW"
        result = smith_waterman(s1, s2, matrix=BLOSUM62, gap_penalty=8.0)

        assert result.seq1_aligned == "PAWHEAE"
        assert result.seq2_aligned == "PAWHEAE"
        assert result.start1 == 4
        assert result.end1 == 11
        assert result.start2 == 4
        assert result.end2 == 11

    def test_local_self_alignment(self):
        result = smith_waterman(SEQ_HEAGAWGHEE, SEQ_HEAGAWGHEE, matrix=BLOSUM62, gap_penalty=8.0)
        assert result.score == EXPECTED_SELF_SCORE_HEAGAWGHEE
        assert result.seq1_aligned == SEQ_HEAGAWGHEE
        assert result.seq2_aligned == SEQ_HEAGAWGHEE


class TestAffineGapPenalties:
    """Tests for Gotoh 3-State Affine Gap DP."""

    def test_affine_identity(self):
        result = needleman_wunsch_affine(SEQ_HEAGAWGHEE, SEQ_HEAGAWGHEE, matrix=BLOSUM62, gap_open=10.0, gap_extend=1.0)
        assert result.score == EXPECTED_SELF_SCORE_HEAGAWGHEE
        assert result.seq1_aligned == SEQ_HEAGAWGHEE
        assert result.seq2_aligned == SEQ_HEAGAWGHEE

    def test_affine_contiguous_gap_cost(self):
        # Gap of length 2 under gap_open=10, gap_extend=1:
        # Cost = -(10 + 1*1) = -11.
        # Under linear penalty of 10, cost would be -20.
        s1 = "HEAWG"
        s2 = "HEA"
        # Match HEA: H(8)+E(5)+A(4) = 17.
        # Gap cost: -11. Expected score = 17 - 11 = 6.0
        result = needleman_wunsch_affine(s1, s2, matrix=BLOSUM62, gap_open=10.0, gap_extend=1.0)
        assert result.score == 6.0
        assert result.seq1_aligned == "HEAWG"
        assert result.seq2_aligned == "HEA--"

    def test_affine_local_motif(self):
        # Flanked by acidic/hydrophobic pairs with negative BLOSUM scores
        s1 = "DDPAWHEAERR"
        s2 = "WWPAWHEAEEE"
        result = smith_waterman_affine(s1, s2, matrix=BLOSUM62, gap_open=10.0, gap_extend=1.0)
        assert result.seq1_aligned == "PAWHEAE"
        assert result.seq2_aligned == "PAWHEAE"


class TestEvaluationAndBenchmark:
    """Tests for PID and Q-score metrics."""

    def test_pid_calculation(self):
        assert compute_pid("HEAG", "HEAP") == 75.0
        assert compute_pid("HEAG", "HEAG") == 100.0
        assert compute_pid("HEA-", "HEAG") == 75.0

    def test_sabrem_reference_file_q_score(self):
        ref_path = Path("data/raw/bench1.0/sabrem/ref/sup_003")
        if ref_path.exists():
            msa = read_msa(ref_path)
            assert len(msa) > 0
            # Q-score of reference against itself must be 1.0 (100% accuracy)
            q = compute_q_score(msa, msa, uppercase_only=True)
            assert pytest.approx(q, 1e-6) == 1.0
