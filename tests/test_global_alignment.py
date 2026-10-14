"""Comprehensive unit tests and property verifications for Needleman–Wunsch global alignment."""

import pytest
from proteinalign.models import AlignmentResult
from proteinalign.scoring import BLOSUM62, PAM250, SubstitutionMatrix
from proteinalign.global_alignment import needleman_wunsch


class TestNeedlemanWunschCore:
    """Core tests for Needleman-Wunsch algorithm under linear gap penalties."""

    def test_identical_sequences(self):
        """Test 1: Identical sequences align with zero gaps and score equal to sum of diagonal substitution scores."""
        seq = "HEAGAWGHEE"
        res = needleman_wunsch(seq, seq, matrix=BLOSUM62, gap_penalty=8.0)

        # Expected BLOSUM62 self-score: H(8)+E(5)+A(4)+G(6)+A(4)+W(11)+G(6)+H(8)+E(5)+E(5) = 62.0
        assert res.score == 62.0
        assert res.seq1_aligned == seq
        assert res.seq2_aligned == seq
        assert res.identity == 100.0
        assert res.start1 == 0 and res.end1 == len(seq)
        assert res.start2 == 0 and res.end2 == len(seq)
        assert res.seq1 == seq and res.seq2 == seq

    def test_single_mismatch(self):
        """Test 2: Sequences with a single mismatch at the end align without gaps when gap penalty > mismatch cost."""
        s1 = "HEAG"
        s2 = "HEAP"
        # Match H(8)+E(5)+A(4) = 17; Mismatch G-P = -2 => Total = 15.0
        res = needleman_wunsch(s1, s2, matrix=BLOSUM62, gap_penalty=8.0)
        assert res.score == 15.0
        assert res.seq1_aligned == "HEAG"
        assert res.seq2_aligned == "HEAP"

    def test_internal_gap_insertion(self):
        """Test 3: Internal gap insertion when residue is deleted."""
        # s1 has 'G' in the middle that s2 lacks
        s1 = "HEAGW"
        s2 = "HEAW"
        # H(8)+E(5)+A(4)+W(11) = 28; Gap penalty = -8 => Total = 20.0
        res = needleman_wunsch(s1, s2, matrix=BLOSUM62, gap_penalty=8.0)
        assert res.score == 20.0
        assert res.seq1_aligned == "HEAGW"
        assert res.seq2_aligned == "HEA-W"

    def test_different_lengths(self):
        """Test 4: Sequences with substantial length differences align correctly."""
        s1 = "HEAGAWGHEE"
        s2 = "HEAW"
        res = needleman_wunsch(s1, s2, matrix=BLOSUM62, gap_penalty=8.0)
        assert len(res.seq1_aligned) == len(res.seq2_aligned)
        assert res.seq1_aligned.replace("-", "") == s1
        assert res.seq2_aligned.replace("-", "") == s2

    def test_empty_versus_non_empty(self):
        """Test 5: Empty sequence against non-empty sequence produces all-gap row."""
        s = "HEA"
        res1 = needleman_wunsch(s, "", gap_penalty=8.0)
        assert res1.seq1_aligned == "HEA"
        assert res1.seq2_aligned == "---"
        assert res1.score == -24.0

        res2 = needleman_wunsch("", s, gap_penalty=8.0)
        assert res2.seq1_aligned == "---"
        assert res2.seq2_aligned == "HEA"
        assert res2.score == -24.0

    def test_both_empty(self):
        """Test 6: Two empty sequences return score 0 and empty strings."""
        res = needleman_wunsch("", "")
        assert res.score == 0.0
        assert res.seq1_aligned == ""
        assert res.seq2_aligned == ""

    def test_single_residue_cases(self):
        """Test 7: Single-residue match, mismatch, and deletion."""
        # Match
        res_match = needleman_wunsch("A", "A", matrix=BLOSUM62)
        assert res_match.score == 4.0
        assert res_match.seq1_aligned == "A" and res_match.seq2_aligned == "A"

        # Mismatch
        res_mismatch = needleman_wunsch("A", "W", matrix=BLOSUM62)
        assert res_mismatch.score == -3.0  # BLOSUM62 A-W = -3
        assert res_mismatch.seq1_aligned == "A" and res_mismatch.seq2_aligned == "W"

        # Deletion
        res_del = needleman_wunsch("A", "", gap_penalty=5.0)
        assert res_del.score == -5.0
        assert res_del.seq1_aligned == "A" and res_del.seq2_aligned == "-"

    def test_repeated_residues_and_homopolymers(self):
        """Test 8: Homopolymer repeats."""
        res = needleman_wunsch("AAAA", "AA", matrix=BLOSUM62, gap_penalty=8.0)
        assert res.seq1_aligned.replace("-", "") == "AAAA"
        assert res.seq2_aligned.replace("-", "") == "AA"
        assert len(res.seq1_aligned) == len(res.seq2_aligned)

    def test_configurable_gap_penalties(self):
        """Test 9: Varying gap penalty shifts preference from gaps to mismatches."""
        s1 = "HEAG"
        s2 = "HEAP"
        # Gap penalty = 1.0: cheaper to insert 2 gaps (score: 17 - 2 = 15) vs mismatch (15.0)
        res_low_gap = needleman_wunsch(s1, s2, matrix=BLOSUM62, gap_penalty=1.0)
        assert res_low_gap.score == 15.0  # Tie-breaking prefers diagonal match

        # High gap penalty = 20.0
        res_high_gap = needleman_wunsch(s1, s2, matrix=BLOSUM62, gap_penalty=20.0)
        assert res_high_gap.score == 15.0
        assert "-" not in res_high_gap.seq1_aligned

    def test_configurable_substitution_scoring(self):
        """Test 10: Comparison across BLOSUM62, PAM250, and simple identity matrix."""
        s1 = "HEAG"
        s2 = "HEAP"

        # 1. Simple identity matrix: match = +2, mismatch = -1
        identity_mat = SubstitutionMatrix.create_identity(match=2.0, mismatch=-1.0)
        res_id = needleman_wunsch(s1, s2, matrix=identity_mat, gap_penalty=3.0)
        # H-H (+2) + E-E (+2) + A-A (+2) + G-P (-1) = 5.0
        assert res_id.score == 5.0

        # 2. PAM250 matrix: H-H (+6) + E-E (+4) + A-A (+2) + G-P (0) = 12.0
        res_pam = needleman_wunsch(s1, s2, matrix=PAM250, gap_penalty=8.0)
        assert res_pam.score == 12.0


class TestNeedlemanWunschInvariants:
    """Formal algebraic invariants that every valid global alignment must satisfy."""

    @pytest.mark.parametrize("s1, s2", [
        ("MKFLILLFNILCLFPVLAAD", "MKFLILLVLFPAAD"),
        ("ACDEFGHIKLMNPQRSTVWY", "ACDEFGHIKLMNPQRSTVWY"),
        ("HEAGAWGHEE", "PAWHEAE"),
        ("W", "CDEF"),
    ])
    def test_reconstruction_and_equal_length(self, s1, s2):
        res = needleman_wunsch(s1, s2, matrix=BLOSUM62, gap_penalty=8.0)
        # Invariant 1: Degap reproduces inputs
        assert res.seq1_aligned.replace("-", "") == s1
        assert res.seq2_aligned.replace("-", "") == s2
        # Invariant 2: Aligned sequences have identical lengths
        assert len(res.seq1_aligned) == len(res.seq2_aligned)
        # Invariant 3: Length constraint max(m, n) <= L <= m + n
        assert max(len(s1), len(s2)) <= res.length <= len(s1) + len(s2)

    def test_score_consistency_recomputation(self):
        """Verifies that the DP score equals the manual summation across returned aligned columns."""
        s1 = "PAWHEAE"
        s2 = "HEAGAWGHEE"
        gap = 8.0
        res = needleman_wunsch(s1, s2, matrix=BLOSUM62, gap_penalty=gap)

        recomputed_score = 0.0
        for c1, c2 in zip(res.seq1_aligned, res.seq2_aligned):
            if c1 == "-" or c2 == "-":
                recomputed_score -= gap
            else:
                recomputed_score += BLOSUM62.score(c1, c2)

        assert pytest.approx(res.score, 1e-6) == recomputed_score
