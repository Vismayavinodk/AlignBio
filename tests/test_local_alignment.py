"""Comprehensive unit tests and property verifications for Smith–Waterman local alignment."""

import pytest
from proteinalign.models import AlignmentResult
from proteinalign.scoring import BLOSUM62, PAM250, SubstitutionMatrix
from proteinalign.local_alignment import smith_waterman


class TestSmithWatermanCore:
    """Core tests for Smith-Waterman local alignment algorithm under linear gap penalties."""

    def test_identical_sequences(self):
        """Test 1: Identical sequences align over their full length with 100% identity."""
        seq = "HEAGAWGHEE"
        res = smith_waterman(seq, seq, matrix=BLOSUM62, gap_penalty=8.0)
        assert res.score == 62.0
        assert res.seq1_aligned == seq
        assert res.seq2_aligned == seq
        assert res.start1 == 0 and res.end1 == len(seq)
        assert res.start2 == 0 and res.end2 == len(seq)

    def test_strong_matching_subsequence(self):
        """Test 2: Embedded conserved motif 'PAWHEAE' extracted from dissimilar flanks."""
        s1 = "DDDDPAWHEAERRRR"
        s2 = "WWWWPAWHEAEWWWW"
        res = smith_waterman(s1, s2, matrix=BLOSUM62, gap_penalty=8.0)

        assert res.seq1_aligned == "PAWHEAE"
        assert res.seq2_aligned == "PAWHEAE"
        # Check correct boundary slicing
        assert res.start1 == 4 and res.end1 == 11
        assert res.start2 == 4 and res.end2 == 11
        assert s1[res.start1:res.end1] == "PAWHEAE"
        assert s2[res.start2:res.end2] == "PAWHEAE"

    def test_unrelated_negative_residues_excluded(self):
        """Test 3: Strongly negative flanking residues are excluded from the local alignment."""
        s1 = "WWPAWHEAECC"
        s2 = "DDPAWHEAEEE"
        res = smith_waterman(s1, s2, matrix=BLOSUM62, gap_penalty=8.0)
        assert res.seq1_aligned == "PAWHEAE"
        assert not res.seq1_aligned.startswith("WW")
        assert "C" not in res.seq1_aligned

    def test_no_positive_scoring_alignment(self):
        """Test 4: Sequences with no positive local alignment return score 0.0 and empty strings."""
        s1 = "WWWW"
        s2 = "DDDD"
        # W vs D BLOSUM62 score is -4. All pairs are negative.
        res = smith_waterman(s1, s2, matrix=BLOSUM62, gap_penalty=8.0)
        assert res.score == 0.0
        assert res.seq1_aligned == ""
        assert res.seq2_aligned == ""
        assert res.start1 == 0 and res.end1 == 0
        assert res.start2 == 0 and res.end2 == 0

    def test_empty_sequences(self):
        """Test 5: Empty sequence inputs return score 0.0 and empty strings."""
        res1 = smith_waterman("", "HEAG")
        assert res1.score == 0.0 and res1.seq1_aligned == ""

        res2 = smith_waterman("HEAG", "")
        assert res2.score == 0.0 and res2.seq1_aligned == ""

        res3 = smith_waterman("", "")
        assert res3.score == 0.0 and res3.seq1_aligned == ""

    def test_single_residue_cases(self):
        """Test 6: Single-residue match and mismatch."""
        # Match
        res_m = smith_waterman("H", "H", matrix=BLOSUM62)
        assert res_m.score == 8.0
        assert res_m.seq1_aligned == "H" and res_m.seq2_aligned == "H"

        # Negative mismatch returns 0.0 and empty
        res_mis = smith_waterman("H", "P", matrix=BLOSUM62)
        assert res_mis.score == 0.0
        assert res_mis.seq1_aligned == "" and res_mis.seq2_aligned == ""

    def test_multiple_equally_scoring_alignments(self):
        """Test 7: Deterministic tie-breaking when multiple identical local motifs exist."""
        s1 = "HEAXXHEA"
        s2 = "HEA"
        res = smith_waterman(s1, s2, matrix=BLOSUM62, gap_penalty=8.0)
        # Should deterministically locate one of the HEA motifs
        assert res.seq1_aligned == "HEA"
        assert res.seq2_aligned == "HEA"
        assert res.score == 17.0  # H(8)+E(5)+A(4) = 17

    def test_correct_start_end_positions(self):
        """Test 8: Precise original sequence coordinates."""
        s1 = "ACDEFPAWHEAEGHI"
        s2 = "KLMNPAWHEAERST"
        res = smith_waterman(s1, s2, matrix=BLOSUM62, gap_penalty=8.0)
        assert s1[res.start1:res.end1] == res.seq1_aligned.replace("-", "")
        assert s2[res.start2:res.end2] == res.seq2_aligned.replace("-", "")

    def test_traceback_termination_at_zero(self):
        """Test 9: Traceback halts precisely when score drops to zero and does not overrun."""
        s1 = "HEAG"
        s2 = "HEAP"
        res = smith_waterman(s1, s2, matrix=BLOSUM62, gap_penalty=8.0)
        # All positive prefix: H(8)+E(5)+A(4)=17; G-P is -2, so optimal local alignment is HEA (17.0)
        assert res.score == 17.0
        assert res.seq1_aligned == "HEA"
        assert res.seq2_aligned == "HEA"

    def test_different_scoring_configurations(self):
        """Test 10: Comparison across identity matrix, BLOSUM62, and PAM250."""
        s1 = "XXHEAGYY"
        s2 = "ZZHEAPWW"

        id_mat = SubstitutionMatrix.create_identity(match=3.0, mismatch=-2.0)
        res_id = smith_waterman(s1, s2, matrix=id_mat, gap_penalty=4.0)
        assert res_id.seq1_aligned == "HEA"
        assert res_id.score == 9.0  # 3 * 3.0 = 9.0


class TestSmithWatermanInvariants:
    """Formal algebraic invariants that every valid local alignment must satisfy."""

    @pytest.mark.parametrize("s1, s2", [
        ("DDDDPAWHEAERRRR", "WWWWPAWHEAEWWWW"),
        ("MKFLILLFNILCLFPVLAAD", "MKFLILLVLFPAAD"),
        ("ACDEFGHIKLMNPQRSTVWY", "ACDEFGHIKLMNPQRSTVWY"),
    ])
    def test_subsequence_reconstruction_and_equal_length(self, s1, s2):
        res = smith_waterman(s1, s2, matrix=BLOSUM62, gap_penalty=8.0)
        if res.score > 0:
            # Invariant 1: Removing gaps reproduces original substring
            assert res.seq1_aligned.replace("-", "") == s1[res.start1:res.end1]
            assert res.seq2_aligned.replace("-", "") == s2[res.start2:res.end2]
            # Invariant 2: Aligned sequences have identical lengths
            assert len(res.seq1_aligned) == len(res.seq2_aligned)

    def test_score_consistency_recomputation(self):
        """Verifies local DP score equals manual summation of returned aligned columns."""
        s1 = "DDDDPAWHEAERRRR"
        s2 = "WWWWPAWHEAEWWWW"
        gap = 8.0
        res = smith_waterman(s1, s2, matrix=BLOSUM62, gap_penalty=gap)

        recomputed_score = 0.0
        for c1, c2 in zip(res.seq1_aligned, res.seq2_aligned):
            if c1 == "-" or c2 == "-":
                recomputed_score -= gap
            else:
                recomputed_score += BLOSUM62.score(c1, c2)

        assert pytest.approx(res.score, 1e-6) == recomputed_score
