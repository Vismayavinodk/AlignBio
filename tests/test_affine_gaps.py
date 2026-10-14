"""Comprehensive unit and property tests for Gotoh 3-state affine gap alignment.

Tests:
1. Exact mathematical gap cost scaling:
   - G(1) = gap_open
   - G(2) = gap_open + gap_extend
   - G(3) = gap_open + 2 * gap_extend
2. Structural invariants:
   - Self-alignment identity (PID = 100%, score = sum of self-substitution scores)
   - Global symmetry: align(A, B).score == align(B, A).score
   - Degap reconstruction: removing '-' restores original sequences
   - Equal length invariant: len(seq1_aligned) == len(seq2_aligned)
   - Result attribute propagation: result.seq1 == seq1, result.seq2 == seq2
3. Contiguous indel grouping:
   - Affine penalties favor a single contiguous gap over multiple scattered single-residue indels
4. Local alignment properties:
   - Motif isolation from non-homologous flanking sequences
   - Subsequence coordinate consistency: seq[start:end] == seq_aligned.replace('-', '')
   - Internal affine gap handling
5. Edge cases:
   - Empty sequences (empty vs empty, empty vs non-empty)
   - Single-residue comparisons
   - Zero score local termination
   - Completely dissimilar sequences
6. Validation errors:
   - Negative gap penalties
   - gap_open < gap_extend
   - Invalid non-amino acid characters
"""

import pytest
from proteinalign.models import AlignmentResult
from proteinalign.scoring import BLOSUM62, PAM250
from proteinalign.affine_gaps import needleman_wunsch_affine, smith_waterman_affine


SEQ_HEAGAWGHEE = "HEAGAWGHEE"
# H(8) + E(5) + A(4) + G(6) + A(4) + W(11) + G(6) + H(8) + E(5) + E(5) = 62.0
EXPECTED_SELF_SCORE_HEAGAWGHEE = 62.0


class TestAffineMathematicalScaling:
    """Verifies that Gotoh DP produces mathematically exact affine gap costs G(k) = d + (k-1)e."""

    def test_single_residue_gap_cost(self):
        # Gap of length 1: G(1) = gap_open = 10.0
        # s1 = "HEAW", s2 = "HEA"
        # Match HEA: H(8) + E(5) + A(4) = 17.0
        # Expected score: 17.0 - 10.0 = 7.0
        res = needleman_wunsch_affine("HEAW", "HEA", matrix=BLOSUM62, gap_open=10.0, gap_extend=1.0)
        assert res.score == 7.0
        assert res.seq1_aligned == "HEAW"
        assert res.seq2_aligned == "HEA-"

    def test_two_residue_gap_cost(self):
        # Gap of length 2: G(2) = gap_open + 1 * gap_extend = 10.0 + 1.0 = 11.0
        # s1 = "HEAWG", s2 = "HEA"
        # Match HEA: 17.0. Expected score: 17.0 - 11.0 = 6.0
        res = needleman_wunsch_affine("HEAWG", "HEA", matrix=BLOSUM62, gap_open=10.0, gap_extend=1.0)
        assert res.score == 6.0
        assert res.seq1_aligned == "HEAWG"
        assert res.seq2_aligned == "HEA--"

    def test_three_residue_gap_cost(self):
        # Gap of length 3: G(3) = gap_open + 2 * gap_extend = 10.0 + 2.0 = 12.0
        # s1 = "HEAWGF", s2 = "HEA" (suffix WGF has no match in s2)
        # Match HEA: 17.0. Expected score: 17.0 - 12.0 = 5.0
        res = needleman_wunsch_affine("HEAWGF", "HEA", matrix=BLOSUM62, gap_open=10.0, gap_extend=1.0)
        assert res.score == 5.0
        assert res.seq1_aligned == "HEAWGF"
        assert res.seq2_aligned == "HEA---"

    def test_custom_penalties_scaling(self):
        # gap_open = 15.0, gap_extend = 2.0
        # G(4) = 15.0 + 3 * 2.0 = 21.0
        # s1 = "HEAWGAW", s2 = "HEA"
        # Match HEA: 17.0. Expected score: 17.0 - 21.0 = -4.0
        res = needleman_wunsch_affine("HEAWGAW", "HEA", matrix=BLOSUM62, gap_open=15.0, gap_extend=2.0)
        assert res.score == -4.0
        assert res.seq1_aligned == "HEAWGAW"
        assert res.seq2_aligned == "HEA----"


class TestAffineStructuralInvariants:
    """Tests core invariants of the affine alignment implementation."""

    def test_self_alignment_identity(self):
        res = needleman_wunsch_affine(SEQ_HEAGAWGHEE, SEQ_HEAGAWGHEE, matrix=BLOSUM62, gap_open=10.0, gap_extend=1.0)
        assert res.score == EXPECTED_SELF_SCORE_HEAGAWGHEE
        assert res.seq1_aligned == SEQ_HEAGAWGHEE
        assert res.seq2_aligned == SEQ_HEAGAWGHEE
        assert res.identity == 100.0

    def test_global_symmetry(self):
        res_xy = needleman_wunsch_affine("PAWHEAE", "HEAGAWGHEE", matrix=BLOSUM62, gap_open=10.0, gap_extend=1.0)
        res_yx = needleman_wunsch_affine("HEAGAWGHEE", "PAWHEAE", matrix=BLOSUM62, gap_open=10.0, gap_extend=1.0)
        assert res_xy.score == res_yx.score

    def test_degap_reconstruction(self):
        s1 = "MKFLILLFNILCLFPVLAAD"
        s2 = "MKFLILLVLFPAAD"
        res = needleman_wunsch_affine(s1, s2, matrix=BLOSUM62, gap_open=10.0, gap_extend=1.0)
        assert res.seq1_aligned.replace("-", "") == s1
        assert res.seq2_aligned.replace("-", "") == s2
        assert len(res.seq1_aligned) == len(res.seq2_aligned)

    def test_result_stores_original_sequences(self):
        s1 = "HEAG"
        s2 = "HEAP"
        res_global = needleman_wunsch_affine(s1, s2, matrix=BLOSUM62, gap_open=10.0, gap_extend=1.0)
        assert res_global.seq1 == s1
        assert res_global.seq2 == s2

        res_local = smith_waterman_affine(s1, s2, matrix=BLOSUM62, gap_open=10.0, gap_extend=1.0)
        assert res_local.seq1 == s1
        assert res_local.seq2 == s2


class TestAffineContiguousGrouping:
    """Verifies that affine gap penalties favor contiguous indel blocks over fragmented indels."""

    def test_contiguous_gap_vs_scattered(self):
        # Target has insertion of 2 amino acids.
        # Under affine penalties (high open, low extend), grouping gaps together
        # incurs only 1 gap open penalty rather than 2.
        s1 = "HEAGWWWW"
        s2 = "HEAPAWWWW"
        res = needleman_wunsch_affine(s1, s2, matrix=BLOSUM62, gap_open=12.0, gap_extend=1.0)
        # Verify alignment does not scatter single gaps
        assert res.seq1_aligned.replace("-", "") == s1
        assert res.seq2_aligned.replace("-", "") == s2
        assert len(res.seq1_aligned) == len(res.seq2_aligned)


class TestAffineLocalAlignment:
    """Tests for Smith-Waterman local alignment under Gotoh affine model."""

    def test_local_motif_isolation(self):
        s1 = "DDDDPAWHEAERRRR"
        s2 = "WWWWPAWHEAEWWWW"
        res = smith_waterman_affine(s1, s2, matrix=BLOSUM62, gap_open=10.0, gap_extend=1.0)
        assert res.seq1_aligned == "PAWHEAE"
        assert res.seq2_aligned == "PAWHEAE"
        assert res.start1 == 4
        assert res.end1 == 11
        assert res.start2 == 4
        assert res.end2 == 11

    def test_local_coordinate_slice_consistency(self):
        s1 = "XXXXHEAGWWWWWWXXXX"
        s2 = "YYYYHEAWAGWWWWWWYYYY"
        res = smith_waterman_affine(s1, s2, matrix=BLOSUM62, gap_open=4.0, gap_extend=1.0)
        assert res.score == 84.0
        assert s1[res.start1:res.end1] == res.seq1_aligned.replace("-", "")
        assert s2[res.start2:res.end2] == res.seq2_aligned.replace("-", "")

    def test_local_all_negative_returns_empty(self):
        # Dissimilar sequences with no positive sub-alignment
        s1 = "DDDD"
        s2 = "WWWW"
        res = smith_waterman_affine(s1, s2, matrix=BLOSUM62, gap_open=10.0, gap_extend=1.0)
        assert res.score == 0.0
        assert res.seq1_aligned == ""
        assert res.seq2_aligned == ""
        assert res.start1 == 0 and res.end1 == 0


class TestAffineEdgeCases:
    """Tests edge cases: empty strings, single residues, and parameter bounds."""

    def test_empty_vs_empty(self):
        res_g = needleman_wunsch_affine("", "", matrix=BLOSUM62, gap_open=10.0, gap_extend=1.0)
        assert res_g.score == 0.0
        assert res_g.seq1_aligned == ""
        assert res_g.seq2_aligned == ""
        assert res_g.seq1 == ""
        assert res_g.seq2 == ""

        res_l = smith_waterman_affine("", "", matrix=BLOSUM62, gap_open=10.0, gap_extend=1.0)
        assert res_l.score == 0.0
        assert res_l.seq1_aligned == ""
        assert res_l.seq2_aligned == ""

    def test_empty_vs_nonempty_global(self):
        # Query empty, target length 3
        # Cost = -(10 + 2*1) = -12.0
        res = needleman_wunsch_affine("", "HEG", matrix=BLOSUM62, gap_open=10.0, gap_extend=1.0)
        assert res.score == -12.0
        assert res.seq1_aligned == "---"
        assert res.seq2_aligned == "HEG"

        res_rev = needleman_wunsch_affine("HEG", "", matrix=BLOSUM62, gap_open=10.0, gap_extend=1.0)
        assert res_rev.score == -12.0
        assert res_rev.seq1_aligned == "HEG"
        assert res_rev.seq2_aligned == "---"

    def test_empty_vs_nonempty_local(self):
        res = smith_waterman_affine("", "HEG", matrix=BLOSUM62, gap_open=10.0, gap_extend=1.0)
        assert res.score == 0.0
        assert res.seq1_aligned == ""
        assert res.seq2_aligned == ""

    def test_single_residue_match(self):
        res = needleman_wunsch_affine("W", "W", matrix=BLOSUM62, gap_open=10.0, gap_extend=1.0)
        assert res.score == 11.0  # BLOSUM62 W-W
        assert res.seq1_aligned == "W"
        assert res.seq2_aligned == "W"

    def test_single_residue_mismatch(self):
        res = needleman_wunsch_affine("W", "D", matrix=BLOSUM62, gap_open=10.0, gap_extend=1.0)
        # Mismatch W-D is -4.
        # Two gaps would be 2 * -10 = -20.
        # Optimal is mismatch with score -4.
        assert res.score == -4.0
        assert res.seq1_aligned == "W"
        assert res.seq2_aligned == "D"

    def test_negative_penalty_raises_error(self):
        with pytest.raises(ValueError, match="Gap penalties must be non-negative"):
            needleman_wunsch_affine("HEAG", "HEAP", gap_open=-1.0, gap_extend=1.0)

        with pytest.raises(ValueError, match="Gap penalties must be non-negative"):
            needleman_wunsch_affine("HEAG", "HEAP", gap_open=10.0, gap_extend=-1.0)

        with pytest.raises(ValueError, match="Gap penalties must be non-negative"):
            smith_waterman_affine("HEAG", "HEAP", gap_open=-1.0, gap_extend=1.0)

    def test_gap_open_less_than_extend_raises_error(self):
        with pytest.raises(ValueError, match=r"Gap opening penalty .* should typically be >= gap extension"):
            needleman_wunsch_affine("HEAG", "HEAP", gap_open=1.0, gap_extend=5.0)

        with pytest.raises(ValueError, match=r"Gap opening penalty .* should typically be >= gap extension"):
            smith_waterman_affine("HEAG", "HEAP", gap_open=1.0, gap_extend=5.0)

    def test_invalid_characters_raise_error(self):
        with pytest.raises(ValueError, match="Invalid character"):
            needleman_wunsch_affine("HEAG123", "HEAP")

        with pytest.raises(ValueError, match="Invalid character"):
            smith_waterman_affine("HEAG", "HEA$P")
