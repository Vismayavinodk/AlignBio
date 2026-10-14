"""Edge cases, boundary conditions, parameter validation, and structural invariants."""

import pytest
from proteinalign import (
    BLOSUM62,
    PAM250,
    needleman_wunsch,
    smith_waterman,
    needleman_wunsch_affine,
    smith_waterman_affine,
)


class TestEmptyAndSingleResidueSequences:
    """Tests for zero-length and single-residue edge cases."""

    def test_both_empty_global(self):
        res = needleman_wunsch("", "")
        assert res.seq1_aligned == ""
        assert res.seq2_aligned == ""
        assert res.score == 0.0

    def test_one_empty_global(self):
        res = needleman_wunsch("HEA", "", gap_penalty=8.0)
        assert res.seq1_aligned == "HEA"
        assert res.seq2_aligned == "---"
        assert res.score == -24.0

    def test_empty_local(self):
        res = smith_waterman("", "HEAG")
        assert res.score == 0.0
        assert res.seq1_aligned == ""
        assert res.seq2_aligned == ""

    def test_single_residue_match(self):
        res = needleman_wunsch("H", "H", matrix=BLOSUM62)
        assert res.score == 8.0
        assert res.seq1_aligned == "H"
        assert res.seq2_aligned == "H"

    def test_single_residue_mismatch(self):
        res = needleman_wunsch("H", "P", matrix=BLOSUM62)
        # BLOSUM62 H-P = -2
        assert res.score == -2.0
        assert res.seq1_aligned == "H"
        assert res.seq2_aligned == "P"


class TestRepeatedResiduesAndLengthDiscrepancy:
    """Tests for homopolymeric tracts and unbalanced lengths."""

    def test_repeated_residues(self):
        res = needleman_wunsch("AAAA", "AA", matrix=BLOSUM62, gap_penalty=8.0)
        assert res.seq1_aligned.replace("-", "") == "AAAA"
        assert res.seq2_aligned.replace("-", "") == "AA"
        assert len(res.seq1_aligned) == len(res.seq2_aligned)

    def test_extreme_length_difference(self):
        short_seq = "A"
        long_seq = "ACDEFGHIKLMNPQRSTVWY"
        res = needleman_wunsch(short_seq, long_seq, matrix=BLOSUM62, gap_penalty=8.0)
        assert res.seq1_aligned.replace("-", "") == short_seq
        assert res.seq2_aligned.replace("-", "") == long_seq
        assert len(res.seq1_aligned) == len(res.seq2_aligned)


class TestCompletelyDissimilarSequences:
    """Tests when no positive local alignment exists."""

    def test_no_positive_local_alignment(self):
        # Two sequences composed of pairs with exclusively negative BLOSUM scores
        # W vs D (-4), C vs D (-3)
        s1 = "WWWW"
        s2 = "DDDD"
        res = smith_waterman(s1, s2, matrix=BLOSUM62, gap_penalty=8.0)
        assert res.score == 0.0
        assert res.seq1_aligned == ""
        assert res.seq2_aligned == ""


class TestParameterValidationAndExceptions:
    """Tests for bad inputs, negative penalties, and unsupported symbols."""

    def test_invalid_amino_acid_symbol(self):
        with pytest.raises(ValueError, match="Invalid character"):
            needleman_wunsch("HEAG123", "HEAP")

        with pytest.raises(ValueError, match="Invalid character"):
            smith_waterman("HEA$G", "HEAP")

    def test_negative_gap_penalty(self):
        with pytest.raises(ValueError, match="Gap penalty must be non-negative"):
            needleman_wunsch("HEAG", "HEAP", gap_penalty=-5.0)

    def test_inconsistent_affine_penalties(self):
        with pytest.raises(ValueError, match="should typically be >="):
            needleman_wunsch_affine("HEAG", "HEAP", gap_open=1.0, gap_extend=10.0)


class TestStructuralInvariants:
    """Formal mathematical invariants every valid alignment must satisfy."""

    @pytest.mark.parametrize("s1, s2", [
        ("MKFLILLFNILCLFPVLAAD", "MKFLILLVLFPAAD"),
        ("ACDEFGHIKLMNPQRSTVWY", "ACDEFGHIKLMNPQRSTVWY"),
        ("HEAGAWGHEE", "PAWHEAE"),
    ])
    def test_degap_and_length_invariants(self, s1, s2):
        for aligner in [needleman_wunsch, needleman_wunsch_affine]:
            res = aligner(s1, s2)
            assert res.seq1_aligned.replace("-", "") == s1
            assert res.seq2_aligned.replace("-", "") == s2
            assert len(res.seq1_aligned) == len(res.seq2_aligned)
