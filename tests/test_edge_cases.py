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


class TestAdversarialExaminerEdgeCases:
    """Adversarial edge-case audit tests evaluating structural and numerical robustness."""

    def test_substring_local_isolation(self):
        """Property: Local SW extracts exact substring without aligning non-homologous flanks."""
        motif = "HEAW"
        target = "DDDDDDPAWHEAWGKLDDDDDD"
        res = smith_waterman(motif, target, matrix=BLOSUM62, gap_penalty=8.0)
        assert res.score == 8.0 + 5.0 + 4.0 + 11.0  # H(8) + E(5) + A(4) + W(11) = 28.0
        assert res.seq1_aligned == "HEAW"
        assert res.seq2_aligned == "HEAW"
        assert res.start1 == 0
        assert res.end1 == 4
        assert target[res.start2:res.end2] == "HEAW"

    def test_leading_and_trailing_gaps_affine(self):
        """Property: Overhanging ends under Gotoh affine are grouped into single gap runs."""
        s1 = "AAAAHEAWAAAA"
        s2 = "HEAW"
        res = needleman_wunsch_affine(s1, s2, matrix=BLOSUM62, gap_open=10.0, gap_extend=1.0)
        assert res.seq1_aligned.replace("-", "") == s1
        assert res.seq2_aligned.replace("-", "") == s2
        assert len(res.seq1_aligned) == len(res.seq2_aligned)
        # s2 must have leading and trailing gaps
        assert res.seq2_aligned.startswith("----")
        assert res.seq2_aligned.endswith("----")

    def test_long_internal_gap_run(self):
        """Property: Indel of length k costs exactly d + (k-1)e under Gotoh DP."""
        k = 20
        s1 = "H" + ("A" * k) + "E"
        s2 = "HE"
        # Match H(8) + E(5) = 13. Gap of length 20: 10 + 19 * 1 = 29. Score = 13 - 29 = -16.0
        res = needleman_wunsch_affine(s1, s2, matrix=BLOSUM62, gap_open=10.0, gap_extend=1.0)
        assert res.score == -16.0
        assert res.seq1_aligned == s1
        assert res.seq2_aligned == "H" + ("-" * k) + "E"

    def test_tied_optimal_alignments_deterministic(self):
        """Property: Tied paths resolve deterministically across repeated runs."""
        s1 = "AA"
        s2 = "TT"
        results = [needleman_wunsch(s1, s2, matrix=BLOSUM62) for _ in range(5)]
        first_aln1 = results[0].seq1_aligned
        first_aln2 = results[0].seq2_aligned
        for r in results[1:]:
            assert r.seq1_aligned == first_aln1
            assert r.seq2_aligned == first_aln2

    def test_ambiguous_amino_acids_allowed_and_rejected(self):
        """Property: IUPAC codes B, Z, X accepted when allow_ambiguous=True, rejected otherwise."""
        seq_amb = "ABZX"
        seq_std = "ADQE"
        res = needleman_wunsch(seq_amb, seq_std, matrix=BLOSUM62, allow_ambiguous=True)
        assert len(res.seq1_aligned) == len(res.seq2_aligned)

        with pytest.raises(ValueError, match="Invalid character"):
            needleman_wunsch(seq_amb, seq_std, matrix=BLOSUM62, allow_ambiguous=False)

    def test_gap_only_profile_column_fallback(self):
        """Property: Gap-only columns fall back to background frequencies without ZeroDivisionError."""
        from proteinalign.profiles import build_profile_from_msa
        msa = {"s1": "A-C", "s2": "A-C"}
        prof = build_profile_from_msa(msa, pseudocount_weight=0.0)
        assert prof.gap_fractions[1] == 1.0
        # Frequencies at gap column fall back to Robinson background
        assert pytest.approx(prof.frequencies[1]["A"], abs=1e-4) == prof.background["A"]

    def test_invalid_background_distribution_raises(self):
        """Property: Malformed background distributions trigger explicit validation errors."""
        from proteinalign.profiles import validate_background_distribution
        with pytest.raises(ValueError, match="missing standard amino acid"):
            validate_background_distribution({"A": 1.0})

        with pytest.raises(ValueError, match="probabilities must sum to 1.0"):
            from proteinalign.scoring import AMINO_ACIDS
            validate_background_distribution({aa: 0.1 for aa in AMINO_ACIDS})

    def test_data_leakage_assertion(self):
        """Property: validate_no_data_leakage raises DataLeakageError when query matches profile sequence."""
        from proteinalign.profiles import build_profile_from_msa, validate_no_data_leakage, DataLeakageError
        msa = {"seqA": "HEAW", "seqB": "HEAG"}
        prof = build_profile_from_msa(msa, sequence_ids=["seqA", "seqB"])
        with pytest.raises(DataLeakageError):
            validate_no_data_leakage(prof, query_id="seqA", allow_leakage=False)

    def test_traceback_subsequence_coordinate_consistency(self):
        """Property: Local alignment returned coordinates strictly slice original sequence strings."""
        s1 = "MKFLILLFNILCLFPVLAAD"
        s2 = "LLFNILCLFPV"
        res = smith_waterman_affine(s1, s2, matrix=BLOSUM62, gap_open=10.0, gap_extend=1.0)
        assert res.score > 0
        extracted1 = s1[res.start1:res.end1]
        extracted2 = s2[res.start2:res.end2]
        assert extracted1 == res.seq1_aligned.replace("-", "")
        assert extracted2 == res.seq2_aligned.replace("-", "")

