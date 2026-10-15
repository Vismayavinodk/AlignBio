"""Unit and property tests for PSSM construction, log-odds scoring, and profile alignment."""

import math
import pytest
from proteinalign.profiles import (
    ROBINSON_FREQUENCIES,
    UNIFORM_FREQUENCIES,
    build_profile_from_msa,
    DataLeakageError,
)
from proteinalign.pssm import PSSM
from proteinalign.scoring import AMINO_ACIDS


TOY_MSA = {
    "prot_1": "HEAGAWG",
    "prot_2": "HEAPAWG",
    "prot_3": "HEAGAWA",
}


class TestPSSMConstructionAndScoring:
    """Verifies analytical log-odds score computation, scaling, and ambiguity handling."""

    def test_hand_calculated_pssm_log_odds(self):
        # Build profile from toy MSA with uniform background for easy hand verification
        bg = dict(UNIFORM_FREQUENCIES)  # q_a = 0.05 for all 20 amino acids
        # Column 0 has 3 'H', N_0 = 3.
        # With alpha = 1.0: P_{0, H} = (3 + 1.0 * 0.05) / (3 + 1.0) = 3.05 / 4.0 = 0.7625
        # Ratio = P_{0, H} / q_H = 0.7625 / 0.05 = 15.25
        # Under scale = 2.0, log_base = 2.0 (half-bits):
        # S_{0, H} = 2.0 * log2(15.25) = 2.0 * 3.9307 = 7.861
        pssm = PSSM.from_msa(
            TOY_MSA,
            pseudocount_weight=1.0,
            background=bg,
            scale=2.0,
            log_base=2.0,
        )
        assert len(pssm) == 7
        assert pssm.consensus.startswith("HEA")

        expected_score_H0 = round(2.0 * math.log2(0.7625 / 0.05), 3)
        assert pytest.approx(pssm.score(0, "H"), abs=1e-3) == expected_score_H0

    def test_ambiguous_residue_scoring(self):
        pssm = PSSM.from_msa(TOY_MSA, scale=2.0, log_base=2.0)
        # B = 0.5 * score(N) + 0.5 * score(D)
        score_B = pssm.score(0, "B")
        score_N = pssm.score(0, "N")
        score_D = pssm.score(0, "D")
        assert pytest.approx(score_B, 1e-4) == 0.5 * score_N + 0.5 * score_D

        # Z = 0.5 * score(Q) + 0.5 * score(E)
        score_Z = pssm.score(0, "Z")
        score_Q = pssm.score(0, "Q")
        score_E = pssm.score(0, "E")
        assert pytest.approx(score_Z, 1e-4) == 0.5 * score_Q + 0.5 * score_E

        # X = expected score under background
        score_X = pssm.score(0, "X")
        expected_bg = sum(pssm.background[aa] * pssm.score(0, aa) for aa in AMINO_ACIDS)
        assert pytest.approx(score_X, 1e-4) == expected_bg

    def test_invalid_parameters_raise_error(self):
        prof = build_profile_from_msa(TOY_MSA)
        with pytest.raises(ValueError, match="scale parameter lambda must be positive"):
            PSSM.from_profile(prof, scale=-1.0)

        with pytest.raises(ValueError, match="Logarithm base must be > 1.0"):
            PSSM.from_profile(prof, log_base=1.0)

    def test_out_of_bounds_and_invalid_residue(self):
        pssm = PSSM.from_msa(TOY_MSA)
        with pytest.raises(IndexError):
            pssm.score(100, "A")

        with pytest.raises(ValueError, match="Unknown amino acid character"):
            pssm.score(0, "1")


class TestProfileToSequenceGlobalAlignment:
    """Tests global Needleman-Wunsch profile-to-sequence alignment (linear and affine)."""

    def test_global_linear_alignment_identity(self):
        pssm = PSSM.from_msa(TOY_MSA, pseudocount_weight=1.0)
        # Align consensus "HEAGAWG" against itself
        res = pssm.align_global("HEAGAWG", gap_penalty=8.0)
        assert res.seq2_aligned == "HEAGAWG"
        assert res.score > 20.0  # Strongly positive match score
        assert pssm.recompute_score(res) == res.score

    def test_global_linear_single_deletion(self):
        pssm = PSSM.from_msa(TOY_MSA, pseudocount_weight=1.0)
        # Query has extra residue at end
        query = "HEAGAWA"
        res = pssm.align_global(query, gap_penalty=8.0)
        assert res.seq1_aligned.replace("-", "") == pssm.consensus
        assert res.seq2_aligned.replace("-", "") == query
        assert len(res.seq1_aligned) == len(res.seq2_aligned)
        assert pytest.approx(pssm.recompute_score(res), 1e-3) == res.score

    def test_global_affine_alignment_and_recomputation(self):
        pssm = PSSM.from_msa(TOY_MSA, pseudocount_weight=1.0)
        # Query has a 2-residue insertion
        query = "HEAGPPAWG"
        res = pssm.align_global(query, gap_open=10.0, gap_extend=1.0)
        assert res.alignment_type == "pssm_global_affine"
        assert len(res.seq1_aligned) == len(res.seq2_aligned)
        assert res.seq1_aligned.replace("-", "") == pssm.consensus
        assert res.seq2_aligned.replace("-", "") == query
        assert pytest.approx(pssm.recompute_score(res), 1e-3) == res.score


class TestProfileToSequenceLocalAlignment:
    """Tests local Smith-Waterman profile-to-sequence alignment (linear and affine)."""

    def test_local_linear_motif_isolation(self):
        # Profile has consensus "HEAGAWG"
        pssm = PSSM.from_msa(TOY_MSA, pseudocount_weight=1.0)
        # Embed consensus in dissimilar flanking residues (D, W, R)
        query = "DDDDHEAGAWGRRRR"
        res = pssm.align_local(query, gap_penalty=8.0)
        assert res.seq2_aligned == "HEAGAWG"
        assert res.start2 == 4
        assert res.end2 == 11
        assert query[res.start2:res.end2] == "HEAGAWG"
        assert pytest.approx(pssm.recompute_score(res), 1e-3) == res.score

    def test_local_affine_internal_gap(self):
        pssm = PSSM.from_msa(TOY_MSA, pseudocount_weight=1.0)
        query = "DDDDHEAGPPAWGRRRR"
        res = pssm.align_local(query, gap_open=4.0, gap_extend=1.0)
        assert res.alignment_type == "pssm_local_affine"
        assert query[res.start2:res.end2] == res.seq2_aligned.replace("-", "")
        assert pytest.approx(pssm.recompute_score(res), 1e-3) == res.score

    def test_local_all_negative_returns_empty(self):
        pssm = PSSM.from_msa(TOY_MSA, pseudocount_weight=1.0)
        # All non-matching residues with negative log-odds
        query = "DDDD"
        res = pssm.align_local(query, gap_penalty=10.0)
        assert res.score == 0.0
        assert res.seq1_aligned == ""
        assert res.seq2_aligned == ""


class TestDataLeakageDetectionInAlignment:
    """Verifies that PSSM alignment refuses to align training sequences without explicit override."""

    def test_alignment_raises_data_leakage_error(self):
        pssm = PSSM.from_msa(TOY_MSA)
        # prot_1 was one of the training sequences
        with pytest.raises(DataLeakageError, match="DATA LEAKAGE DETECTED"):
            pssm.align_global("HEAGAWG", query_id="prot_1")

        with pytest.raises(DataLeakageError, match="DATA LEAKAGE DETECTED"):
            pssm.align_local("HEAGAWG", query_id="prot_2")

    def test_allow_leakage_permits_alignment(self):
        pssm = PSSM.from_msa(TOY_MSA)
        res = pssm.align_global("HEAGAWG", query_id="prot_1", allow_leakage=True)
        assert res.score > 0.0


class TestPSSMEdgeCases:
    """Verifies boundary conditions, empty sequences, and single-column profiles."""

    def test_empty_query_sequence(self):
        pssm = PSSM.from_msa(TOY_MSA)
        res_g = pssm.align_global("", gap_penalty=8.0)
        assert res_g.score == -len(pssm) * 8.0
        assert res_g.seq1_aligned == pssm.consensus
        assert res_g.seq2_aligned == "-" * len(pssm)

        res_l = pssm.align_local("", gap_penalty=8.0)
        assert res_l.score == 0.0
        assert res_l.seq1_aligned == ""
        assert res_l.seq2_aligned == ""

    def test_single_column_profile(self):
        msa = ["H", "H", "H"]
        pssm = PSSM.from_msa(msa)
        assert len(pssm) == 1
        res = pssm.align_global("H", gap_penalty=8.0)
        assert res.seq1_aligned == "H"
        assert res.seq2_aligned == "H"
        assert res.score == pssm.score(0, "H")

    def test_determinism_across_runs(self):
        pssm = PSSM.from_msa(TOY_MSA)
        query = "HEAPAWG"
        res1 = pssm.align_global(query, gap_open=10.0, gap_extend=1.0)
        res2 = pssm.align_global(query, gap_open=10.0, gap_extend=1.0)
        assert res1.score == res2.score
        assert res1.seq1_aligned == res2.seq1_aligned
        assert res1.seq2_aligned == res2.seq2_aligned
