"""Unit and property tests for profile construction, frequency estimation, and data leakage safeguards."""

import pytest
from proteinalign.profiles import (
    ROBINSON_FREQUENCIES,
    UNIFORM_FREQUENCIES,
    Profile,
    extract_alignment_columns,
    compute_column_counts,
    compute_profile_frequencies,
    build_profile_from_msa,
    split_msa_train_test,
    validate_no_data_leakage,
    validate_background_distribution,
    DataLeakageError,
)
from proteinalign.scoring import AMINO_ACIDS


TOY_MSA = [
    "HEAG",
    "HEAP",
    "HE--",
]


class TestAlignmentColumnExtraction:
    """Tests alignment transpose and sequence validation."""

    def test_valid_column_extraction(self):
        cols = extract_alignment_columns(TOY_MSA)
        assert len(cols) == 4
        assert cols[0] == ["H", "H", "H"]
        assert cols[1] == ["E", "E", "E"]
        assert cols[2] == ["A", "A", "-"]
        assert cols[3] == ["G", "P", "-"]

    def test_empty_sequence_collection_raises_error(self):
        with pytest.raises(ValueError, match="empty sequence collection"):
            extract_alignment_columns([])

    def test_inconsistent_lengths_raise_error(self):
        with pytest.raises(ValueError, match="Inconsistent sequence lengths"):
            extract_alignment_columns(["HEAG", "HEA"])

    def test_invalid_characters_raise_error(self):
        with pytest.raises(ValueError, match="Invalid character '1'"):
            extract_alignment_columns(["HEA1", "HEAG"])

        with pytest.raises(ValueError, match="Invalid character '\\$'"):
            extract_alignment_columns(["HEA$", "HEAG"])


class TestColumnCounts:
    """Tests observed count computation and gap policy."""

    def test_hand_calculated_counts(self):
        cols = extract_alignment_columns(TOY_MSA)
        counts = compute_column_counts(cols)
        assert len(counts) == 4

        # Col 0: H=3, others=0
        assert counts[0]["H"] == 3
        assert sum(counts[0].values()) == 3

        # Col 1: E=3, others=0
        assert counts[1]["E"] == 3
        assert sum(counts[1].values()) == 3

        # Col 2: A=2, gap '-' ignored
        assert counts[2]["A"] == 2
        assert sum(counts[2].values()) == 2

        # Col 3: G=1, P=1, gap '-' ignored
        assert counts[3]["G"] == 1
        assert counts[3]["P"] == 1
        assert sum(counts[3].values()) == 2

    def test_ambiguous_residue_handling(self):
        # 'B' -> Asn (N) + Asp (D); 'Z' -> Gln (Q) + Glu (E)
        msa = ["HB", "HZ"]
        cols = extract_alignment_columns(msa)
        counts = compute_column_counts(cols)
        # Col 1: B and Z -> N(0.5) + D(0.5) + Q(0.5) + E(0.5)
        # Rounded to 1 or 0 integers
        assert sum(counts[1].values()) >= 2


class TestProfileFrequencies:
    """Tests pseudocount regularization, background distributions, and normalization."""

    def test_hand_calculated_frequency_regularization(self):
        # Let column have 3 'H', 0 others, with alpha = 1.0
        counts = [{"H": 3, **{aa: 0 for aa in AMINO_ACIDS if aa != "H"}}]
        bg = dict(ROBINSON_FREQUENCIES)
        q_H = bg["H"]
        q_A = bg["A"]

        # Recurrence: P_{0, a} = (c_{0, a} + alpha * q_a) / (N_0 + alpha)
        # N_0 = 3, alpha = 1.0, denom = 4.0
        expected_H = (3.0 + 1.0 * q_H) / 4.0
        expected_A = (0.0 + 1.0 * q_A) / 4.0
        freqs = compute_profile_frequencies(counts, pseudocount_weight=1.0, background=bg)
        assert pytest.approx(freqs[0]["H"], 1e-6) == expected_H
        assert pytest.approx(freqs[0]["A"], 1e-6) == expected_A

    def test_frequencies_sum_to_one(self):
        cols = extract_alignment_columns(TOY_MSA)
        counts = compute_column_counts(cols)
        freqs = compute_profile_frequencies(counts, pseudocount_weight=1.0)
        for col_f in freqs:
            assert pytest.approx(sum(col_f.values()), 1e-6) == 1.0

    def test_all_gap_column_defaults_to_background(self):
        # Column with only gaps
        counts = [{aa: 0 for aa in AMINO_ACIDS}]
        bg = dict(ROBINSON_FREQUENCIES)
        freqs = compute_profile_frequencies(counts, pseudocount_weight=1.0, background=bg)
        for aa in AMINO_ACIDS:
            assert pytest.approx(freqs[0][aa], 1e-4) == bg[aa]

    def test_negative_pseudocount_weight_raises_error(self):
        counts = [{aa: 0 for aa in AMINO_ACIDS}]
        with pytest.raises(ValueError, match="Pseudocount weight must be non-negative"):
            compute_profile_frequencies(counts, pseudocount_weight=-1.0)

    def test_invalid_background_raises_error(self):
        bad_bg = {"A": 0.5, "R": 0.2}  # Missing other 18 amino acids
        with pytest.raises(ValueError, match="Background distribution missing"):
            validate_background_distribution(bad_bg)

        bad_sum = {aa: 0.1 for aa in AMINO_ACIDS}  # Sum = 2.0 != 1.0
        with pytest.raises(ValueError, match="probabilities must sum to 1.0"):
            validate_background_distribution(bad_sum)


class TestProfileClass:
    """Tests the high-level Profile data structure and constructor."""

    def test_build_profile_from_msa(self):
        prof = build_profile_from_msa(TOY_MSA, name="test_fam", pseudocount_weight=1.0)
        assert prof.name == "test_fam"
        assert prof.alignment_length == 4
        assert prof.num_sequences == 3
        assert len(prof) == 4
        assert prof.consensus == "HEAG"  # Col 3: G has frequency > P or ties
        assert prof.gap_fractions == [0.0, 0.0, 1.0 / 3.0, 1.0 / 3.0]

    def test_profile_dict_input(self):
        msa_dict = {
            "seq_alpha": "HEAG",
            "seq_beta": "HEAP",
        }
        prof = build_profile_from_msa(msa_dict, name="dict_msa")
        assert prof.sequence_ids == ["seq_alpha", "seq_beta"]
        assert prof.num_sequences == 2

    def test_column_accessors(self):
        prof = build_profile_from_msa(TOY_MSA)
        c0 = prof.get_column_counts(0)
        assert c0["H"] == 3
        f0 = prof.get_column_frequencies(0)
        assert f0["H"] > 0.7
        assert prof.get_gap_fraction(0) == 0.0

        with pytest.raises(IndexError):
            prof.get_column_counts(10)

        with pytest.raises(IndexError):
            prof.get_column_frequencies(-1)


class TestDataLeakageSafeguards:
    """Tests train/test MSA splitting and data leakage detection."""

    def test_split_msa_train_test(self):
        msa = {
            "target_1": "HE--AG",
            "ref_A": "HEAPAG",
            "ref_B": "HEAWAG",
        }
        train_msa, test_seq = split_msa_train_test(msa, test_id="target_1")
        assert "target_1" not in train_msa
        assert "ref_A" in train_msa and "ref_B" in train_msa
        # Gaps removed from held-out sequence
        assert test_seq == "HEAG"

    def test_split_missing_target_raises_error(self):
        msa = {"s1": "HEA", "s2": "HEA"}
        with pytest.raises(KeyError, match="Held-out test target 'unknown' not found"):
            split_msa_train_test(msa, test_id="unknown")

    def test_split_too_few_sequences_raises_error(self):
        msa = {"only_one": "HEA"}
        with pytest.raises(ValueError, match="fewer than 2 sequences"):
            split_msa_train_test(msa, test_id="only_one")

    def test_validate_no_data_leakage_detects_overlap(self):
        prof = build_profile_from_msa(
            {"prot_A": "HEAG", "prot_B": "HEAP"},
            name="leak_test",
        )
        # Attempting to evaluate on prot_A which was in training set
        with pytest.raises(DataLeakageError, match="DATA LEAKAGE DETECTED"):
            validate_no_data_leakage(prof, query_id="prot_A")

        # Held-out sequence prot_C succeeds
        validate_no_data_leakage(prof, query_id="prot_C")

        # Explicit override succeeds
        validate_no_data_leakage(prof, query_id="prot_A", allow_leakage=True)
