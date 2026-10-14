"""Unit tests for substitution matrices, validation, and residue scoring."""

import pytest
from proteinalign.scoring import (
    AMINO_ACIDS,
    BLOSUM62,
    PAM250,
    SubstitutionMatrix,
    is_valid_amino_acid,
    validate_sequence,
)


class TestSimpleIdentityScoring:
    """Tests for simple match/mismatch substitution matrix."""

    def test_match_and_mismatch(self):
        mat = SubstitutionMatrix.create_identity(match=3.0, mismatch=-2.0)
        assert mat.score("A", "A") == 3.0
        assert mat.score("a", "A") == 3.0
        assert mat.score("A", "G") == -2.0
        assert mat.score("W", "C") == -2.0

    def test_custom_values(self):
        mat = SubstitutionMatrix.create_identity(match=5.0, mismatch=-4.0)
        assert mat.score("H", "H") == 5.0
        assert mat.score("H", "P") == -4.0


class TestBLOSUM62Scoring:
    """Tests for NCBI standard BLOSUM62 matrix."""

    def test_diagonal_self_scores(self):
        assert BLOSUM62.score("W", "W") == 11.0
        assert BLOSUM62.score("C", "C") == 9.0
        assert BLOSUM62.score("H", "H") == 8.0
        assert BLOSUM62.score("P", "P") == 7.0
        assert BLOSUM62.score("A", "A") == 4.0

    def test_conservative_substitutions(self):
        # Chemically similar pairs should have positive log-odds
        assert BLOSUM62.score("K", "R") == 2.0  # Basic
        assert BLOSUM62.score("I", "V") == 3.0  # Aliphatic / hydrophobic
        assert BLOSUM62.score("F", "Y") == 3.0  # Aromatic

    def test_non_conservative_substitutions(self):
        assert BLOSUM62.score("W", "C") == -2.0
        assert BLOSUM62.score("D", "R") == -2.0
        assert BLOSUM62.score("G", "W") == -2.0

    def test_symmetry(self):
        for a in AMINO_ACIDS[:10]:
            for b in AMINO_ACIDS[:10]:
                assert BLOSUM62.score(a, b) == BLOSUM62.score(b, a)


class TestPAM250Scoring:
    """Tests for Dayhoff PAM250 matrix."""

    def test_pam250_scores(self):
        assert PAM250.score("W", "W") == 17.0
        assert PAM250.score("C", "C") == 12.0
        assert PAM250.score("Y", "Y") == 10.0
        assert PAM250.score("A", "A") == 2.0

    def test_pam250_symmetry(self):
        for a in AMINO_ACIDS[:10]:
            for b in AMINO_ACIDS[:10]:
                assert PAM250.score(a, b) == PAM250.score(b, a)


class TestAmbiguousResiduesAndPolicy:
    """Tests for handling unknown residues (X), ambiguous residues (B, Z), and stop codons."""

    def test_unknown_residue_x(self):
        # Aligning with unknown residue X incurs mild default penalty (-1.0)
        assert BLOSUM62.score("A", "X") == -1.0
        assert BLOSUM62.score("X", "W") == -1.0

    def test_ambiguous_b_asx(self):
        # B is Aspartate (D) or Asparagine (N). Score is average of D and N.
        # BLOSUM62: A-D = -2, A-N = -2 => average = -2.0
        assert BLOSUM62.score("A", "B") == -2.0

    def test_ambiguous_z_glx(self):
        # Z is Glutamate (E) or Glutamine (Q). Score is average of E and Q.
        # BLOSUM62: A-E = -1, A-Q = -1 => average = -1.0
        assert BLOSUM62.score("A", "Z") == -1.0


class TestResidueValidation:
    """Tests for residue character validation and exception handling."""

    def test_valid_characters(self):
        assert is_valid_amino_acid("A")
        assert is_valid_amino_acid("w")
        assert is_valid_amino_acid("X", allow_ambiguous=True)
        assert not is_valid_amino_acid("1")
        assert not is_valid_amino_acid("$")
        assert not is_valid_amino_acid("AA")  # Length != 1

    def test_sequence_validation_exceptions(self):
        validate_sequence("ACDEFGHIKLMNPQRSTVWY")  # Valid, no exception

        with pytest.raises(ValueError, match="Invalid character '1'"):
            validate_sequence("ACD1EFG")

        with pytest.raises(ValueError, match="Invalid character '!'"):
            validate_sequence("HEAG!P")
