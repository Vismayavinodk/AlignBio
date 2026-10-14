"""Substitution scoring matrices and sequence validation for protein alignment."""

from typing import Dict, List, Optional, Set, Tuple, Union

# Standard 20 canonical amino acid residues
STANDARD_AMINO_ACIDS: str = "ARNDCQEGHILKMFPSTWYV"
AMINO_ACIDS: str = STANDARD_AMINO_ACIDS

# Ambiguous amino acid symbols recognized by IUPAC
# B: Aspartate or Asparagine (D or N)
# Z: Glutamate or Glutamine (E or Q)
# J: Leucine or Isoleucine (L or I)
# X: Any / unknown amino acid
# *: Translation stop codon
AMBIGUOUS_AMINO_ACIDS: str = "BZJX*"

ALL_VALID_AMINO_ACIDS: Set[str] = set(STANDARD_AMINO_ACIDS + AMBIGUOUS_AMINO_ACIDS)


def is_valid_amino_acid(char: str, allow_ambiguous: bool = True) -> bool:
    """Checks whether a character is a valid amino acid residue symbol.

    Args:
        char: Single character string.
        allow_ambiguous: If True, recognizes IUPAC ambiguous codes (B, Z, J, X, *).

    Returns:
        True if character is a recognized residue symbol, False otherwise.
    """
    if len(char) != 1:
        return False
    u = char.upper()
    if u in STANDARD_AMINO_ACIDS:
        return True
    return allow_ambiguous and (u in AMBIGUOUS_AMINO_ACIDS)


def validate_sequence(seq: str, allow_ambiguous: bool = True) -> None:
    """Validates that a sequence consists exclusively of valid residue characters.

    Args:
        seq: Protein sequence string.
        allow_ambiguous: Whether ambiguous IUPAC codes are permitted.

    Raises:
        ValueError: If sequence contains invalid non-amino acid characters.
    """
    for idx, char in enumerate(seq):
        if char == "-":
            continue
        if not is_valid_amino_acid(char, allow_ambiguous=allow_ambiguous):
            raise ValueError(
                f"Invalid character '{char}' at index {idx} in sequence. "
                f"Expected standard IUPAC amino acid residue."
            )


# Standard NCBI BLOSUM62 log-odds values
_BLOSUM62_RAW: Dict[Tuple[str, str], int] = {
    ("A", "A"): 4, ("A", "R"): -1, ("A", "N"): -2, ("A", "D"): -2, ("A", "C"): 0, ("A", "Q"): -1, ("A", "E"): -1, ("A", "G"): 0, ("A", "H"): -2, ("A", "I"): -1, ("A", "L"): -1, ("A", "K"): -1, ("A", "M"): -1, ("A", "F"): -2, ("A", "P"): -1, ("A", "S"): 1, ("A", "T"): 0, ("A", "W"): -3, ("A", "Y"): -2, ("A", "V"): 0,
    ("R", "R"): 5, ("R", "N"): 0, ("R", "D"): -2, ("R", "C"): -3, ("R", "Q"): 1, ("R", "E"): 0, ("R", "G"): -2, ("R", "H"): 0, ("R", "I"): -3, ("R", "L"): -2, ("R", "K"): 2, ("R", "M"): -1, ("R", "F"): -3, ("R", "P"): -2, ("R", "S"): -1, ("R", "T"): -1, ("R", "W"): -3, ("R", "Y"): -2, ("R", "V"): -3,
    ("N", "N"): 6, ("N", "D"): 1, ("N", "C"): -3, ("N", "Q"): 0, ("N", "E"): 0, ("N", "G"): 0, ("N", "H"): 1, ("N", "I"): -3, ("N", "L"): -3, ("N", "K"): 0, ("N", "M"): -2, ("N", "F"): -3, ("N", "P"): -2, ("N", "S"): 1, ("N", "T"): 0, ("N", "W"): -4, ("N", "Y"): -2, ("N", "V"): -3,
    ("D", "D"): 6, ("D", "C"): -3, ("D", "Q"): 0, ("D", "E"): 2, ("D", "G"): -1, ("D", "H"): -1, ("D", "I"): -3, ("D", "L"): -4, ("D", "K"): -1, ("D", "M"): -3, ("D", "F"): -3, ("D", "P"): -1, ("D", "S"): 0, ("D", "T"): -1, ("D", "W"): -4, ("D", "Y"): -3, ("D", "V"): -3,
    ("C", "C"): 9, ("C", "Q"): -3, ("C", "E"): -4, ("C", "G"): -3, ("C", "H"): -3, ("C", "I"): -1, ("C", "L"): -1, ("C", "K"): -3, ("C", "M"): -1, ("C", "F"): -2, ("C", "P"): -3, ("C", "S"): -1, ("C", "T"): -1, ("C", "W"): -2, ("C", "Y"): -2, ("C", "V"): -1,
    ("Q", "Q"): 5, ("Q", "E"): 2, ("Q", "G"): -2, ("Q", "H"): 0, ("Q", "I"): -3, ("Q", "L"): -2, ("Q", "K"): 1, ("Q", "M"): 0, ("Q", "F"): -3, ("Q", "P"): -1, ("Q", "S"): 0, ("Q", "T"): -1, ("Q", "W"): -2, ("Q", "Y"): -1, ("Q", "V"): -2,
    ("E", "E"): 5, ("E", "G"): -2, ("E", "H"): 0, ("E", "I"): -3, ("E", "L"): -3, ("E", "K"): 1, ("E", "M"): -2, ("E", "F"): -3, ("E", "P"): -1, ("E", "S"): 0, ("E", "T"): -1, ("E", "W"): -3, ("E", "Y"): -2, ("E", "V"): -2,
    ("G", "G"): 6, ("G", "H"): -2, ("G", "I"): -4, ("G", "L"): -4, ("G", "K"): -2, ("G", "M"): -3, ("G", "F"): -3, ("G", "P"): -2, ("G", "S"): 0, ("G", "T"): -2, ("G", "W"): -2, ("G", "Y"): -3, ("G", "V"): -3,
    ("H", "H"): 8, ("H", "I"): -3, ("H", "L"): -3, ("H", "K"): -1, ("H", "M"): -2, ("H", "F"): -1, ("H", "P"): -2, ("H", "S"): -1, ("H", "T"): -2, ("H", "W"): -2, ("H", "Y"): 2, ("H", "V"): -3,
    ("I", "I"): 4, ("I", "L"): 2, ("I", "K"): -3, ("I", "M"): 1, ("I", "F"): 0, ("I", "P"): -3, ("I", "S"): -2, ("I", "T"): -1, ("I", "W"): -3, ("I", "Y"): -1, ("I", "V"): 3,
    ("L", "L"): 4, ("L", "K"): -2, ("L", "M"): 2, ("L", "F"): 0, ("L", "P"): -3, ("L", "S"): -2, ("L", "T"): -1, ("L", "W"): -2, ("L", "Y"): -1, ("L", "V"): 1,
    ("K", "K"): 5, ("K", "M"): -1, ("K", "F"): -3, ("K", "P"): -1, ("K", "S"): 0, ("K", "T"): -1, ("K", "W"): -3, ("K", "Y"): -2, ("K", "V"): -2,
    ("M", "M"): 5, ("M", "F"): 0, ("M", "P"): -2, ("M", "S"): -1, ("M", "T"): -1, ("M", "W"): -1, ("M", "Y"): -1, ("M", "V"): 1,
    ("F", "F"): 6, ("F", "P"): -4, ("F", "S"): -2, ("F", "T"): -2, ("F", "W"): 1, ("F", "Y"): 3, ("F", "V"): -1,
    ("P", "P"): 7, ("P", "S"): -1, ("P", "T"): -1, ("P", "W"): -4, ("P", "Y"): -3, ("P", "V"): -2,
    ("S", "S"): 4, ("S", "T"): 1, ("S", "W"): -3, ("S", "Y"): -2, ("S", "V"): -2,
    ("T", "T"): 5, ("T", "W"): -2, ("T", "Y"): -2, ("T", "V"): 0,
    ("W", "W"): 11, ("W", "Y"): 2, ("W", "V"): -3,
    ("Y", "Y"): 7, ("Y", "V"): -1,
    ("V", "V"): 4,
}

# Standard Dayhoff PAM250 log-odds values
_PAM250_RAW: Dict[Tuple[str, str], int] = {
    ("A", "A"): 2, ("A", "R"): -2, ("A", "N"): 0, ("A", "D"): 0, ("A", "C"): -2, ("A", "Q"): 0, ("A", "E"): 0, ("A", "G"): 1, ("A", "H"): -1, ("A", "I"): -1, ("A", "L"): -2, ("A", "K"): -1, ("A", "M"): -1, ("A", "F"): -3, ("A", "P"): 1, ("A", "S"): 1, ("A", "T"): 1, ("A", "W"): -6, ("A", "Y"): -3, ("A", "V"): 0,
    ("R", "R"): 6, ("R", "N"): 0, ("R", "D"): -1, ("R", "C"): -4, ("R", "Q"): 1, ("R", "E"): -1, ("R", "G"): -3, ("R", "H"): 2, ("R", "I"): -2, ("R", "L"): -3, ("R", "K"): 3, ("R", "M"): 0, ("R", "F"): -4, ("R", "P"): 0, ("R", "S"): 0, ("R", "T"): -1, ("R", "W"): 2, ("R", "Y"): -4, ("R", "V"): -2,
    ("N", "N"): 2, ("N", "D"): 2, ("N", "C"): -4, ("N", "Q"): 1, ("N", "E"): 1, ("N", "G"): 0, ("N", "H"): 2, ("N", "I"): -2, ("N", "L"): -3, ("N", "K"): 1, ("N", "M"): -2, ("N", "F"): -3, ("N", "P"): -1, ("N", "S"): 1, ("N", "T"): 0, ("N", "W"): -4, ("N", "Y"): -2, ("N", "V"): -2,
    ("D", "D"): 4, ("D", "C"): -5, ("D", "Q"): 2, ("D", "E"): 3, ("D", "G"): 1, ("D", "H"): 1, ("D", "I"): -2, ("D", "L"): -4, ("D", "K"): 0, ("D", "M"): -3, ("D", "F"): -6, ("D", "P"): -1, ("D", "S"): 0, ("D", "T"): 0, ("D", "W"): -7, ("D", "Y"): -4, ("D", "V"): -2,
    ("C", "C"): 12, ("C", "Q"): -5, ("C", "E"): -5, ("C", "G"): -3, ("C", "H"): -3, ("C", "I"): -2, ("C", "L"): -6, ("C", "K"): -5, ("C", "M"): -5, ("C", "F"): -4, ("C", "P"): -3, ("C", "S"): 0, ("C", "T"): -2, ("C", "W"): -8, ("C", "Y"): 0, ("C", "V"): -2,
    ("Q", "Q"): 4, ("Q", "E"): 2, ("Q", "G"): -1, ("Q", "H"): 3, ("Q", "I"): -2, ("Q", "L"): -2, ("Q", "K"): 1, ("Q", "M"): -1, ("Q", "F"): -5, ("Q", "P"): 0, ("Q", "S"): -1, ("Q", "T"): -1, ("Q", "W"): -5, ("Q", "Y"): -4, ("Q", "V"): -2,
    ("E", "E"): 4, ("E", "G"): 0, ("E", "H"): 1, ("E", "I"): -2, ("E", "L"): -3, ("E", "K"): 0, ("E", "M"): -2, ("E", "F"): -5, ("E", "P"): -1, ("E", "S"): 0, ("E", "T"): -1, ("E", "W"): -7, ("E", "Y"): -4, ("E", "V"): -2,
    ("G", "G"): 5, ("G", "H"): -2, ("G", "I"): -3, ("G", "L"): -4, ("G", "K"): -2, ("G", "M"): -3, ("G", "F"): -5, ("G", "P"): 0, ("G", "S"): 1, ("G", "T"): 0, ("G", "W"): -7, ("G", "Y"): -5, ("G", "V"): -1,
    ("H", "H"): 6, ("H", "I"): -2, ("H", "L"): -2, ("H", "K"): 0, ("H", "M"): -2, ("H", "F"): -2, ("H", "P"): 0, ("H", "S"): -1, ("H", "T"): -1, ("H", "W"): -3, ("H", "Y"): 0, ("H", "V"): -2,
    ("I", "I"): 5, ("I", "L"): 2, ("I", "K"): -2, ("I", "M"): 2, ("I", "F"): 1, ("I", "P"): -2, ("I", "S"): -1, ("I", "T"): 0, ("I", "W"): -5, ("I", "Y"): -1, ("I", "V"): 4,
    ("L", "L"): 6, ("L", "K"): -3, ("L", "M"): 4, ("L", "F"): 2, ("L", "P"): -3, ("L", "S"): -3, ("L", "T"): -2, ("L", "W"): -2, ("L", "Y"): -1, ("L", "V"): 2,
    ("K", "K"): 5, ("K", "M"): 0, ("K", "F"): -5, ("K", "P"): -1, ("K", "S"): 0, ("K", "T"): 0, ("K", "W"): -3, ("K", "Y"): -4, ("K", "V"): -2,
    ("M", "M"): 6, ("M", "F"): 0, ("M", "P"): -2, ("M", "S"): -2, ("M", "T"): -1, ("M", "W"): -4, ("M", "Y"): -2, ("M", "V"): 2,
    ("F", "F"): 9, ("F", "P"): -5, ("F", "S"): -3, ("F", "T"): -3, ("F", "W"): 0, ("F", "Y"): 7, ("F", "V"): -1,
    ("P", "P"): 6, ("P", "S"): 1, ("P", "T"): 0, ("P", "W"): -6, ("P", "Y"): -5, ("P", "V"): -1,
    ("S", "S"): 2, ("S", "T"): 1, ("S", "W"): -2, ("S", "Y"): -3, ("S", "V"): -1,
    ("T", "T"): 3, ("T", "W"): -5, ("T", "Y"): -3, ("T", "V"): 0,
    ("W", "W"): 17, ("W", "Y"): 0, ("W", "V"): -6,
    ("Y", "Y"): 10, ("Y", "V"): -2,
    ("V", "V"): 4,
}


def _build_symmetric_matrix(raw: Dict[Tuple[str, str], int]) -> Dict[Tuple[str, str], float]:
    """Generates a complete bidirectional symmetric score dictionary."""
    full = {}
    for (a, b), val in raw.items():
        full[(a, b)] = float(val)
        full[(b, a)] = float(val)
    return full


class SubstitutionMatrix:
    """Represents a pairwise amino acid substitution scoring matrix.

    Provides a clean, uniform interface for retrieving substitution scores,
    with built-in validation for ambiguous residues (e.g. X, B, Z) and unsupported characters.
    """

    def __init__(
        self,
        matrix_dict: Dict[Tuple[str, str], float],
        name: str = "custom",
        default_score: float = -4.0,
        unknown_score: float = -1.0,
    ) -> None:
        """Initializes a substitution matrix.

        Args:
            matrix_dict: Mapping of (residue1, residue2) tuples to numerical scores.
            name: Descriptive identifier (e.g. 'BLOSUM62', 'PAM250').
            default_score: Score applied to unlisted non-identical standard residue pairs.
            unknown_score: Score applied when aligning with unknown/ambiguous residue 'X'.
        """
        self.matrix = matrix_dict
        self.name = name
        self.default_score = default_score
        self.unknown_score = unknown_score

    def score(self, a: str, b: str) -> float:
        """Returns the substitution score between residues a and b.

        Case-insensitive. Handles ambiguous residues (X, B, Z) and unknown symbols
        according to the documented policy.

        Args:
            a: First amino acid character.
            b: Second amino acid character.

        Returns:
            Numerical substitution score.
        """
        ua = a.upper()
        ub = b.upper()

        # Check explicit matrix entries
        if (ua, ub) in self.matrix:
            return self.matrix[(ua, ub)]

        # Identical residue fallback
        if ua == ub:
            return 1.0

        # Unknown / ambiguous residue policy:
        # X represents an undetermined amino acid; alignment with X incurs a mild penalty (-1.0).
        if ua == "X" or ub == "X":
            return self.unknown_score

        # Ambiguous B (Asx: D or N) -> average of D and N scores
        if ua == "B" or ub == "B":
            other = ub if ua == "B" else ua
            s_d = self.matrix.get((other, "D"), self.matrix.get(("D", other), self.default_score))
            s_n = self.matrix.get((other, "N"), self.matrix.get(("N", other), self.default_score))
            return (s_d + s_n) / 2.0

        # Ambiguous Z (Glx: E or Q) -> average of E and Q scores
        if ua == "Z" or ub == "Z":
            other = ub if ua == "Z" else ua
            s_e = self.matrix.get((other, "E"), self.matrix.get(("E", other), self.default_score))
            s_q = self.matrix.get((other, "Q"), self.matrix.get(("Q", other), self.default_score))
            return (s_e + s_q) / 2.0

        return self.default_score

    @classmethod
    def load_blosum62(cls) -> "SubstitutionMatrix":
        """Loads the standard NCBI BLOSUM62 substitution matrix."""
        return cls(_build_symmetric_matrix(_BLOSUM62_RAW), name="BLOSUM62", default_score=-4.0, unknown_score=-1.0)

    @classmethod
    def load_pam250(cls) -> "SubstitutionMatrix":
        """Loads the standard Dayhoff PAM250 substitution matrix."""
        return cls(_build_symmetric_matrix(_PAM250_RAW), name="PAM250", default_score=-4.0, unknown_score=-1.0)

    @classmethod
    def create_identity(cls, match: float = 2.0, mismatch: float = -1.0) -> "SubstitutionMatrix":
        """Creates a simple match reward and mismatch penalty matrix for test verification.

        Args:
            match: Score awarded for identical residues (positive float).
            mismatch: Penalty incurred for non-identical residues (negative float).

        Returns:
            SubstitutionMatrix instance.
        """
        mat = {}
        for a in STANDARD_AMINO_ACIDS:
            for b in STANDARD_AMINO_ACIDS:
                mat[(a, b)] = match if a == b else mismatch
        return cls(mat, name=f"Identity({match:+g}/{mismatch:+g})", default_score=mismatch, unknown_score=mismatch)

    @classmethod
    def from_biopython(cls, matrix_name: str = "BLOSUM62") -> "SubstitutionMatrix":
        """Loads a substitution matrix from Biopython's data directory if available.

        Args:
            matrix_name: Name of standard matrix (e.g. 'BLOSUM62', 'PAM250').

        Returns:
            SubstitutionMatrix instance.
        """
        try:
            from Bio.Align import substitution_matrices
            bio_mat = substitution_matrices.load(matrix_name)
            mat_dict = {}
            for k, val in bio_mat.items():
                mat_dict[(k[0].upper(), k[1].upper())] = float(val)
                mat_dict[(k[1].upper(), k[0].upper())] = float(val)
            return cls(mat_dict, name=f"Biopython_{matrix_name}", default_score=-4.0)
        except Exception:
            if matrix_name.upper() == "BLOSUM62":
                return cls.load_blosum62()
            elif matrix_name.upper() == "PAM250":
                return cls.load_pam250()
            raise


# Global singleton instances
BLOSUM62 = SubstitutionMatrix.load_blosum62()
PAM250 = SubstitutionMatrix.load_pam250()
