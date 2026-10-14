from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class AlignmentResult:
    """Encapsulates the complete result of a pairwise sequence alignment.

    Attributes:
        seq1_aligned: Aligned representation of query sequence (including gaps '-').
        seq2_aligned: Aligned representation of target sequence (including gaps '-').
        score: Numerical score computed by dynamic programming.
        start1: 0-indexed start position in seq1 (for local alignments; 0 for global).
        end1: 0-indexed end position (exclusive) in seq1 (len(seq1) for global).
        start2: 0-indexed start position in seq2 (for local alignments; 0 for global).
        end2: 0-indexed end position (exclusive) in seq2 (len(seq2) for global).
        alignment_type: Alignment algorithm type ('global' or 'local').
        gap_open: Linear gap penalty or affine opening penalty used.
        gap_extend: Gap extension penalty used.
        matrix_name: Name of the substitution matrix used (e.g., 'BLOSUM62', 'PAM250').
        seq1: Original unaligned query sequence.
        seq2: Original unaligned target sequence.
    """

    seq1_aligned: str
    seq2_aligned: str
    score: float
    start1: int
    end1: int
    start2: int
    end2: int
    alignment_type: str = "global"
    gap_open: float = 8.0
    gap_extend: float = 8.0
    matrix_name: str = "BLOSUM62"
    seq1: str = ""
    seq2: str = ""

    @property
    def length(self) -> int:
        """Returns the length of the aligned sequences (number of alignment columns)."""
        return len(self.seq1_aligned)

    @property
    def identity(self) -> float:
        """Computes percentage of identical aligned residues across all alignment columns."""
        if self.length == 0:
            return 0.0
        matches = sum(
            1
            for a, b in zip(self.seq1_aligned, self.seq2_aligned)
            if a == b and a != "-"
        )
        return (matches / self.length) * 100.0

    @property
    def similarity(self) -> float:
        """Computes fraction of aligned positions that do not contain a gap."""
        if self.length == 0:
            return 0.0
        aligned_pairs = sum(
            1
            for a, b in zip(self.seq1_aligned, self.seq2_aligned)
            if a != "-" and b != "-"
        )
        return (aligned_pairs / self.length) * 100.0
