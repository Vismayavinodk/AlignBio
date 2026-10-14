"""ProteinAlign: Pairwise and Multiple Sequence Alignment Benchmark Toolkit."""

from .models import AlignmentResult
from .scoring import BLOSUM62, PAM250, SubstitutionMatrix
from .global_alignment import needleman_wunsch
from .local_alignment import smith_waterman
from .affine_gaps import needleman_wunsch_affine, smith_waterman_affine
from .profiles import (
    ROBINSON_FREQUENCIES,
    UNIFORM_FREQUENCIES,
    extract_alignment_columns,
    compute_column_counts,
    compute_profile_frequencies,
)
from .pssm import PSSM
from .io import read_fasta, write_fasta, read_msa
from .evaluation import compute_pid, compute_q_score, compute_sp_score
from .visualization import format_alignment, alignment_to_html

__version__ = "0.1.0"

__all__ = [
    "AlignmentResult",
    "SubstitutionMatrix",
    "BLOSUM62",
    "PAM250",
    "needleman_wunsch",
    "smith_waterman",
    "needleman_wunsch_affine",
    "smith_waterman_affine",
    "ROBINSON_FREQUENCIES",
    "UNIFORM_FREQUENCIES",
    "extract_alignment_columns",
    "compute_column_counts",
    "compute_profile_frequencies",
    "PSSM",
    "read_fasta",
    "write_fasta",
    "read_msa",
    "compute_pid",
    "compute_q_score",
    "compute_sp_score",
    "format_alignment",
    "alignment_to_html",
]
