"""Independent validation against Biopython's PairwiseAligner."""

import pytest

try:
    from Bio import Align
    from Bio.Align import substitution_matrices
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False

from proteinalign import (
    BLOSUM62,
    PAM250,
    needleman_wunsch,
    smith_waterman,
    needleman_wunsch_affine,
    smith_waterman_affine,
)


@pytest.mark.skipif(not BIOPYTHON_AVAILABLE, reason="Biopython is not installed")
class TestBiopythonValidation:
    """Compares dynamic programming alignment scores against Biopython reference engine."""

    def test_global_linear_score_matches_biopython(self):
        s1 = "HEAGAWGHEE"
        s2 = "PAWHEAE"

        aligner = Align.PairwiseAligner()
        aligner.mode = "global"
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
        aligner.open_gap_score = -8.0
        aligner.extend_gap_score = -8.0

        expected_score = aligner.score(s1, s2)
        my_result = needleman_wunsch(s1, s2, matrix=BLOSUM62, gap_penalty=8.0)

        assert pytest.approx(my_result.score, 1e-4) == expected_score

    def test_local_linear_score_matches_biopython(self):
        s1 = "DDDDPAWHEAERRRR"
        s2 = "WWWWPAWHEAEWWWW"

        aligner = Align.PairwiseAligner()
        aligner.mode = "local"
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
        aligner.open_gap_score = -8.0
        aligner.extend_gap_score = -8.0

        expected_score = aligner.score(s1, s2)
        my_result = smith_waterman(s1, s2, matrix=BLOSUM62, gap_penalty=8.0)

        assert pytest.approx(my_result.score, 1e-4) == expected_score

    def test_global_affine_score_matches_biopython(self):
        s1 = "HEAGAWGHEE"
        s2 = "PAWHEAE"

        aligner = Align.PairwiseAligner()
        aligner.mode = "global"
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
        aligner.open_gap_score = -10.0
        aligner.extend_gap_score = -1.0

        expected_score = aligner.score(s1, s2)
        my_result = needleman_wunsch_affine(s1, s2, matrix=BLOSUM62, gap_open=10.0, gap_extend=1.0)

        assert pytest.approx(my_result.score, 1e-4) == expected_score

    def test_global_affine_pam250_matches_biopython(self):
        s1 = "MKFLILLFNILCLFPVLAAD"
        s2 = "MKFLILLVLFPAAD"

        aligner = Align.PairwiseAligner()
        aligner.mode = "global"
        aligner.substitution_matrix = substitution_matrices.load("PAM250")
        aligner.open_gap_score = -11.0
        aligner.extend_gap_score = -1.0

        expected_score = aligner.score(s1, s2)
        my_result = needleman_wunsch_affine(s1, s2, matrix=PAM250, gap_open=11.0, gap_extend=1.0)

        assert pytest.approx(my_result.score, 1e-4) == expected_score

    def test_local_affine_score_matches_biopython(self):
        s1 = "DDDDPAWHEAERRRR"
        s2 = "WWWWPAWHEAEWWWW"

        aligner = Align.PairwiseAligner()
        aligner.mode = "local"
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
        aligner.open_gap_score = -10.0
        aligner.extend_gap_score = -1.0

        expected_score = aligner.score(s1, s2)
        my_result = smith_waterman_affine(s1, s2, matrix=BLOSUM62, gap_open=10.0, gap_extend=1.0)

        assert pytest.approx(my_result.score, 1e-4) == expected_score

    def test_local_affine_with_internal_gap_matches_biopython(self):
        s1 = "XXXXHEAGWWWWWWXXXX"
        s2 = "YYYYHEAWAGWWWWWWYYYY"

        aligner = Align.PairwiseAligner()
        aligner.mode = "local"
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
        aligner.open_gap_score = -4.0
        aligner.extend_gap_score = -1.0

        expected_score = aligner.score(s1, s2)
        my_result = smith_waterman_affine(s1, s2, matrix=BLOSUM62, gap_open=4.0, gap_extend=1.0)

        assert pytest.approx(my_result.score, 1e-4) == expected_score
