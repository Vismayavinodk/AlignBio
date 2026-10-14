# ProteinAlign Lab

**Mathematical Modelling and Comparative Evaluation of Protein Sequence Alignment Using Dynamic Programming and Position-Specific Scoring Matrices.**

## Overview
ProteinAlign Lab is an extensible, modular Python research toolkit designed for rigorous comparative analysis of protein sequence alignment algorithms, substitution scoring schemes, gap penalty formulations, and position-specific scoring models.

Algorithms are benchmarked against **MU-SABRE**, a multiple sequence alignment benchmark derived from **SABmark v1.65** re-aligned using the **MUSTANG** structural superposition algorithm.

---

## Key Modules and Features
- **Global Alignment**: Needleman–Wunsch algorithm with linear gap penalties (`global_alignment.py`).
- **Local Alignment**: Smith–Waterman algorithm with zero-floor cutoff (`local_alignment.py`).
- **Affine Gap Penalties**: Gotoh (1982) three-state dynamic programming decomposing gap costs into opening ($d$) and extension ($e$) (`affine_gaps.py`).
- **Scoring Matrices**: Standard BLOSUM62 and PAM250 log-odds matrices (`scoring.py`).
- **Profile & PSSM Modeling**: Column count extraction, pseudocount regularization, and profile-to-sequence dynamic programming (`profiles.py`, `pssm.py`).
- **Benchmark Evaluation**: Reference core column extraction ($Q$-score/developer score), Sum-of-Pairs (SP-score), and percent identity (`evaluation.py`).
- **Independent Validation**: Cross-validated against Biopython (`Bio.Align.PairwiseAligner`).

---

## Documentation
- [`docs/mathematics.md`](docs/mathematics.md): Complete mathematical derivations, recurrence relations, and assumptions.
- [`docs/methodology.md`](docs/methodology.md): Benchmark provenance, pairwise decomposition, train/val/test splits, and data leakage safeguards.
- [`docs/architecture.md`](docs/architecture.md): Software design, module responsibilities, function signatures, and interface contracts.
- [`docs/testing.md`](docs/testing.md): Testing strategy, algebraic invariants, and edge case specifications.
- [`docs/phase3_plan.md`](docs/phase3_plan.md): Phase 3 tasks, acceptance criteria, milestones, and dependencies.
- [`report/dataset_inspection_report.md`](report/dataset_inspection_report.md): Quantitative audit of the 303 SABRE benchmark files.

---

## Quick Start

### Installation
Clone the repository and install in editable mode:
```bash
git clone <repo-url>
cd ProteinAlign-Lab
pip install -e .
```

### Running Tests
Execute the comprehensive test suite:
```bash
python -m pytest tests/ -v
```

### Basic Usage Example
```python
from proteinalign import needleman_wunsch, needleman_wunsch_affine, BLOSUM62, PAM250

seq1 = "HEAGAWGHEE"
seq2 = "PAWHEAE"

# Global alignment with linear gap penalty
result_linear = needleman_wunsch(seq1, seq2, matrix=BLOSUM62, gap_penalty=8.0)
print(f"Linear NW Score: {result_linear.score}")
print(result_linear.seq1_aligned)
print(result_linear.seq2_aligned)

# Global alignment with Gotoh affine gap penalties
result_affine = needleman_wunsch_affine(seq1, seq2, matrix=BLOSUM62, gap_open=10.0, gap_extend=1.0)
print(f"Affine Gotoh Score: {result_affine.score}")
```

---

## License
MIT License. See [LICENSE](LICENSE) for details.
