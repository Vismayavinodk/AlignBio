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
- [`docs/independent_review.md`](docs/independent_review.md): Phase 8 external audit and independent verification report.
- [`docs/submission_readiness.md`](docs/submission_readiness.md): Final readiness assessment across 11 submission categories.
- [`report/project_defense.md`](report/project_defense.md): Professor-facing viva defense guide answering 15 core questions.
- [`report/phase6_report.md`](report/phase6_report.md): Comprehensive empirical evaluation, statistical findings, and scaling analysis.
- [`docs/architecture.md`](docs/architecture.md): Software design, module responsibilities, function signatures, and interface contracts.
- [`docs/testing.md`](docs/testing.md): Testing strategy, algebraic invariants, and edge case specifications.

---

## Empirical Benchmark Findings Summary

All experiments evaluated on the **MU-SABRE** structural benchmark (303 reference alignments, $N=288$ paired comparisons, $N=183$ held-out profile targets):

| Experiment | Comparison | Key Metric | Result | Statistical Significance |
| :--- | :--- | :--- | :--- | :--- |
| **RQ1 (Gap Models)** | Linear ($d=8$) vs Gotoh Affine ($d=10, e=1$) | Developer $Q$-score | Affine: **0.6150** vs Linear: **0.4082** | Wilcoxon **$p = 7.19 \times 10^{-16}$** |
| **RQ1 (Indels)** | Linear vs Gotoh Affine | Indel Openings / Mean Length | Openings: **8.73** vs **20.06** / Length: **5.17** vs **1.56** aa | Wilcoxon **$p = 3.96 \times 10^{-45}$** |
| **RQ2 (Matrices)** | BLOSUM62 vs PAM250 (15–25% ID) | Developer $Q$-score | BLOSUM62: **0.5000** vs PAM250: **0.3750** | Wilcoxon **$p = 8.28 \times 10^{-4}$** |
| **RQ3 (Modes)** | Global vs Local under length asymmetry | Query Coverage | Local isolates core: **22.2%** vs Global: **100%** | Local prevents terminal gap penalty |
| **RQ4 (PSSM LOO)** | PSSM Profile vs Fixed BLOSUM62 | Core Placement $Q$-score | PSSM: **0.6000** (Mean 0.5349, zero leakage) | $p = 1.97 \times 10^{-5}$ |
| **RQ6 (Scaling)** | DP Runtime vs Length Product $m \cdot n$ | Power-Law Exponent $\beta$ | Affine: **$\beta = 1.096$** ($R^2=0.9914$) | Matches theoretical $\mathcal{O}(mn)$ |
| **RQ6 (Memory)** | Peak Heap Memory Allocation | Memory Ratio (Affine / Linear) | Ratio = **$3.00\times$** ($L=50$ to $L=1000$) | Exactly reflects 3 DP matrices |

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
Execute the comprehensive test suite (148 tests):
```bash
pytest -v
```

### Reproducing Benchmark Experiments
```bash
# 1. Dataset validation and stratified train/val/test split
python experiments/validate_dataset.py

# 2. Rapid smoke test demonstration (< 2 seconds)
python experiments/run_benchmarks.py --smoke

# 3. Full test split benchmark evaluation (RQ1 - RQ4)
python experiments/run_benchmarks.py --config experiments/configs/benchmark_config.json

# 4. Parameter sensitivity grid search (RQ5)
python experiments/run_sensitivity.py --config experiments/configs/sensitivity_config.json

# 5. Computational runtime and memory scaling benchmarks (RQ6)
python experiments/run_scaling.py --config experiments/configs/scaling_config.json

# 6. Statistical tests and publication figure generation
python experiments/analyze_results.py
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

# PSSM Profile-to-Sequence Alignment
from proteinalign import build_profile_from_msa, PSSM

msa = {"s1": "HEAGAWG", "s2": "HEAPAWG"}
profile = build_profile_from_msa(msa, name="my_fam", pseudocount_weight=1.0)
pssm = PSSM.from_profile(profile, scale=2.0)
res_pssm = pssm.align_global("HEAGAWA", gap_open=10.0, gap_extend=1.0)
print(f"PSSM Profile Score: {res_pssm.score}")
```

---

## License
MIT License. See [LICENSE](LICENSE) for details.
