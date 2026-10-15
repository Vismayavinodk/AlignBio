# ProteinAlign Lab — Experimental Evaluation & Statistical Summary

**Phase 6: Rigorous Empirical Evaluation of Dynamic Programming Alignment & Scoring Models**

## RQ1: Linear vs Gotoh Affine Gap Penalties
- **Sample Size:** 288 paired sequence comparisons
- **Q-score (Accuracy):** Linear Median = 0.4082 (IQR=0.8589) vs Affine Median = 0.6150 (IQR=0.7720) | **Wilcoxon p = 7.19e-16**
- **Indel Fragmentation (Gap Openings):** Linear Mean = 20.1 vs Affine Mean = 8.7 | **Wilcoxon p = 3.96e-45**
- **Mean Indel Length:** Linear Mean = 1.56 residues vs Affine Mean = 5.17 residues | **Wilcoxon p = 1.14e-47**
- **Computational Cost:** Linear Median = 15.07 ms vs Affine Median = 24.60 ms (approx 1.5–2.0x overhead for 3-layer DP)

## RQ2: Substitution Matrices (BLOSUM62 vs PAM250)
| Homology Tier | N Pairs | BLOSUM62 Median Q (IQR) | PAM250 Median Q (IQR) | Wilcoxon p-value |
| :--- | :--- | :--- | :--- | :--- |
| Twilight (<15%) | 36 | 0.0000 (0.1679) | 0.0000 (0.4516) | 9.53e-01 |
| Low-Mid (15-25%) | 111 | 0.5000 (0.5902) | 0.3750 (0.7044) | 8.28e-04 |
| High (>25%) | 141 | 0.8750 (0.4318) | 0.8654 (0.5000) | 4.92e-03 |

## RQ3: Alignment Mode (Global Needleman-Wunsch vs Local Smith-Waterman)
- **High-Disparity Pairs (|L1-L2|/max > 0.25):** N = 52
  - Local Alignment Query Coverage Median: 22.2% (isolates conserved core domain, avoiding penalty for non-homologous extensions)
  - Accuracy: Global Median Q = 0.3673 vs Local Median Q = 0.0000
- **Low-Disparity Pairs (<= 0.25):** N = 236 | Global Median Q = 0.7054 vs Local Median Q = 0.0000

## RQ4: PSSM Profile vs Fixed Matrix Baseline (Leave-One-Out)
- **Held-Out Test Sequences:** N = 183 (strictly zero data leakage)
- **PSSM Profile Placement Q-score:** Median = 0.6000, Mean = 0.5349 (IQR=0.7716)
- **Fixed BLOSUM62 Consensus Q-score:** Median = 0.6932, Mean = 0.5990 (IQR=0.6674)
- **Statistical Test:** Wilcoxon signed-rank p = 1.97e-05

## RQ6: Computational Efficiency and Asymptotic Scaling
| Dynamic Programming Algorithm | Empirical Scaling Exponent (beta) | R-squared | Theoretical Complexity |
| :--- | :--- | :--- | :--- |
| Gotoh Global Affine | 1.096 | 0.9914 | O(m * n) (beta = 1.0) |
| Gotoh Local Affine | 1.094 | 0.9912 | O(m * n) (beta = 1.0) |
| Needleman-Wunsch Linear | 1.083 | 0.9907 | O(m * n) (beta = 1.0) |
| PSSM Global Affine | 1.105 | 0.9922 | O(m * n) (beta = 1.0) |
| Smith-Waterman Linear | 1.084 | 0.9913 | O(m * n) (beta = 1.0) |

*Note: All scaling experiments conducted on Windows 64-bit Python 3.12 with wall-clock timing via time.perf_counter_ns and heap memory tracked via tracemalloc.*