# ProteinAlign Lab — Phase 5 Implementation and Research Validation Report

**Project**: ProteinAlign Lab: Mathematical Modelling and Comparative Evaluation of Protein Sequence Alignment  
**Milestone**: Phase 5 — Position-Specific Scoring Matrices (PSSMs) and Profile-Based Protein Alignment  
**Date**: October 14, 2026  
**Status**: 100% Complete & Verified (**127/127 Automated Tests Passing**)

---

## 1. Executive Summary

In Phase 5, we designed, implemented, and mathematically verified **position-specific scoring matrices (PSSMs)** and **profile-to-sequence dynamic programming alignment**:
1. **Pre-Implementation Audit**: Evaluated the repository's profile and scoring infrastructure, identifying 6 key limitations (such as lack of input validation, unvalidated background distributions, missing local/affine profile alignment, and absence of data-leakage controls), and resolved all of them.
2. **Profile Construction Layer** (`src/proteinalign/profiles.py`):
   - Created the immutable `Profile` abstraction capturing observed counts, regularized frequencies, gap fractions, consensus strings, and contributing sequence IDs.
   - Formalized an explicit gap policy (excluding gaps from amino acid count denominators $N_i$ while recording gap fractions $f_{\text{gap}, i}$) and an ambiguous residue distribution policy.
   - Built pseudocount regularization under Laplace/Dirichlet mixture priors with Robinson & Robinson (1991) and uniform background models.
3. **PSSM Log-Odds Scoring Engine** (`src/proteinalign/pssm.py`):
   - Implemented position-specific log-odds scoring:
     $$S_{i, a} = \lambda \cdot \log_b \left( \frac{P_{i, a}}{q_a} \right)$$
     with explicit, configurable scaling factors $\lambda$ and logarithm bases $b$ (e.g. half-bits).
   - Formulated mathematically unbiased ambiguous residue scoring ($X$ evaluates to the expected log-odds under background: $\mathbb{E}_q[S_{i, \cdot}]$).
4. **Profile-to-Sequence Dynamic Programming**:
   - Implemented both global (Needleman–Wunsch) and local (Smith–Waterman) profile alignment with linear and Gotoh 3-state affine gap penalties.
   - Built an independent score recomputation validator (`pssm.recompute_score(res)`) that verifies dynamic programming matrix scores against reconstructed alignments.
5. **Scientific Rigor & Data Leakage Prevention**:
   - Built `split_msa_train_test` for leave-one-out benchmark partitioning.
   - Implemented `validate_no_data_leakage` and `DataLeakageError` to prevent evaluating on sequences used in profile construction.
6. **Automated Testing Suite**: Added 33 new automated tests across `tests/test_profiles.py` and `tests/test_pssm.py`, bringing the total test suite to **127 passing tests**.

---

## 2. PART A — Pre-Implementation Audit Report

Before making changes for Phase 5, a comprehensive audit of the existing codebase was conducted. The table below outlines the identified issues and remediations:

| # | Issue | Severity | Location | Consequence | Recommended Fix | Applied Remediation |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **1** | **Missing Input Validation in Profile Extraction** | **High** | `profiles.py:17-26` | Invalid characters (digits, punctuation) or malformed inputs silently passed into column counts where they were silently dropped, distorting denominators. | Validate that sequences are non-empty, equal in length, and contain valid amino acids or gaps. | Added validation raising descriptive `ValueError` in `extract_alignment_columns`. |
| **2** | **Missing Profile Data Structure** | **High** | `profiles.py` | Profile frequencies and counts existed only as loose dictionaries, dropping metadata like sequence IDs, consensus, and gap fractions. | Encapsulate full profile metadata in a typed, frozen `Profile` dataclass. | Built `Profile` dataclass with accessors and metadata in `profiles.py`. |
| **3** | **Missing Local & Affine Profile Alignment** | **High** | `pssm.py` | Prototype only had basic linear global alignment; no local Smith–Waterman or Gotoh affine gap DP existed for profiles. | Implement full global/local linear/affine profile dynamic programming. | Implemented `align_global` and `align_local` supporting both linear and Gotoh affine gap models in `pssm.py`. |
| **4** | **Hardcoded Log Base & Fallback Score** | **Medium** | `pssm.py:35-73` | Logarithm base was hardcoded to `log2` and unknown residues received an arbitrary `-4.0` score. | Allow configurable scale ($\lambda$) and log base ($b$); compute expected background score for ambiguous residue $X$. | Implemented configurable parameters and expected background score calculation in `PSSM.score()`. |
| **5** | **Lack of Background Distribution Validation** | **Medium** | `profiles.py` | Raw `ROBINSON_FREQUENCIES` summed to 0.999 rather than 1.0, causing subtle frequency distortion; unvalidated user distributions could cause mathematical corruption. | Normalize `ROBINSON_FREQUENCIES` and add `validate_background_distribution()`. | Normalized Robinson frequencies to sum to 1.0; added strict positivity and sum validation. |
| **6** | **Absence of Data Leakage Safeguards** | **Medium** | `pssm.py`, `profiles.py` | No mechanism existed to prevent evaluating a profile against a sequence included in its training set, violating benchmark integrity. | Track sequence IDs in `Profile`/`PSSM`; add `split_msa_train_test` and `DataLeakageError`. | Implemented `split_msa_train_test`, `validate_no_data_leakage`, and `DataLeakageError`. |

---

## 3. Mathematical Formulation and Algorithmic Conventions

### 3.1 Observed Counts and Effective Depth
Given an MSA of $M$ aligned sequences of uniform length $L$:
- For column $i \in \{0, \dots, L-1\}$ and amino acid $a \in \Sigma$:
  $$c_{i, a} = \sum_{k=1}^M \mathbb{I}(s_{k, i} == a)$$
- Effective non-gap observation depth:
  $$N_i = \sum_{a \in \Sigma} c_{i, a}$$
- Gap fraction:
  $$f_{\text{gap}, i} = \frac{\sum_{k=1}^M \mathbb{I}(s_{k, i} \in \{'-', '.'\})}{M}$$

### 3.2 Pseudocount Regularization (Dirichlet / Laplace Prior)
$$P_{i, a} = \frac{c_{i, a} + \alpha \cdot q_a}{N_i + \alpha}$$
where:
- $\alpha \ge 0$ is the pseudocount weight (default $\alpha = 1.0$).
- $q_a > 0$ is the background proteome probability ($\sum_{a \in \Sigma} q_a = 1.0$).
- Sum property: $\sum_{a \in \Sigma} P_{i, a} = 1.0$ for all columns $i$.
- For gap-only columns ($N_i = 0$), $P_{i, a} = q_a$.

### 3.3 Log-Odds Scoring
$$S_{i, a} = \lambda \cdot \log_b \left( \frac{P_{i, a}}{q_a} \right)$$
- When $b = 2.0$ and $\lambda = 2.0$, scores are measured in **half-bits**.
- When $b = 2.0$ and $\lambda = 3.0$, scores are measured in **third-bits**.
- When $b = e$, scores are measured in **nats**.

### 3.4 Profile-to-Sequence Dynamic Programming
When aligning query sequence $Y = (y_1, \dots, y_n)$ against a PSSM profile of length $m$:
- The pairwise score $S(x_i, y_j)$ in Needleman–Wunsch or Smith–Waterman is substituted with the position-specific score:
  $$s(i, y_j) = S_{i-1, y_j}$$
- **Gotoh Affine Recurrence**:
  $$M(i, j) = S_{i-1, y_j} + \max\left( M(i-1, j-1), \, I_x(i-1, j-1), \, I_y(i-1, j-1) \right)$$
  $$I_x(i, j) = \max\left( M(i-1, j) - g_{\text{open}}, \, I_x(i-1, j) - g_{\text{extend}} \right) \quad \text{[Profile deletion / gap in query]}$$
  $$I_y(i, j) = \max\left( M(i, j-1) - g_{\text{open}}, \, I_y(i, j-1) - g_{\text{extend}} \right) \quad \text{[Query insertion / gap in profile]}$$

---

## 4. Worked Viva Example: Hand-Calculated Profile & PSSM

Consider an MSA with $M = 3$ sequences at column $i = 0$:
- Observed residues: `H`, `H`, `H` ($c_{0, H} = 3$, all other $c_{0, a} = 0$).
- Background: Uniform distribution $q_a = 0.05$ for all $a \in \Sigma$.
- Parameters: $\alpha = 1.0$, scale $\lambda = 2.0$, base $b = 2.0$.

### Hand Calculation:
1. **Effective depth**: $N_0 = 3$.
2. **Histidine probability**:
   $$P_{0, H} = \frac{3 + 1.0 \times 0.05}{3 + 1.0} = \frac{3.05}{4.0} = 0.7625$$
3. **Unobserved residue probability (e.g. Alanine $A$)**:
   $$P_{0, A} = \frac{0 + 1.0 \times 0.05}{3 + 1.0} = \frac{0.05}{4.0} = 0.0125$$
4. **Log-odds scores**:
   $$S_{0, H} = 2.0 \cdot \log_2 \left( \frac{0.7625}{0.05} \right) = 2.0 \cdot \log_2(15.25) = 2.0 \times 3.9307 = +7.861 \text{ half-bits}$$
   $$S_{0, A} = 2.0 \cdot \log_2 \left( \frac{0.0125}{0.05} \right) = 2.0 \cdot \log_2(0.25) = 2.0 \times (-2.0) = -4.000 \text{ half-bits}$$

Conserved Histidine receives a reward of $+7.86$ half-bits, while unobserved Alanine receives a penalty of $-4.00$ half-bits.

---

## 5. Demonstration & Score Recomputation

Running [`demo.py`](file:///d:/tcs/packages/sem9/MM/ProteinAlign-Lab/demo.py) section [5] demonstrates the complete workflow:
```text
[5] Position-Specific Scoring Matrix (PSSM) & Profile Alignment
----------------------------------------------------------------------
  [Note: The toy MSA below is illustrative and demonstrates mathematical mechanics]
Training MSA Homologs : ['homolog_1', 'homolog_2', 'homolog_3']
Held-Out Query Target : 'target_X' -> 'HEAPAWG'
Profile Consensus     : HEAGAWG (Length: 7)
PSSM Score for 'H' at pos 0 (conserved)  : +10.21 half-bits
PSSM Score for 'D' at pos 0 (unobserved) : -4.00 half-bits

Profile-to-Sequence Global Score: 52.67
  Profile : HEAGAWG
  Query   : HEAPAWG
  Independent Score Recomputation : 52.67 (Verified match: True)
  Data Leakage Safeguard: Successfully prevented evaluating training ID 'homolog_1'!
```

The independently recomputed score ($52.67$) matches the dynamic programming matrix score ($52.67$) to floating-point precision.

---

## 6. Automated Test Suite Status

The automated test suite now comprises **127 tests** across 9 test files:

```bash
pytest tests/ -v
```

```text
============================= test session starts =============================
platform win32 -- Python 3.12.2, pytest-8.3.5, pluggy-1.5.0
rootdir: D:\tcs\packages\sem9\MM\ProteinAlign-Lab
configfile: pyproject.toml
collected 127 items

tests\test_affine_gaps.py ....................                           [ 15%]
tests\test_alignment.py ............                                     [ 25%]
tests\test_biopython_compat.py ......                                    [ 29%]
tests\test_edge_cases.py ..............                                  [ 40%]
tests\test_global_alignment.py ...............                           [ 52%]
tests\test_local_alignment.py ..............                             [ 63%]
tests\test_profiles.py ..................                                [ 77%]
tests\test_pssm.py ...............                                       [ 89%]
tests\test_scoring.py .............                                      [100%]

============================= 127 passed in 0.28s =============================
```

---

## 7. Known Limitations & Phase 6 Readiness

### Known Limitations:
1. **Unweighted Sequences**: Profile construction currently uses unweighted sequence observations ($c_{i, a} = \sum \mathbb{I}$). For highly redundant MSAs, Henikoff position-based sequence weighting can further mitigate representation bias.
2. **Fixed Gap Penalties**: Profile alignment uses constant gap opening ($d$) and extension ($e$) penalties across all positions. Full profile HMMs (e.g. HMMER) allow position-specific gap transition probabilities.

### Phase 6 Readiness Checklist:
- [x] Needleman–Wunsch global alignment (linear & affine)
- [x] Smith–Waterman local alignment (linear & affine)
- [x] BLOSUM62 and PAM250 substitution matrices
- [x] Gotoh (1982) 3-state affine gap dynamic programming
- [x] Profile extraction and pseudocount regularization
- [x] PSSM generation and profile-to-sequence DP
- [x] Data leakage safeguards (`split_msa_train_test`, `validate_no_data_leakage`)
- [x] 127 passing automated unit and integration tests
- [ ] Phase 6: Large-scale benchmark experiments on MU-SABRE (303 reference alignments)
- [ ] Phase 6: Parameter sensitivity analysis, runtime/memory scaling, and final research report
