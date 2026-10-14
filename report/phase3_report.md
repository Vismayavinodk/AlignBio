# ProteinAlign Lab — Phase 3 Implementation and Testing Report

**Project**: ProteinAlign Lab: Mathematical Modelling and Comparative Evaluation of Protein Sequence Alignment  
**Milestone**: Phase 3 — Core Linear-Gap Global & Local Alignment Implementation  
**Date**: October 14, 2026  
**Status**: 100% Complete & Verified (70/70 Automated Tests Passing)

---

## 1. Executive Summary

In Phase 3, we implemented and rigorously validated the core dynamic programming algorithms for pairwise protein sequence alignment under a linear gap penalty model:
1. **Needleman–Wunsch Global Alignment** (`src/proteinalign/global_alignment.py`).
2. **Smith–Waterman Local Alignment** (`src/proteinalign/local_alignment.py`).
3. **Pluggable Substitution Scoring Layer** (`src/proteinalign/scoring.py`), featuring standard BLOSUM62, PAM250, custom identity/mismatch matrices, and ambiguous residue policies.
4. **Comprehensive Test Suite** (`tests/test_global_alignment.py`, `tests/test_local_alignment.py`, `tests/test_scoring.py`, `tests/test_edge_cases.py`, `tests/test_biopython_compat.py`).
5. **Cross-Validation with Biopython**: Verified that our independent implementation yields exact mathematical score equivalence with Biopython's `Bio.Align.PairwiseAligner`.

---

## 2. Mathematical Mapping to Code

### 2.1 Needleman–Wunsch (Global Alignment)
- **Recurrence**:
  $$F(i, j) = \max\left( F(i-1, j-1) + S(x_i, y_j), \, F(i-1, j) - d, \, F(i, j-1) - d \right)$$
- **Code Mapping**:
  ```python
  diag = F[i - 1][j - 1] + matrix.score(char1, char2)
  up   = F[i - 1][j] - gap_penalty
  left = F[i][j - 1] - gap_penalty
  F[i][j] = max(diag, up, left)
  ```
- **Boundary Initialization**:
  $F(0, 0) = 0$, $F(i, 0) = -i \cdot d$, $F(0, j) = -j \cdot d$.
- **Complexity**: Time $\mathcal{O}(m \cdot n)$, Space $\mathcal{O}(m \cdot n)$.

### 2.2 Smith–Waterman (Local Alignment)
- **Recurrence**:
  $$H(i, j) = \max\left( 0, \, H(i-1, j-1) + S(x_i, y_j), \, H(i-1, j) - d, \, H(i, j-1) - d \right)$$
- **Code Mapping**:
  ```python
  diag = H[i - 1][j - 1] + matrix.score(char1, char2)
  up   = H[i - 1][j] - gap_penalty
  left = H[i][j - 1] - gap_penalty
  score = max(0.0, diag, up, left)
  H[i][j] = score
  ```
- **Boundary Initialization**: $H(i, 0) = 0$, $H(0, j) = 0$.
- **Traceback**: Starts at $\arg\max_{i, j} H(i, j)$ and terminates immediately when $H(i, j) = 0$.

---

## 3. Worked Traceback Example (Viva Preparation)

Consider a small toy example:
- Query $X = \text{"HEAG"}$, Target $Y = \text{"HEAP"}$
- Scoring: Simple Identity Matrix ($+2$ match, $-1$ mismatch), Linear Gap $d = 3.0$.

### Forward DP Matrix $F(i, j)$:

| | $\epsilon$ | H | E | A | P |
|---|---|---|---|---|---|
| **$\epsilon$** | **0.0** | -3.0 | -6.0 | -9.0 | -12.0 |
| **H** | -3.0 | **2.0** | -1.0 | -4.0 | -7.0 |
| **E** | -6.0 | -1.0 | **4.0** | 1.0 | -2.0 |
| **A** | -9.0 | -4.0 | 1.0 | **6.0** | 3.0 |
| **G** | -12.0 | -7.0 | -2.0 | 3.0 | **5.0** |

### Step-by-Step Traceback:
1. Start at $(4, 4)$ where $F(4, 4) = 5.0$:
   - Diagonal predecessor $F(3, 3) + S(\text{G}, \text{P}) = 6.0 + (-1.0) = 5.0$.
   - Equality holds! Emit column $(\text{G}, \text{P})$. Step to $(3, 3)$.
2. At $(3, 3)$ where $F(3, 3) = 6.0$:
   - Diagonal predecessor $F(2, 2) + S(\text{A}, \text{A}) = 4.0 + 2.0 = 6.0$.
   - Emit column $(\text{A}, \text{A})$. Step to $(2, 2)$.
3. At $(2, 2)$ where $F(2, 2) = 4.0$:
   - Diagonal predecessor $F(1, 1) + S(\text{E}, \text{E}) = 2.0 + 2.0 = 4.0$.
   - Emit column $(\text{E}, \text{E})$. Step to $(1, 1)$.
4. At $(1, 1)$ where $F(1, 1) = 2.0$:
   - Diagonal predecessor $F(0, 0) + S(\text{H}, \text{H}) = 0.0 + 2.0 = 2.0$.
   - Emit column $(\text{H}, \text{H})$. Step to $(0, 0)$.
5. Terminate at $(0, 0)$. Reverse emitted columns:
   ```text
   Seq1: H E A G
   Seq2: H E A P
   Score: 5.0
   ```

---

## 4. Test Suite Execution Results

Running `pytest` executed **70 automated tests** across 6 test suites:

```text
tests/test_global_alignment.py ..................                         [ 25%]
tests/test_local_alignment.py ..................                          [ 51%]
tests/test_scoring.py ................                                    [ 74%]
tests/test_edge_cases.py ..............                                   [ 94%]
tests/test_biopython_compat.py ..                                         [ 97%]
tests/test_alignment.py ..                                                [100%]

============================== 70 passed in 0.44s ==============================
```

### Coverage Breakdown
- **Global Alignment (18 tests)**: Identical pairs, single mismatches, internal indels, homopolymers, variable gap penalties, parameter sensitivity, score recomputation consistency, and degap recovery.
- **Local Alignment (18 tests)**: Strong matching subsequences, exclusion of negative flanks, completely dissimilar strings (score 0), single-residue motifs, start/end slice indexing.
- **Scoring (16 tests)**: Match/mismatch matrices, BLOSUM62/PAM250 log-odds values, symmetry property, ambiguous residues (`X`, `B`, `Z`), residue validation exceptions.
- **Edge Cases & Invariants (14 tests)**: Zero-length inputs, single-residue inputs, extreme length discrepancies (1 vs 20), input validation exceptions, non-negative gap constraints.
- **Biopython Cross-Validation (2 tests)**: Verified that global and local alignment scores agree with `Bio.Align.PairwiseAligner` within $10^{-4}$.

---

## 5. Acceptance Criteria Checklist

- [x] Needleman–Wunsch returns a correct optimal global alignment and score.
- [x] Smith–Waterman returns a correct optimal local alignment and score.
- [x] Both algorithms reconstruct valid alignments through traceback.
- [x] Both algorithms accept configurable scoring and linear gap penalties.
- [x] BLOSUM62 and PAM250 work through the common scoring interface.
- [x] Edge cases and invalid inputs are handled with clear exceptions.
- [x] Unit tests cover ordinary cases, boundary cases, and traceback correctness.
- [x] Scores are consistent with returned alignments.
- [x] Independent comparisons with Biopython conducted and validated.
- [x] Package installs in editable mode and imports cleanly.
- [x] Usage documentation reflects the actual implementation.
- [x] Repository is clean, organized, and ready for Phase 4.

---

## 6. Transition to Phase 4

With Phase 3 complete and verified, the foundation is set for **Phase 4: Affine Gap Penalties and Gotoh's Three-State Dynamic Programming**. In Phase 4, we will:
1. Deepen the Gotoh 3-state DP engine ($M, I_x, I_y$).
2. Analyze biological consequences of separating gap opening ($d$) from gap extension ($e$).
3. Compare gap distribution patterns (contiguous blocks vs scattered single indels).
