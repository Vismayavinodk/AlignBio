# ProteinAlign Lab — Testing Strategy and Verification Protocol

This document outlines the testing methodology, structural invariants, edge case taxonomy, and independent cross-validation protocols for the `proteinalign` toolkit.

---

## 1. Testing Philosophy and Standards

The testing strategy adheres to strict biomedical software engineering practices:
1. **Mathematical Verifiability**: Small toy examples must have analytical, hand-calculated dynamic programming scores to verify recurrence logic independently of code.
2. **Structural Invariants**: Every valid pairwise alignment must satisfy fundamental algebraic invariants regardless of sequence composition.
3. **Fail-Fast Parameter Validation**: Algorithms must reject invalid inputs (e.g. non-amino acid characters, negative penalties) at the boundary rather than producing corrupted matrices.
4. **Independent Validation**: Alignment scores must be cross-validated against an independent, trusted reference implementation (**Biopython**).
5. **Separation of Concerns**: Unit and property tests are kept strictly separate from long-running benchmark experiments.

---

## 2. Invariant Specifications

Every alignment $\mathcal{A}(X, Y) = \begin{pmatrix} X' \\ Y' \end{pmatrix}$ produced by any algorithm must satisfy:

### Invariant 1: Degap Sequence Recoverability
Removing all gap characters `'-'` must reconstruct the original input sequences exactly:
$$\pi_{\Sigma}(X') = X \quad \text{and} \quad \pi_{\Sigma}(Y') = Y$$

### Invariant 2: Equal Length Matrix Invariant
The aligned strings must have identical lengths:
$$|X'| = |Y'| = L, \quad \text{where } \max(|X|, |Y|) \le L \le |X| + |Y|$$

### Invariant 3: Score Consistency Invariant
The dynamic programming score returned by the algorithm must equal the manual summation of substitution scores and gap penalties computed across columns:
$$\text{Score} = \sum_{k=1}^L \sigma(x'_k, y'_k)$$

### Invariant 4: Symmetry Invariant (for Pairwise Matrices)
For any symmetric substitution matrix $S(a, b) = S(b, a)$:
$$\text{score}(X, Y) = \text{score}(Y, X)$$

---

## 3. Test Suites Taxonomy

### 3.1 Unit & Hand-Calculated Tests (`tests/test_alignment.py`)
- **Self-Alignment Identity**: Aligning sequence $X = \text{"HEAGAWGHEE"}$ to itself must produce 100% sequence identity and a score equal to the sum of self-substitution scores in BLOSUM62:
  $$8(\text{H}) + 5(\text{E}) + 4(\text{A}) + 6(\text{G}) + 4(\text{A}) + 11(\text{W}) + 6(\text{G}) + 8(\text{H}) + 5(\text{E}) + 5(\text{E}) = 62.0$$
- **Single Substitution**: `"HEAG"` vs `"HEAP"` aligns with zero gaps, yielding:
  $$8(\text{H-H}) + 5(\text{E-E}) + 4(\text{A-A}) - 2(\text{G-P}) = 15.0$$
- **Single Deletion**: `"HEA"` vs `"HEAG"` under linear gap penalty $d=8$:
  $$8 + 5 + 4 - 8 = 9.0$$
- **Local Motif Isolation**: An identical motif `"PAWHEAE"` embedded within dissimilar flanking sequences must be extracted with exact motif coordinates and self-score $44.0$.
- **Affine Gap Length Scaling**: Gotoh DP penalizes a 2-residue gap under $d=10, e=1$ as $-(10 + 1) = -11.0$, preferring a single contiguous block over multiple isolated gaps.

### 3.2 Edge Cases and Exception Handling (`tests/test_edge_cases.py`)
- **Empty Sequences**:
  - Global alignment of `("", "")` $\to$ score 0.0, output `("", "")`.
  - Global alignment of `("HEA", "")` $\to$ score $-3d$, output `("HEA", "---")`.
  - Local alignment with any empty string $\to$ score 0.0.
- **Single-Residue Sequences**: Match (`"H"`, `"H"`), mismatch (`"H"`, `"P"`).
- **Homopolymeric Repeats**: Aligning repetitive tracts (`"AAAA"`, `"AA"`).
- **Extreme Length Discrepancy**: 1 residue vs 20 residues.
- **Completely Dissimilar Local Alignment**: Strings with only strongly negative BLOSUM pairs (e.g. `"WWWW"` vs `"DDDD"`) yield score 0.0 and empty string.
- **Parameter Validation**:
  - Invalid characters (e.g. `"HEAG123"`, `"HEA$G"`) raise `ValueError`.
  - Negative gap penalties raise `ValueError`.
  - Inconsistent affine parameters ($gap\_open < gap\_extend$) raise `ValueError`.

### 3.3 Affine Gap Dynamic Programming Tests (`tests/test_affine_gaps.py`)
- **Mathematical Gap Cost Scaling**:
  - $G(1) = g_{\text{open}} = 10.0$
  - $G(2) = g_{\text{open}} + g_{\text{extend}} = 11.0$
  - $G(3) = g_{\text{open}} + 2 \times g_{\text{extend}} = 12.0$
  - Custom penalty scaling ($g_{\text{open}} = 15.0, g_{\text{extend}} = 2.0$)
- **Contiguous Indel Grouping**:
  - Confirms that Gotoh DP groups indels into contiguous blocks (e.g. `----`) rather than multiple scattered single-residue indels (`--`, `--`), reflecting biological loop mutations.
- **Local Affine Alignment**:
  - Motif extraction from divergent flanking regions.
  - Coordinate slice invariant: `seq[start:end] == seq_aligned.replace('-', '')`.
  - All-negative sequences return score 0.0, empty alignment strings, and zero coordinates.
- **Parameter Validation & Edge Cases**:
  - Rejection of invalid non-IUPAC characters with `ValueError`.
  - Strict enforcement of non-negative penalties ($g_{\text{open}} \ge 0, g_{\text{extend}} \ge 0$).
  - Enforcement of $g_{\text{open}} \ge g_{\text{extend}}$.
  - Retention of input sequences: `res.seq1 == seq1` and `res.seq2 == seq2`.

### 3.4 Independent Cross-Validation (`tests/test_biopython_compat.py`)
Alignment results are validated against `Bio.Align.PairwiseAligner`:
- **Linear Global & Local Alignment**: Score matches Biopython within $10^{-4}$ under BLOSUM62.
- **Affine Global Alignment**: Score matches Biopython within $10^{-4}$ under BLOSUM62 ($g_{\text{open}} = 10.0, g_{\text{extend}} = 1.0$).
- **Affine Global Alignment under PAM250**: Score matches Biopython within $10^{-4}$ ($g_{\text{open}} = 11.0, g_{\text{extend}} = 1.0$).
- **Affine Local Alignment**: Score matches Biopython within $10^{-4}$ for motif isolation and internal gap insertion.

### 3.5 Profile Construction and Regularization Tests (`tests/test_profiles.py`)
- **Alignment Transpose & Validation**:
  - Rejection of empty alignments, length mismatches, and invalid non-amino acid characters.
- **Observed Column Counts & Gap Policy**:
  - Gaps (`'-'`, `'.'`) are excluded from amino acid count denominators ($N_i$).
  - Ambiguous residues (`'B'`, `'Z'`, `'X'`) are distributed fractionally according to documented policy.
- **Pseudocount Regularization & Normalization**:
  - Strict probability normalization ($\sum_a P_{i, a} = 1.0$).
  - Analytical verification of regularized probabilities: $P_{i, a} = (c_{i, a} + \alpha q_a) / (N_i + \alpha)$.
  - All-gap columns default smoothly to the background prior $q_a$.
- **Data Leakage Safeguards**:
  - `split_msa_train_test`: Validates held-out target extraction and un-gapping.
  - `validate_no_data_leakage`: Detects query ID overlap and raises `DataLeakageError`.

### 3.6 PSSM Log-Odds Scoring and Profile Alignment Tests (`tests/test_pssm.py`)
- **Analytical Log-Odds Score Verification**:
  - Hand-calculated scores under uniform and Robinson background distributions.
  - Ambiguous residue scoring (`B`, `Z`, `X`).
  - Score scaling ($\lambda = 1.0, 2.0, 3.0$) and logarithm bases ($2.0, e, 10.0$).
- **Profile-to-Sequence Global Alignment**:
  - Global linear and Gotoh 3-state affine dynamic programming.
  - Reconstruction invariants and length constraints.
- **Profile-to-Sequence Local Alignment**:
  - Smith–Waterman motif isolation from divergent flanking sequences.
  - Affine internal indel handling and coordinate slice matching.
- **Score Recomputation Invariant**:
  - `pssm.recompute_score(res) == res.score` directly verifies DP matrix scores against reconstructed alignments.
- **Data Leakage Prevention**:
  - Alignment engine raises `DataLeakageError` when query sequence ID belongs to training MSA without `allow_leakage=True`.

---

## 4. Test Execution Instructions

Run the complete test suite (127 passing tests):
```bash
pytest tests/ -v
```

Run only profile construction tests:
```bash
pytest tests/test_profiles.py -v
```

Run only PSSM scoring and profile alignment tests:
```bash
pytest tests/test_pssm.py -v
```

Run Biopython independent cross-validation:
```bash
pytest tests/test_biopython_compat.py -v
```
