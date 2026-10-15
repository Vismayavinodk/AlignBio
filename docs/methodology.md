# ProteinAlign Lab — Experimental Methodology and Benchmark Design

This document details the dataset provenance, pairwise experiment generation, benchmark validation protocol, data leakage prevention, and experimental evaluation methodology.

---

## 1. Dataset Provenance and Citation

### 1.1 Benchmark Origin: MU-SABRE
The benchmark dataset used in ProteinAlign Lab is located at `data/raw/bench1.0/sabrem/ref`. 

As documented in `data/raw/bench1.0/BENCH_README.txt`, this dataset was assembled by **Dr. Robert C. Edgar** (author of MUSCLE) and consists of the **SABRE** multiple sequence alignments derived from **SABmark v1.65**, re-aligned using the **MUSTANG** structural alignment tool:
- **Primary Source**: SABmark (Sequence Alignment Benchmark) version 1.65.
- **Structural Re-alignment**: MUSTANG (MUltiple STructural Alignment AlGorithm).
- **Benchmark Name**: `MU-SABRE` (`sabrem`).

### 1.2 Academic Citations
1. **SABmark Benchmark**:
   > Van Walle, I., Lasters, I., & Wyns, L. (2005). SABmark—a benchmark for sequence alignment that covers the entire known fold space. *Bioinformatics*, 21(7), 1267–1268.
2. **MUSTANG Structural Alignment**:
   > Konagurthu, A. S., Whisstock, J. C., Henschel, A., & Lesk, A. M. (2006). MUSTANG: A multiple structural alignment algorithm. *Proteins: Structure, Function, and Bioinformatics*, 64(3), 559–574.
3. **Benchmark Compilation and MUSCLE Evaluation**:
   > Edgar, R. C. (2010). Quality measures for protein alignment benchmarks. *Nucleic Acids Research*, 38(7), 2145–2153.
   > Edgar, R. C. (2004). MUSCLE: multiple sequence alignment with high accuracy and high throughput. *Nucleic Acids Research*, 32(5), 1792–1797.

### 1.3 Benchmark File Format & Core Conservation Convention
- **File Format**: Standard FASTA text files.
- **Case Convention**:
  - **Uppercase characters** (`A-Z`): Represent conserved structural core residues aligned with high geometric confidence by MUSTANG based on known 3D PDB coordinates. These residues constitute the ground-truth standard for accuracy assessment.
  - **Lowercase characters** (`a-z`): Represent loops, unconserved termini, and structurally variable segments. By benchmark convention, lowercase residues are excluded from core assessment.
  - **Hyphens (`-`)**: Represent structural insertions and deletions.

---

## 2. Converting Reference MSAs to Pairwise Experimental Data

Each file in `sabrem/ref` is a Multiple Sequence Alignment containing $K$ homologous sequences ($3 \le K \le 24$). For pairwise evaluation:

### 2.1 Pairwise Extraction
For an MSA containing sequences $\{S_1, S_2, \dots, S_K\}$:
1. Generate all unique unordered pairs $(S_i, S_j)$ where $1 \le i < j \le K$, yielding $\binom{K}{2}$ pairs per file.
2. For each pair $(S_i, S_j)$, strip all gap characters `'-'` to recover the raw unaligned protein sequences $(X, Y) = (\pi_{\Sigma}(S_i), \pi_{\Sigma}(S_j))$.
3. Record the ground-truth reference alignment $\mathcal{A}_{\text{ref}}(i, j)$ formed by the projected columns of $(S_i, S_j)$ in the MSA.

### 2.2 Preserving Core Residue Correspondences
Let $(u, v)$ represent an aligned pair where residue $u$ of sequence $X$ is aligned to residue $v$ of sequence $Y$:
- A correspondence $(u, v)$ is defined as **core-conserved** if and only if both residues appear in uppercase in the reference alignment:
  $$\mathcal{C}_{\text{ref}}(X, Y) = \left\{ (u, v) \mid S_i[k] = X[u] \text{ is uppercase} \land S_j[k] = Y[v] \text{ is uppercase} \right\}$$
- When evaluating a predicted test alignment $\mathcal{A}_{\text{test}}(X, Y)$ with aligned residue pairs $\mathcal{C}_{\text{test}}(X, Y)$, we compute accuracy strictly on the core set:
  $$Q = \frac{|\mathcal{C}_{\text{test}}(X, Y) \cap \mathcal{C}_{\text{ref}}(X, Y)|}{|\mathcal{C}_{\text{ref}}(X, Y)|}$$

---

## 3. Train, Validation, and Test Splitting Policy

To ensure strictly reproducible and unbiased evaluation:
1. **Family-Level Partitioning**: Splitting must occur at the **file / superfamily level**, not by randomly splitting pairwise combinations derived from the same file. Randomly splitting pairs from the same MSA would cause severe data leakage across folds.
2. **Split Proportions**:
   - **Training Set (60%)**: 181 alignment files (used for tuning gap penalty parameters and pseudocount priors).
   - **Validation Set (20%)**: 61 alignment files (used for hyperparameter selection and model comparison).
   - **Test Set (20%)**: 61 alignment files (held out for final evaluation and comparative reporting).
3. **Random Seed**: Fixed random seed (`seed = 42`) with stratified binning based on sequence divergence (% identity).

---

## 4. Preventing Data Leakage in PSSM Construction

When evaluating profile-based alignment (PSSM) against a target query sequence $Y$:
1. **Leave-One-Out Exclusion (LOO)**: The query sequence $Y$ **must never** be included in the MSA used to construct the PSSM. Including the query would artificially inflate profile frequencies for rare residues present only in the query.
   - Programmatically enforced by `proteinalign.split_msa_train_test(msa, test_id)`.
2. **Automated Leakage Detection**: The profile alignment engine tracks all training sequence headers in `profile.sequence_ids`. If an evaluation query ID matches any training identifier, `validate_no_data_leakage` immediately aborts execution by raising a `DataLeakageError`.
3. **Independent Background**: Proteome background frequencies ($p_a$) are derived from global databases (UniProt/Swiss-Prot via Robinson & Robinson, 1991), completely independent of the test family.

---

---

## 5. Formal Evaluation Metrics

### 5.1 Developer Score ($Q$-Score / Recall)
The $Q$-score measures the fraction of reference structural core residue pairs correctly predicted in the test alignment:
$$Q = \frac{|\mathcal{C}_{\text{test}}(X, Y) \cap \mathcal{C}_{\text{ref}}(X, Y)|}{|\mathcal{C}_{\text{ref}}(X, Y)|}$$
where $\mathcal{C}_{\text{ref}}$ is the set of aligned residue index pairs $(u, v)$ where both $X[u]$ and $Y[v]$ are uppercase (conserved core identified by MUSTANG from 3D PDB coordinates).

### 5.2 Modeler Score ($M$ / Precision)
The Modeler score measures the fraction of test aligned residue pairs that are present in the reference core:
$$M = \frac{|\mathcal{C}_{\text{test}}(X, Y) \cap \mathcal{C}_{\text{ref}}(X, Y)|}{|\mathcal{C}_{\text{test}}(X, Y)|}$$

### 5.3 Balanced $F_1$-Score
The harmonic mean of sensitivity ($Q$) and precision ($M$):
$$F_1 = 2 \cdot \frac{Q \cdot M}{Q + M}$$

### 5.4 Gap Topology Metrics
1. **Total Gap Characters:** Count of `'-'` characters across aligned sequences.
2. **Gap Openings Count:** Count of contiguous indel blocks ($\text{regex } `-+`$).
3. **Mean Indel Run Length:** $\bar{L}_{\text{gap}} = \frac{\text{Total Gap Characters}}{\text{Gap Openings}}$.

### 5.5 Score Recomputation and Invariant Auditing
Every alignment score is independently verified by direct recomputation from the aligned string characters:
$$\text{Score}_{\text{linear}} = \sum_{c \in \text{matches}} S(a_c, b_c) - d \cdot (\text{gaps}_1 + \text{gaps}_2)$$
$$\text{Score}_{\text{affine}} = \sum_{c \in \text{matches}} S(a_c, b_c) - \sum_{k \in \text{runs}_1} (d + (k-1)e) - \sum_{k \in \text{runs}_2} (d + (k-1)e)$$
Any discrepancy exceeding $10^{-4}$ triggers a validity failure.

---

## 6. Experimental Design and Controlled Comparisons

### 6.1 RQ1 — Gap Penalties: Linear vs Gotoh Affine
- **Hypothesis:** Linear gap penalties penalize every gap character equally, resulting in severe indel fragmentation (many single-residue gaps). Affine penalties group insertions/deletions into biologically realistic contiguous blocks, improving alignment accuracy.
- **Controlled Variables:** Identical sequence pairs, identical substitution matrix (BLOSUM62), global alignment mode.
- **Empirical Findings ($N=288$ paired comparisons across 61 test families):**
  - **Accuracy ($Q$-score):** Linear Median = 0.4082 vs Affine Median = 0.6150 (**Wilcoxon $p = 7.19 \times 10^{-16}$**, statistically significant superiority).
  - **Gap Openings:** Linear Mean = 20.06 vs Affine Mean = 8.73 (**Wilcoxon $p = 3.96 \times 10^{-45}$**).
  - **Mean Indel Run Length:** Linear Mean = 1.56 residues vs Affine Mean = 5.17 residues (**Wilcoxon $p = 1.14 \times 10^{-47}$**).
  - **Computational Cost:** Linear Median = 15.07 ms vs Affine Median = 24.60 ms (approx 1.63x constant-factor overhead for 3-layer DP).

### 6.2 RQ2 — Substitution Matrices: BLOSUM62 vs PAM250
- **Hypothesis:** BLOSUM62 (derived from local blocks with $\ge 62\%$ identity) provides sharper discrimination for moderate-to-high sequence homology, whereas PAM250 (derived from Dayhoff's Markov evolutionary model) was parameterized for deep divergences.
- **Empirical Findings ($N=288$ pairs):**
  - **Twilight Zone ($<15\%$ identity, $N=36$):** BLOSUM62 Median $Q = 0.0000$ vs PAM250 Median $Q = 0.0000$ ($p = 0.953$, no significant difference; both struggle without profile information).
  - **Low-Mid Homology ($15\text{--}25\%$ identity, $N=111$):** BLOSUM62 Median $Q = 0.5000$ vs PAM250 Median $Q = 0.3750$ (**Wilcoxon $p = 8.28 \times 10^{-4}$**, BLOSUM62 significantly superior).
  - **High Homology ($>25\%$ identity, $N=141$):** BLOSUM62 Median $Q = 0.8750$ vs PAM250 Median $Q = 0.8654$ (**Wilcoxon $p = 4.92 \times 10^{-3}$**).

### 6.3 RQ3 — Alignment Mode: Global Needleman-Wunsch vs Local Smith-Waterman
- **Objective:** Evaluate behavior under sequence length asymmetry and domain boundaries.
- **Empirical Findings:**
  - On high-length-disparity pairs ($|L_1 - L_2| / \max(L_1, L_2) > 0.25$, $N=52$): Smith-Waterman local alignment isolates the conserved catalytic/structural domain (median query coverage = 22.2%), whereas global alignment forces end-to-end matching, incurring severe terminal gap penalties.
  - On length-matched homologous pairs ($\le 0.25$, $N=236$): Global alignment achieves superior overall core recovery (Median $Q = 0.7054$).

### 6.4 RQ4 — Position-Specific Scoring Matrices (PSSM) under Leave-One-Out
- **Protocol:** Strictly leave-one-out (LOO) training: for each held-out target sequence $Y$, the profile is trained strictly on $\mathcal{M} \setminus \{Y\}$ with automated data leakage verification.
- **Empirical Findings ($N=183$ held-out targets):**
  - Profile core residue placement accuracy achieves Median $Q = 0.6000$ (Mean = 0.5349).
  - Fixed-matrix consensus alignment achieves Median $Q = 0.6932$ (Mean = 0.5990).
  - Position-specific profile scoring retains higher conservation specificity on variable loop regions, though fixed matrices benefit when the consensus sequence closely matches the target.

### 6.5 RQ5 — Parameter Sensitivity Surfaces
- **Grid:** Gap-opening $d \in [6.0, 16.0]$, gap-extension $e \in [0.5, 3.0]$ ($N=450$ parameter evaluations).
- **Pseudocount Sweep:** Regularization weight $\alpha \in [0.1, 10.0]$ ($N=60$ evaluations).
- **Optimal Region:** Stable plateau identified around $(d^*, e^*) = (10.0, 1.0)$ to $(12.0, 1.0)$.
- **Regularization:** PSSM accuracy remains stable for $\alpha \in [0.5, 2.0]$, with moderate degradation when unregularized ($\alpha \to 0$) or over-smoothed toward background ($\alpha \ge 5.0$).

### 6.6 RQ6 — Computational Scaling and Asymptotic Growth
- **Benchmark:** Wall-clock execution time (`time.perf_counter_ns`) and peak heap memory (`tracemalloc`) across lengths $L \in [50, 100, 200, 300, 500, 750, 1000]$ (length products $m \cdot n \in [2500, 1000000]$).
- **Empirical Power-Law Scaling:** $\log_{10}(\text{Runtime}) = \alpha + \beta \log_{10}(m \cdot n)$
  - Gotoh Global Affine: $\beta = 1.096$ ($R^2 = 0.9914$)
  - Gotoh Local Affine: $\beta = 1.094$ ($R^2 = 0.9912$)
  - Needleman-Wunsch Linear: $\beta = 1.083$ ($R^2 = 0.9907$)
  - PSSM Global Affine: $\beta = 1.105$ ($R^2 = 0.9922$)
  - Smith-Waterman Linear: $\beta = 1.084$ ($R^2 = 0.9913$)
- **Memory Scaling:** Gotoh 3-state affine DP consumes exactly $3.00\times$ the heap memory of linear single-matrix DP across all sequence lengths (e.g. at $L=1000$: NW = 31.4 MB, Gotoh Affine = 94.1 MB).

