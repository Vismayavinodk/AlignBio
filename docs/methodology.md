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
2. **Redundancy Filtering**: If any homolog in the family MSA shares $\ge 90\%$ sequence identity with the query, that homolog is also pruned from profile construction to avoid memorization.
3. **Independent Background**: Proteome background frequencies ($p_a$) must be derived from global databases (e.g. UniProt/Swiss-Prot via Robinson & Robinson), not from the test family.

---

## 5. Experimental Plan and Comparative Evaluation

The Phase 3 benchmark experiments will evaluate five orthogonal dimensions:

### 5.1 Alignment Accuracy ($Q$-Score and Modeler Score)
Compare the alignment accuracy of:
1. Needleman–Wunsch (Linear Gap) vs. Gotoh (Affine Gap).
2. BLOSUM62 vs. PAM250 vs. Unit Identity Matrix.
3. Pairwise sequence-sequence vs. Profile-sequence (PSSM).

### 5.2 Computational Complexity and Runtime Benchmarking
- Measure wall-clock runtime using high-precision timers (`time.perf_counter_ns`).
- Benchmark runtime scaling as a function of sequence length product $m \cdot n$ across lengths from 50 to 1,000 residues.
- Validate the theoretical $\mathcal{O}(m \cdot n)$ quadratic scaling law.

### 5.3 Memory Profiling
- Measure peak heap memory allocation using Python's `tracemalloc`.
- Compare full-matrix traceback storage ($\mathcal{O}(m \cdot n)$) versus score-only DP ($\mathcal{O}(\min(m, n))$).

### 5.4 Sensitivity to Gap Penalty Parameters
- Grid search over gap opening penalties $d \in [6.0, 16.0]$ in steps of 2.0.
- Grid search over gap extension penalties $e \in [0.5, 3.0]$ in steps of 0.5.
- Generate contour / heat map plots illustrating $Q$-score sensitivity surfaces across protein identity tiers ($<20\%$, $20-40\%$, $>40\%$).
