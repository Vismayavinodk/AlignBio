# Dataset Inspection and Provenance Report

**Project**: ProteinAlign Lab  
**Date**: October 14, 2026  
**Auditor**: Automated Verification Pipeline (`summarize_dataset.py`)  
**Status**: Verified & Validated (303/303 files valid)

---

## 1. Executive Summary
A comprehensive audit and integrity check was performed on the benchmark archive `bench.tar.gz` and its extracted dataset tree located at `data/raw/bench1.0/sabrem/ref`. The dataset was confirmed to be **MU-SABRE**, a multiple sequence alignment benchmark derived from **SABmark v1.65** re-aligned using the structural alignment algorithm **MUSTANG**. All 303 reference files were validated for format compliance, sequence integrity, and matrix consistency. No corrupted or malformed files were found.

---

## 2. Dataset Provenance and Attribution

### 2.1 Benchmark Identity
- **Archive File**: `bench.tar.gz` (21,414,017 bytes)
- **Extracted Root**: `data/raw/bench1.0/`
- **Target Evaluation Subdirectory**: `data/raw/bench1.0/sabrem/ref/`
- **Compiler**: Dr. Robert C. Edgar (author of MUSCLE)
- **Benchmark Name**: `MU-SABRE` (`sabrem`)
- **Underlying Benchmark**: SABmark v1.65 (Sequence Alignment Benchmark)
- **Structural Re-aligner**: MUSTANG (Konagurthu et al., 2006)

### 2.2 Academic Citations
1. **SABmark**: Van Walle, I., Lasters, I., & Wyns, L. (2005). SABmark—a benchmark for sequence alignment that covers the entire known fold space. *Bioinformatics*, 21(7), 1267–1268.
2. **MUSTANG**: Konagurthu, A. S., et al. (2006). MUSTANG: A multiple structural alignment algorithm. *Proteins*, 64(3), 559–574.
3. **Benchmark Set**: Edgar, R. C. (2010). Quality measures for protein alignment benchmarks. *Nucleic Acids Research*, 38(7), 2145–2153.

---

## 3. Dataset-Wide Verification Statistics

The dataset was parsed using `summarize_dataset.py`. The quantitative findings are summarized below:

| Metric | Measured Value | Validation Status |
| :--- | :--- | :--- |
| **Total Reference Files** | 303 | Complete |
| **Valid Alignment Files** | 303 (100.0%) | Passed |
| **Malformed / Invalid Files** | 0 (0.0%) | Passed |
| **Total Aligned Sequences** | 1,650 | Complete |
| **Sequences per Alignment (Min)** | 3 | Valid |
| **Sequences per Alignment (Max)** | 24 | Valid |
| **Sequences per Alignment (Mean)** | 5.45 | Valid |
| **Alignment Length (Min)** | 42 columns | Valid |
| **Alignment Length (Max)** | 1,076 columns | Valid |
| **Alignment Length (Mean)** | 264.69 columns | Valid |

---

## 4. File Structure and Annotation Conventions

1. **FASTA Alignment Matrix**: Each file contains $K$ sequences with equal character lengths including gap characters (`'-'`).
2. **Residue Case Distinction**:
   - **Uppercase (`A-Z`)**: Denotes structurally conserved core positions identified by MUSTANG based on superposition of known 3D PDB crystal structures. These residues serve as the gold-standard ground truth for $Q$-score evaluation.
   - **Lowercase (`a-z`)**: Denotes variable loops and terminal extensions where structural alignment is ambiguous. By convention, lowercase positions must be excluded from core accuracy calculations.
3. **Gap Symbol**: Hyphens (`'-'`) represent structural deletions and insertions.

---

## 5. Known Limitations and Preprocessing Guidelines

1. **Pairwise Decomposition**: Because the dataset consists of multiple sequence alignments, pairwise experiments require extracting all $\binom{K}{2}$ pairs from each file and degapping them to form unaligned input pairs.
2. **Divergence Range (Twilight Zone)**: SABmark contains twilight-zone protein families with $<20\%$ sequence identity, presenting a rigorous test for substitution matrices and gap penalties.
3. **Data Leakage Safeguard**: PSSMs constructed from an MSA must strictly exclude the query sequence and its nearest homologs ($>90\%$ identity) during evaluation.
