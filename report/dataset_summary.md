# Dataset Summary: SABmark / SABRE Benchmark (`sabrem`)

## Overview
The dataset analyzed for the **ProteinAlign Lab** benchmark evaluation is located at `data/raw/bench1.0/sabrem/ref`. This benchmark corresponds to the **SABRE** multiple sequence alignments derived from **SABmark v1.65**, re-aligned using the **MUSTANG** structural alignment tool (`MU-SABRE`).

## Execution & Verification
The dataset was verified using [`summarize_dataset.py`](file:///d:/tcs/packages/sem9/MM/ProteinAlign-Lab/summarize_dataset.py) on October 14, 2026.

### Summary Metrics
| Metric | Value |
| :--- | :--- |
| **Benchmark Directory** | `data/raw/bench1.0/sabrem/ref` |
| **Total Files Processed** | 303 |
| **Valid Alignment Files** | 303 (100%) |
| **Invalid Files** | 0 (0%) |
| **Total Aligned Sequences** | 1,650 |
| **Sequences per Alignment (Min / Max / Mean)** | 3 / 24 / 5.45 |
| **Alignment Length (Min / Max / Mean)** | 42 / 1,076 / 264.69 columns |

## File Format & Case Convention
Each file in `sabrem/ref/` is formatted in FASTA format:
1. **Equal Row Length**: In every file, each sequence row has an identical number of characters (confirming a valid aligned matrix).
2. **Uppercase vs. Lowercase Residues**:
   - **Uppercase characters**: Denote biologically/structurally conserved core positions reliably identified by MUSTANG. These columns form the ground-truth core for benchmark assessment (e.g., $Q$-score / developer score).
   - **Lowercase characters**: Denote unconserved flanking or loop regions where alignment is uncertain; by benchmark convention, lowercase residues are ignored during core accuracy evaluation.
3. **Gap Characters (`-`)**: Represent deletions and insertions necessary to align the homologous structural elements across the protein family.
