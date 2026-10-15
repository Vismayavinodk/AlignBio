# ProteinAlign Lab — Professor-Facing Project Defense & Viva Guide

**Project Title:** ProteinAlign Lab: Mathematical Modelling and Comparative Evaluation of Protein Sequence Alignment Using Dynamic Programming and Position-Specific Scoring Matrices  
**Candidate Preparation Guide:** Oral Examination, Viva Defense, and Academic Presentation  
**Date:** October 2026  

---

## 1. What is the exact mathematical problem?

Protein sequence alignment is an optimization problem over discrete combinatorial path spaces. Given two sequences $X = (x_1, \dots, x_m)$ and $Y = (y_1, \dots, y_n)$ over the standard 20-letter amino acid alphabet $\Sigma$, an alignment $\mathcal{A}$ is a pair of equal-length strings $(X', Y')$ of length $L \ge \max(m, n)$ over $\Sigma \cup \{'-'\}$ such that removing all gap characters `'-'` restores $X$ and $Y$ respectively, with the constraint that no column contains two gaps (`'-'`, `'-'`).

The mathematical problem is to find an alignment $\mathcal{A}^* \in \Omega(X, Y)$ that maximizes the objective scoring function:
$$\mathcal{A}^* = \arg\max_{\mathcal{A} \in \Omega(X, Y)} \mathcal{S}(\mathcal{A})$$
where the total alignment score $\mathcal{S}(\mathcal{A})$ decomposes into residue-to-residue substitution scores and gap insertion/deletion costs:
$$\mathcal{S}(\mathcal{A}) = \sum_{c \in \text{aligned columns}} S(X'[c], Y'[c]) - \sum_{k \in \text{gap runs}} \gamma(k)$$
The size of the unconstrained alignment space $|\Omega(X, Y)|$ grows as $\binom{m + n}{m} \approx \frac{2^{2n}}{\sqrt{\pi n}}$, making brute-force enumeration computationally impossible even for modest lengths ($m=n=100 \implies \approx 10^{59}$ alignments).

---

## 2. Why is dynamic programming suitable?

Dynamic programming (DP) is suitable because protein sequence alignment exhibits:
1. **Optimal Substructure (Bellman's Principle of Optimality):** The optimal alignment of prefix $X[1..i]$ and $Y[1..j]$ is composed of the optimal alignment of an immediate sub-problem ($X[1..i-1], Y[1..j-1]$, $X[1..i-1], Y[1..j]$, or $X[1..i], Y[1..j-1]$) plus the transition cost of the final column $(x_i, y_j)$, $(x_i, -)$, or $(-, y_j)$.
2. **Overlapping Subproblems:** Evaluating alignments involves repeatedly evaluating identical prefix alignments $(X[1..i], Y[1..j])$. Storing sub-problem solutions in a 2D table of size $(m+1) \times (n+1)$ reduces exponential combinatorial complexity $\mathcal{O}(2^{m+n})$ to quadratic polynomial time $\mathcal{O}(m \cdot n)$.

---

## 3. What does each recurrence optimize?

- **Needleman–Wunsch Global Alignment:** Optimizes the end-to-end alignment score across the entire length of both sequences:
  $$F(i, j) = \max \begin{cases} F(i-1, j-1) + S(x_i, y_j) & \text{(match / substitution)} \\ F(i-1, j) - d & \text{(gap in } Y \text{ / deletion)} \\ F(i, j-1) - d & \text{(gap in } X \text{ / insertion)} \end{cases}$$
  Boundary conditions $F(i, 0) = -i \cdot d$ and $F(0, j) = -j \cdot d$ penalize terminal overhanging residues, forcing complete global correspondence.
- **Smith–Waterman Local Alignment:** Optimizes the maximum score of any pair of contiguous subsequences:
  $$H(i, j) = \max \begin{cases} 0 & \text{(local restart cutoff)} \\ H(i-1, j-1) + S(x_i, y_j) \\ H(i-1, j) - d \\ H(i, j-1) - d \end{cases}$$
  The zero floor discards negative-scoring prefix alignments, allowing the DP to isolate high-scoring local domains without penalizing flanking non-homologous regions. The global optimum is $\max_{i, j} H(i, j)$.
- **Gotoh (1982) Affine Gap Dynamic Programming:** Optimizes alignment score under a two-parameter affine gap cost model $G(k) = d + (k-1)e$ using three coupled recurrence layers:
  - $M(i, j)$: Optimal alignment of prefixes $X[1..i]$ and $Y[1..j]$ ending with a matched/substituted residue pair.
  - $I_x(i, j)$: Optimal alignment ending with a gap in $Y$ (deletion of $x_i$).
  - $I_y(i, j)$: Optimal alignment ending with a gap in $X$ (insertion of $y_j$).

---

## 4. Why are gap-opening and gap-extension penalties different?

From molecular biology, an insertion or deletion (indel) event arises from a single mutational event—such as DNA replication slippage or unequal crossing-over. Once the protein loop opens to accommodate an indel, accommodating additional consecutive residues requires relatively little additional free energy disruption compared to the initial structural disruption of opening the backbone.

- **Linear Gap Model:** $\gamma_{\text{linear}}(k) = k \cdot d$. Every gap character costs $d$. A gap of length 5 costs $5d$, while 5 isolated single-residue gaps also cost $5d$. To maximize substitution scores, DP fragments indels into single-residue gaps across secondary structure elements.
- **Affine Gap Model:** $\gamma_{\text{affine}}(k) = d + (k-1)e$ where $d \gg e$ (e.g. $d=10, e=1$). Opening a gap incurs a severe penalty $d=10$, but extending it costs only $e=1$. Under this model, 5 contiguous gaps cost $10 + 4(1) = 14$, whereas 5 isolated gaps cost $5 \times 10 = 50$. This strongly favors grouping gaps into biologically realistic contiguous loop blocks.
- **Empirical Validation (RQ1):** On MU-SABRE ($N=288$), Gotoh affine gaps reduced mean gap openings from 20.06 to 8.73 ($p = 3.96 \times 10^{-45}$) and increased mean gap length from 1.56 to 5.17 residues ($p = 1.14 \times 10^{-47}$), significantly increasing developer sensitivity ($Q$-score 0.6150 vs 0.4082, $p = 7.19 \times 10^{-16}$).

---

## 5. How do BLOSUM62 and PAM250 differ conceptually?

| Property | PAM250 (Dayhoff et al., 1978) | BLOSUM62 (Henikoff & Henikoff, 1992) |
| :--- | :--- | :--- |
| **Philosophical Basis** | Evolutionary Markov model of point mutations. | Direct empirical observation of conserved protein blocks. |
| **Source Data** | 71 closely related protein families ($\ge 85\%$ identity). | BLOCKS database of ungapped local alignments. |
| **Divergence Model** | Derived PAM1 from closely related proteins and extrapolated by matrix multiplication: $M_{250} = (M_1)^{250}$. | Clustered all sequences sharing $\ge 62\%$ identity to prevent over-representation, counting observed substitutions directly. |
| **Target Divergence** | Designed for deep evolutionary distance (250 PAM $\approx 80\%$ sequence divergence). | Parameterized for intermediate evolutionary distances ($20\text{--}60\%$ identity). |
| **Empirical Finding (RQ2)** | Under identical Gotoh penalties ($d=10, e=1$), PAM250 achieved median $Q=0.3750$ in the $15\text{--}25\%$ identity tier, compared to $Q=0.5000$ for BLOSUM62 ($p = 8.28 \times 10^{-4}$). | Significantly outperformed PAM250 across mid- and high-homology tiers. In the twilight zone ($<15\%$), neither fixed matrix succeeded without profiles. |

---

## 6. What does a PSSM capture that a fixed substitution matrix cannot?

A fixed substitution matrix treats every position in a protein identically: replacing a Cysteine with a Serine always receives score $-1$, regardless of whether that Cysteine is a catalytic active-site nucleophile or a non-essential surface loop residue.

A **Position-Specific Scoring Matrix (PSSM)** builds an evolutionary profile where each column corresponds to a specific structural position in a protein family:
- **Functionally Critical Positions:** Positions that are strictly conserved across a family (e.g. active-site Histidine) exhibit an observed probability $P_{i, H} \approx 1.0$. The log-odds score $\lambda \log_2(P_{i, H} / q_H)$ becomes a large positive value ($+8\text{--}+10$), while all substitutions are severely penalized.
- **Tolerant Positions:** Surface loop positions that naturally vary between Leucine, Isoleucine, and Valine receive moderate scores for all hydrophobic residues.
This position-dependent evolutionary conservation signature allows PSSMs to detect structural homology where pairwise sequence similarity has eroded.

---

## 7. How are profile probabilities estimated and smoothed?

Given a multiple sequence alignment of depth $K$, observed counts $c_{i, a}$ at column $i$ for amino acid $a$ are computed:
$$c_{i, a} = \sum_{k=1}^K \mathbb{I}(S_k[i] = a) + \text{fractional weights for ambiguous codes } (B, Z, X)$$

To prevent zero-frequency artifacts (where an unobserved residue would produce $P_{i, a} = 0 \implies \log(0) = -\infty$), we apply **Dirichlet-style pseudocount regularization** using the Robinson & Robinson (1991) proteome background distribution $q_a$:
$$P_{i, a} = \frac{c_{i, a} + \alpha q_a}{N_i + \alpha}$$
where $N_i = \sum_{a' \in \Sigma} c_{i, a'}$ is the total non-gap count at column $i$, and $\alpha \ge 0$ is the pseudocount smoothing parameter.

Log-odds scores are computed in half-bits ($\lambda = 2.0, b = 2.0$):
$$S_{i, a} = 2.0 \log_2 \left( \frac{P_{i, a}}{q_a} \right)$$
When a residue occurs exactly at background probability ($P_{i, a} = q_a$), $S_{i, a} = 2.0 \log_2(1.0) = 0.0$.

---

## 8. What is the computational complexity of each implemented method?

| Algorithm | State Space (Matrices) | Time Complexity | Auxiliary Space Complexity | Empirical Memory ($L=1000$) |
| :--- | :--- | :--- | :--- | :--- |
| **Needleman–Wunsch Linear** | 1 matrix ($F$) | $\mathcal{O}(m \cdot n)$ | $\mathcal{O}(m \cdot n)$ | 31.4 MB |
| **Smith–Waterman Linear** | 1 matrix ($H$) | $\mathcal{O}(m \cdot n)$ | $\mathcal{O}(m \cdot n)$ | 31.4 MB |
| **Gotoh Global Affine** | 3 matrices ($M, I_x, I_y$) | $\mathcal{O}(m \cdot n)$ | $\mathcal{O}(m \cdot n)$ | 94.1 MB ($3.00\times$) |
| **Gotoh Local Affine** | 3 matrices ($M, I_x, I_y$) | $\mathcal{O}(m \cdot n)$ | $\mathcal{O}(m \cdot n)$ | 94.1 MB ($3.00\times$) |
| **PSSM Global Affine** | 3 matrices ($M, I_x, I_y$) | $\mathcal{O}(L_{\text{prof}} \cdot n)$ | $\mathcal{O}(L_{\text{prof}} \cdot n)$ | 94.1 MB |

Empirical log-log regressions confirmed runtime power-law exponents $\beta \in [1.083, 1.105]$ ($R^2 > 0.99$), demonstrating agreement with theoretical quadratic complexity.

---

## 9. How was correctness tested independently?

Correctness was established across multiple independent validation tiers:
1. **Mathematical Invariants:** Self-alignment identity ($PID = 100\%$, score = sum of self-substitutions), symmetry ($align(A, B).score == align(B, A).score$), and degapping reconstruction ($X_{\text{aln}}.replace('-', '') == X$).
2. **Independent Score Recomputation:** An independent function [`recompute_score`](file:///d:/tcs/packages/sem9/MM/ProteinAlign-Lab/src/proteinalign/evaluation.py#L150-L195) recalculates alignment scores directly from output strings using regex indel parsing, verifying 100% of benchmark records (`assert abs(score - recomputed_score) < 1e-4`).
3. **External Reference Cross-Validation:** Compared outputs against Biopython (`Bio.Align.PairwiseAligner`) in [`tests/test_biopython_compat.py`](file:///d:/tcs/packages/sem9/MM/ProteinAlign-Lab/tests/test_biopython_compat.py).
4. **Hand-Calculated Benchmarks:** Exact hand-derived examples in [`tests/test_experiments.py`](file:///d:/tcs/packages/sem9/MM/ProteinAlign-Lab/tests/test_experiments.py) for $Q$-score, Modeler score, $F_1$, and Gotoh indels ($k=1, 2, 20$).

---

## 10. What is the dataset's provenance and what are its limitations?

- **Provenance:** The dataset is **MU-SABRE**, compiled by Dr. Robert C. Edgar from **SABmark v1.65** and re-aligned using the **MUSTANG** structural superposition algorithm (Konagurthu et al., 2006). It consists of 303 FASTA files in `data/raw/bench1.0/sabrem/ref`.
- **Limitations:**
  1. *Structural Core Only:* Ground-truth alignment is defined only for uppercase core residues. Lowercase loop residues cannot be evaluated for correctness.
  2. *Finite Sample Depth:* Family depths range from 3 to 24 sequences (median 4). Profile construction in shallow families ($K=3$) is constrained by limited evolutionary sampling.
  3. *Reference Re-alignment Bias:* Structural alignments depend on MUSTANG's coordinate superposition heuristics and parameterization.

---

## 11. How was evaluation leakage prevented?

1. **Superfamily-Level Split:** Data partitioning was performed strictly at the **file / superfamily level** (181 Train, 61 Val, 61 Test), never by splitting sequence pairs from the same family.
2. **Leave-One-Out (LOO) PSSM Construction:** When evaluating profile alignment against target sequence $Y$, $Y$ is strictly removed from the MSA: $\mathcal{M}_{\text{train}} = \mathcal{M} \setminus \{Y\}$.
3. **Automated Leakage Assertion:** `validate_no_data_leakage` checks $Y$'s identifier against `profile.sequence_ids`, raising `DataLeakageError` upon any overlap.
4. **Independent Background:** Background amino acid frequencies are drawn from global Swiss-Prot proteome frequencies (Robinson & Robinson, 1991), completely external to the benchmark set.

---

## 12. What are the strongest experimentally supported findings?

1. **Affine Penalty Superiority:** Gotoh affine gaps significantly outperform linear penalties in developer sensitivity ($Q = 0.6150$ vs $0.4082$, Wilcoxon $p = 7.19 \times 10^{-16}$) while reducing indel fragmentation by $56.5\%$ ($p = 3.96 \times 10^{-45}$).
2. **BLOSUM62 vs PAM250 Differentiation:** BLOSUM62 is statistically superior in the $15\text{--}25\%$ identity tier ($p = 8.28 \times 10^{-4}$) and $>25\%$ tier ($p = 4.92 \times 10^{-3}$).
3. **Local Alignment Domain Isolation:** Under severe length asymmetry ($>25\%$), Smith–Waterman isolates local domains (median coverage 22.2%) without incurring terminal overhanging gap penalties.
4. **Memory Scaling Constant:** Peak heap memory for Gotoh 3-state DP is exactly $3.00\times$ linear DP across all lengths ($L=50$ to $L=1000$).

---

## 13. What does the project not establish?

1. **Biological Truth of Variable Loops:** The project evaluates developer sensitivity only on MUSTANG-defined core residues; it does not establish the accuracy of loop alignments where 3D coordinate superposition is ambiguous.
2. **Universal Superiority of PSSMs in Deep Twilight:** In shallow families ($K=3$), profile construction on $K-1=2$ sequences provides limited evolutionary depth, meaning PSSMs do not outperform fixed consensus alignment across all twilight sets without deeper sequence databases (e.g. UniRef).
3. **Hardware-Independent Absolute Run Times:** Millisecond execution times depend on the host CPU clock frequency and cache structure.

---

## 14. What is the main limitation?

The primary algorithmic limitation is **quadratic $\mathcal{O}(m \cdot n)$ space complexity**. Storing the full $(m+1) \times (n+1)$ dynamic programming matrix requires 94.1 MB of heap memory at $L=1000$ residues. For large multi-domain proteins ($L > 5000$) or genomic sequences, full-matrix traceback causes memory exhaustion unless space-efficient algorithms (e.g. Hirschberg's divide-and-conquer) are deployed.

---

## 15. If given two additional weeks, what would be the most valuable next step?

The single most valuable technical extension would be implementing **Hirschberg's (1975) linear-space dynamic programming algorithm** for Gotoh affine gaps. Hirschberg's divide-and-conquer algorithm computes optimal global alignments in $\mathcal{O}(m \cdot n)$ time while reducing auxiliary space complexity from $\mathcal{O}(m \cdot n)$ to $\mathcal{O}(\min(m, n))$, slashing peak memory at $L=1000$ from 94 MB to under 100 KB.
