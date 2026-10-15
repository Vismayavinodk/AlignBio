# ProteinAlign Lab — Mathematical Foundations and Algorithmic Formulations

This document provides the formal mathematical derivations, state equations, boundary conditions, probabilistic foundations, and algorithmic specifications for pairwise and profile protein sequence alignment.

---

## 1. Mathematical Definition of Protein Sequence Alignment

### 1.1 Alphabet and Sequence Space
Let $\Sigma$ denote the standard IUPAC amino acid alphabet of size $|\Sigma| = 20$:
$$\Sigma = \{\text{A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V}\}$$
Let $\Sigma' = \Sigma \cup \{'-'\}$ be the alphabet augmented with the gap character $'-'$, where $|'-'| = 1$.

Given two primary protein sequences:
- Query sequence: $X = (x_1, x_2, \dots, x_m) \in \Sigma^m$ of length $|X| = m$.
- Target sequence: $Y = (y_1, y_2, \dots, y_n) \in \Sigma^n$ of length $|Y| = n$.

### 1.2 Alignment Representation
An alignment $\mathcal{A}(X, Y)$ is defined as a $2 \times L$ matrix $\begin{pmatrix} X' \\ Y' \end{pmatrix}$ over $\Sigma'$:
$$X' = (x'_1, x'_2, \dots, x'_L), \quad Y' = (y'_1, y'_2, \dots, y'_L)$$
satisfying three fundamental structural properties:
1. **Length Constraint**: $\max(m, n) \le L \le m + n$.
2. **Reconstruction (Degap) Invariant**: Removing all gap symbols $'-'$ from $X'$ yields $X$, and removing all $'-'$ from $Y'$ yields $Y$:
   $$\pi_{\Sigma}(X') = X \quad \text{and} \quad \pi_{\Sigma}(Y') = Y$$
   where $\pi_{\Sigma}$ is the projection homomorphism erasing $'-'$.
3. **Non-Degeneracy (No Double Gaps)**: No column consists entirely of gaps:
   $$\forall k \in \{1, \dots, L\}, \quad \neg (x'_k = \text{'-'} \land y'_k = \text{'-'})$$

### 1.3 Optimization Objective
Let $\sigma: \Sigma' \times \Sigma' \to \mathbb{R}$ denote a scoring function over pairs of aligned characters. The total score $\mathcal{S}(\mathcal{A})$ of an alignment $\mathcal{A}$ under an additive scoring scheme is:
$$\mathcal{S}(\mathcal{A}) = \sum_{k=1}^L \sigma(x'_k, y'_k)$$
The sequence alignment problem is to determine an optimal alignment $\mathcal{A}^*$ that maximizes the objective score:
$$\mathcal{A}^* = \arg\max_{\mathcal{A} \in \Omega(X, Y)} \mathcal{S}(\mathcal{A})$$
where $\Omega(X, Y)$ is the complete combinatorial space of all valid alignments between $X$ and $Y$.

---

## 2. Global Alignment: The Needleman–Wunsch Algorithm

### 2.1 Problem Formulation
Global alignment seeks an end-to-end alignment covering the entirety of both sequences from $x_1$ to $x_m$ and $y_1$ to $y_n$.

### 2.2 Dynamic Programming Recurrence Relation
Let $F(i, j)$ denote the maximum alignment score between the prefix $X_{1..i} = (x_1, \dots, x_i)$ and prefix $Y_{1..j} = (y_1, \dots, y_j)$ for $0 \le i \le m$ and $0 \le j \le n$.

By Bellman's principle of optimality, an optimal alignment of prefixes $X_{1..i}$ and $Y_{1..j}$ must terminate in one of three mutually exclusive column configurations:
1. $x_i$ aligned with $y_j$ (substitution/match): preceding score is $F(i-1, j-1) + S(x_i, y_j)$.
2. $x_i$ aligned with a gap in $Y$ (deletion in $Y$ / insertion in $X$): preceding score is $F(i-1, j) - d$.
3. $y_j$ aligned with a gap in $X$ (insertion in $Y$ / deletion in $X$): preceding score is $F(i, j-1) - d$.

Thus, the recurrence relation is:
$$F(i, j) = \max \begin{cases}
  F(i-1, j-1) + S(x_i, y_j) & \text{[Diagonal: substitution of } x_i \text{ with } y_j \text{]} \\
  F(i-1, j) - d            & \text{[Vertical: gap in } Y \text{ / deletion]} \\
  F(i, j-1) - d            & \text{[Horizontal: gap in } X \text{ / insertion]}
\end{cases}$$
where $S(a, b)$ is the substitution score for residue pair $(a, b)$ and $d > 0$ is the linear gap penalty.

### 2.3 Boundary Conditions and Initialization
The base cases represent aligning an empty sequence against a sequence of length $k$, which incurs $k$ gap penalties:
$$F(0, 0) = 0$$
$$F(i, 0) = -i \cdot d \quad \forall i \in \{1, 2, \dots, m\}$$
$$F(0, j) = -j \cdot d \quad \forall j \in \{1, 2, \dots, n\}$$

### 2.4 Matrix Filling and Traceback
- **Filling Order**: Row-major or column-major order; cell $(i, j)$ depends strictly on $(i-1, j-1)$, $(i-1, j)$, and $(i, j-1)$.
- **Termination / Global Score**: The optimal global alignment score is stored in the bottom-right terminal cell:
  $$\text{Score}_{\text{global}} = F(m, n)$$
- **Traceback**: A traceback pointer matrix $T(i, j) \in \{\text{Diag}, \text{Up}, \text{Left}\}$ records the maximizing argument at each step. Traceback begins at $(m, n)$ and terminates when $(0, 0)$ is reached.

---

## 3. Local Alignment: The Smith–Waterman Algorithm

### 3.1 Problem Formulation
Local alignment identifies the pair of substrings $X[i_1..i_2]$ and $Y[j_1..j_2]$ ($1 \le i_1 \le i_2 \le m$, $1 \le j_1 \le j_2 \le n$) whose mutual alignment achieves the highest score, without penalizing non-homologous flanking terminal residues.

### 3.2 Dynamic Programming Recurrence Relation
Let $H(i, j)$ denote the maximum score of a local alignment terminating precisely at $x_i$ and $y_j$:
$$H(i, j) = \max \begin{cases}
  0                         & \text{[Reset: terminate path if score drops below zero]} \\
  H(i-1, j-1) + S(x_i, y_j) & \text{[Diagonal: match / mismatch]} \\
  H(i-1, j) - d            & \text{[Vertical: gap in } Y \text{]} \\
  H(i, j-1) - d            & \text{[Horizontal: gap in } X \text{]}
\end{cases}$$

### 3.3 Boundary Conditions
Because any local alignment can begin at any position without penalty:
$$H(i, 0) = 0 \quad \forall i \in \{0, 1, \dots, m\}$$
$$H(0, j) = 0 \quad \forall j \in \{0, 1, \dots, n\}$$

### 3.4 Termination and Traceback
- **Local Optimal Score**: Unlike Needleman–Wunsch, the optimal score is not located at $(m, n)$, but at the global matrix maximum:
  $$\text{Score}_{\text{local}} = \max_{\substack{0 \le i \le m \\ 0 \le j \le n}} H(i, j) = H(i^*, j^*)$$
- **Traceback**: Starts at the maximal coordinate $(i^*, j^*)$ and steps backward until a cell with value $H(i, j) = 0$ is encountered. The starting coordinate $(i_{\text{start}}, j_{\text{start}})$ and ending coordinate $(i^*, j^*)$ define the local aligned motif boundaries.

---

## 4. Comparison of Global vs. Local Alignment

| Criterion | Needleman–Wunsch (Global) | Smith–Waterman (Local) |
| :--- | :--- | :--- |
| **Alignment Scope** | Full length of both sequences ($1..m$ and $1..n$) | Optimal sub-segments ($i_1..i_2$ and $j_1..j_2$) |
| **End Gap Penalty** | Penalized at boundary ($F(i, 0) = -i \cdot d$) | Zero penalty ($H(i, 0) = 0$) |
| **Zero-Floor Cutoff** | No floor; scores can be negative | Lower bounded by $0$ |
| **Optimal Score Location** | Fixed at bottom-right cell: $F(m, n)$ | Matrix maximum: $\max_{i, j} H(i, j)$ |
| **Traceback Termination** | At origin $(0, 0)$ | Upon reaching any cell with score $0$ |
| **Biological Use Case** | Closely related proteins of similar length and conserved domain architecture across the entire sequence. | Multi-domain proteins sharing a conserved domain, remote homolog detection, or locating active sites within divergent flanking sequences. |

---

## 5. Substitution Scoring Matrices: BLOSUM62 and PAM250

### 5.1 Log-Odds Ratio Derivation
A substitution matrix provides a log-odds measure evaluating whether an observed aligned residue pair $(a, b)$ is more likely under a hypothesis of homology ($\mathcal{H}_1$) versus random chance ($\mathcal{H}_0$):
- **Target Frequency $q_{ab}$**: Probability of residue $a$ aligning with $b$ in homologous proteins.
- **Background Frequencies $p_a, p_b$**: Independent occurrence frequencies in the proteome.
- **Log-Odds Score**:
  $$S(a, b) = \frac{1}{\lambda} \log \left( \frac{q_{ab}}{p_a \cdot p_b} \right)$$
where $\lambda$ is a scaling factor (typically $\lambda = \frac{\ln 2}{2}$ for half-bit units).

### 5.2 Mathematical Requirement for Local Alignment
For Smith–Waterman local alignment to produce finite, biologically meaningful motifs rather than spanning the entire sequence length, the expected substitution score under the background distribution must be strictly negative:
$$\mathbb{E}[S] = \sum_{a \in \Sigma} \sum_{b \in \Sigma} p_a \cdot p_b \cdot S(a, b) < 0$$
Simultaneously, at least one pair must satisfy $S(a, b) > 0$.

### 5.3 BLOSUM vs. PAM Matrices
- **PAM (Point Accepted Mutation) Series (Dayhoff et al.)**:
  - Derived from an explicit Markov model of evolutionary divergence based on closely related sequences ($>85\%$ identity).
  - PAM$k$ denotes $k$ accepted point mutations per 100 residues.
  - **PAM250**: Represents deep evolutionary distance (~20% sequence identity, ~2.5 mutations per site). Suitable for distantly related sequences.
- **BLOSUM (Blocks Substitution Matrix) Series (Henikoff & Henikoff)**:
  - Derived from observed ungapped local blocks of aligned protein families directly without an evolutionary extrapolation model.
  - Sequences with identity $\ge X\%$ are clustered to eliminate bias.
  - **BLOSUM62**: Clusters sequences with $\ge 62\%$ identity. It is the gold standard default for general-purpose pairwise alignment, balancing sensitivity for medium-to-distantly related proteins.

---

## 6. Gap Penalty Models

### 6.1 Linear Gap Penalty
A linear penalty assigns a uniform cost proportional to the number of gap characters:
$$\gamma_{\text{linear}}(k) = -k \cdot d \quad (d > 0)$$
**Limitation**: A single evolutionary insertion or deletion (indel) event often involves multiple contiguous amino acids due to loop insertions or DNA replication slippage. Penalizing a 5-residue gap as 5 independent indels ($-5d$) is excessively punitive and unnaturally fragments alignments.

### 6.2 Affine Gap Penalty
The affine model decomposes gap costs into an initial **gap opening penalty** ($d$ or $g_o$) and an incremental **gap extension penalty** ($e$ or $g_e$):
$$W(k) = d + (k - 1) \cdot e \quad \text{for } k \ge 1, \quad \text{where } d > e > 0$$
$$\gamma_{\text{affine}}(k) = -W(k) = -d - (k - 1) \cdot e$$
- **Gap opening ($d$)**: High penalty reflecting the thermodynamic and evolutionary barrier to initiating an insertion/deletion.
- **Gap extension ($e$)**: Lower penalty reflecting that extending an established indel is structurally easier (e.g. extending an existing loop).

---

## 7. Gotoh (1982) Three-State Dynamic Programming for Affine Gaps

Naive dynamic programming with affine gap penalties requires evaluating all possible predecessor gap lengths $k \in [1, \max(m, n)]$ at every cell, resulting in $\mathcal{O}(m \cdot n \cdot \max(m, n)) \approx \mathcal{O}(n^3)$ computational time. **Gotoh (1982)** proved that affine penalties can be solved in optimal $\mathcal{O}(m \cdot n)$ time and space by decomposing the dynamic programming state space into three coupled matrices representing the alignment states.

### 7.1 Gap Penalty Convention and Sign Duality
In ProteinAlign Lab, gap penalties are parameterized as **positive penalty magnitudes** subtracted during score computation:
$$G(k) = d + (k - 1) \cdot e \quad (k \ge 1, \; d > 0, \; e > 0)$$
- $d = \text{gap\_open}$: Penalty incurred to initiate an indel of length 1.
- $e = \text{gap\_extend}$: Incremental penalty incurred for each additional residue beyond the first.

> [!NOTE]
> **Sign Convention Note**:
> In textbook literature and libraries like Biopython's `PairwiseAligner`, gap scores are often represented as negative values added to the score ($S + g$). The two conventions are mathematically isomorphic:
> $$\text{ProteinAlign Lab:} \quad S - d - (k - 1) \cdot e \iff \text{Biopython:} \quad S + g_{\text{open}} + (k - 1) \cdot g_{\text{extend}} \quad (g \le 0)$$
> In ProteinAlign Lab, users provide positive magnitudes ($d, e > 0$), and any negative input is strictly rejected with a descriptive `ValueError`.

### 7.2 State Space Definitions
At each cell $(i, j)$ for prefix $X_{1..i}$ and prefix $Y_{1..j}$:
- $M(i, j)$: Optimal alignment score given that residue $x_i$ is aligned with residue $y_j$ (match / mismatch state).
- $I_x(i, j)$: Optimal alignment score given that residue $x_i$ is aligned with a gap in $Y$ (deletion of $x_i$, represented by column $(x_i, \text{'-'})$).
- $I_y(i, j)$: Optimal alignment score given that residue $y_j$ is aligned with a gap in $X$ (insertion of $y_j$, represented by column $(\text{'-'}, y_j)$).

### 7.3 Recurrence Relations (Global Gotoh Alignment)
For $1 \le i \le m$ and $1 \le j \le n$:

$$M(i, j) = S(x_i, y_j) + \max \begin{cases}
  M(i-1, j-1) \\
  I_x(i-1, j-1) \\
  I_y(i-1, j-1)
\end{cases}$$

$$I_x(i, j) = \max \begin{cases}
  M(i-1, j) - d   & \text{[Open gap in } Y \text{ from match state]} \\
  I_x(i-1, j) - e & \text{[Extend existing gap in } Y \text{]}
\end{cases}$$

$$I_y(i, j) = \max \begin{cases}
  M(i, j-1) - d   & \text{[Open gap in } X \text{ from match state]} \\
  I_y(i, j-1) - e & \text{[Extend existing gap in } X \text{]}
\end{cases}$$

### 7.4 Justification of Opposing Gap Transition Exclusion
In Gotoh's classic formulation, direct transitions between $I_x$ and $I_y$ are omitted:
$$I_x(i, j) \not\leftarrow I_y(i-1, j) \quad \text{and} \quad I_y(i, j) \not\leftarrow I_x(i, j-1)$$
**Mathematical and Biological Justification**:
1. **Redundancy with Substitution**: Aligning residue $x_i$ to a gap $(x_i, \text{'-'})$ immediately adjacent to a gap aligning with $y_j$ $(\text{'-'}, y_j)$ represents a sequence mutation where $x_i$ and $y_j$ overlap. In an alignment matrix, this would consume two columns:
   $$\begin{pmatrix} x_i & \text{'-'} \\ \text{'-'} & y_j \end{pmatrix}$$
2. **Suboptimality**: Under affine penalties, creating an insertion adjacent to a deletion incurs two gap openings: $2 \cdot d$. Because $2d \gg 0$ (e.g. $2 \times 10 = 20$), this is virtually always strictly inferior to aligning $(x_i, y_j)$ directly as a mismatch in state $M$ with score $S(x_i, y_j)$ (typically $\ge -4$ in BLOSUM62).
3. **Elimination of Degenerate Tracebacks**: Disallowing $I_x \leftrightarrow I_y$ prevents the DP search space from generating arbitrary, biologically meaningless permutations of adjacent insertions and deletions.

### 7.5 Boundary Conditions and Impossible State Initialization
Let $-\infty$ denote an absorbing negative sentinel value ($-\infty = -10^9$) that cannot participate in any optimal path:
$$M(0, 0) = 0, \quad I_x(0, 0) = -\infty, \quad I_y(0, 0) = -\infty$$
For the boundary columns ($j = 0$, aligning prefix $X_{1..i}$ against empty $Y$):
$$M(i, 0) = -\infty, \quad I_y(i, 0) = -\infty, \quad I_x(i, 0) = -d - (i - 1) \cdot e \quad (1 \le i \le m)$$
For the boundary rows ($i = 0$, aligning empty $X$ against prefix $Y_{1..j}$):
$$M(0, j) = -\infty, \quad I_x(0, j) = -\infty, \quad I_y(0, j) = -d - (j - 1) \cdot e \quad (1 \le j \le n)$$

### 7.6 Optimal Global Score and State-Aware Traceback
- **Optimal Global Score**:
  $$\text{Score}_{\text{affine-global}} = \max\left( M(m, n), I_x(m, n), I_y(m, n) \right)$$
- **Tie-Breaking Priority at Termination**:
  $$\text{State Priority}: \quad M(m, n) > I_x(m, n) > I_y(m, n)$$
- **Traceback Invariants**:
  The traceback tracks the current DP layer $\text{state} \in \{M, I_x, I_y\}$.
  - When in state $M$: emit $(x_i, y_j)$, decrement $i \gets i - 1, j \gets j - 1$, and transition to whichever preceding state $\{M, I_x, I_y\}$ at $(i-1, j-1)$ achieved $M(i, j) - S(x_i, y_j)$.
  - When in state $I_x$: emit $(x_i, \text{'-'})$, decrement $i \gets i - 1$, and transition to $M$ if $I_x(i, j) == M(i-1, j) - d$, else remain in $I_x$.
  - When in state $I_y$: emit $(\text{'-'}, y_j)$, decrement $j \gets j - 1$, and transition to $M$ if $I_y(i, j) == M(i, j-1) - d$, else remain in $I_y$.
  - Boundary safeguard: If $i == 0$, only $I_y$ moves are valid; if $j == 0$, only $I_x$ moves are valid.

### 7.7 Affine Local Alignment (Smith–Waterman with Gotoh DP)
For local sequence alignment, all three layers incorporate the zero-restart cutoff:
$$M(i, j) = \max \begin{cases}
  0 \\
  S(x_i, y_j) + \max\left(0, M(i-1, j-1), I_x(i-1, j-1), I_y(i-1, j-1)\right)
\end{cases}$$
$$I_x(i, j) = \max \begin{cases}
  0 \\
  M(i-1, j) - d \\
  I_x(i-1, j) - e
\end{cases}$$
$$I_y(i, j) = \max \begin{cases}
  0 \\
  M(i, j-1) - d \\
  I_y(i, j-1) - e
\end{cases}$$

#### Theorem: Local Affine Optimal Alignment Terminates in State $M$
**Statement**: For any optimal local alignment under positive affine gap penalties ($d > 0, e > 0$), the peak score is attained in the match state:
$$\max_{\substack{0 \le i \le m \\ 0 \le j \le n}} \max(M(i, j), I_x(i, j), I_y(i, j)) = \max_{\substack{0 \le i \le m \\ 0 \le j \le n}} M(i, j)$$
**Proof**:
Suppose for contradiction that an optimal local alignment path terminates in a gap state, without loss of generality $I_x(i, j)$.
By the recurrence relation:
$$I_x(i, j) \le \max(M(i-1, j) - d, I_x(i-1, j) - e)$$
Since $d > 0$ and $e > 0$, any terminal gap column $(x_i, \text{'-'})$ strictly subtracts penalty $d$ or $e$ from the alignment without contributing any positive score. Truncating the trailing gap produces a shorter local sub-alignment ending at $(i-1, j)$ with score:
$$\text{Score}' \ge I_x(i, j) + \min(d, e) > I_x(i, j)$$
Thus, trailing gaps strictly degrade local alignment scores. The global maximum must therefore be achieved at a match/mismatch cell:
$$(i^*, j^*) = \arg\max_{i, j} M(i, j)$$
Traceback commences at $(i^*, j^*)$ in state $M$ and halts as soon as the score drops to $\le 0$.

### 7.8 Computational Complexity
- **Time Complexity**: Each of the $(m+1) \times (n+1)$ cells in three matrices requires $\mathcal{O}(1)$ operations (one addition and two or three comparisons). Total time complexity is strictly $\mathcal{O}(m \cdot n)$.
- **Space Complexity**: Storing the three score matrices $M, I_x, I_y$ requires $3 \times (m+1) \times (n+1) \times 8 \text{ bytes} \in \mathcal{O}(m \cdot n)$ storage. Traceback reconstruction operates in linear time $\mathcal{O}(m + n)$ on the filled DP matrices.

---

## 8. Position-Specific Scoring Matrices (PSSMs) and Profiles

### 8.1 Motivation: Context-Dependent Conservation
Standard substitution matrices (e.g. BLOSUM62, PAM250) treat amino acid mutability as spatially uniform: a tryptophan-to-arginine mismatch ($W \to R$) receives score $-3$ everywhere in every protein. In real proteins, however, evolutionary pressure is highly position-dependent:
- In catalytic active sites and the hydrophobic core, mutations disrupt stability or function and are eliminated by purifying selection.
- In solvent-exposed surface loops, mutations are neutral and accumulate freely.
A **Position-Specific Scoring Matrix (PSSM)** captures site-specific evolutionary preferences across a homologous protein family.

### 8.2 Profile Construction from Multiple Sequence Alignment
Let an alignment of $M$ homologous sequences over length $L$ be given. For column $i \in \{0, \dots, L-1\}$ and standard amino acid $a \in \Sigma$:

#### 1. Observed Residue Counts and Gap Policy
- $c_{i, a} = \sum_{k=1}^M \mathbb{I}(s_{k, i} == a)$
- $N_i = \sum_{a \in \Sigma} c_{i, a}$ (effective non-gap observation depth)
- $f_{\text{gap}, i} = \frac{\sum_{k=1}^M \mathbb{I}(s_{k, i} \in \{'-', '.'\})}{M}$

> [!IMPORTANT]
> **Gap Policy**: Gap characters are indel events, not amino acid mutations. They are recorded in $f_{\text{gap}, i}$ for structural analysis but excluded from the amino acid count denominator $N_i$. If an entire column consists of gaps ($N_i = 0$), the frequency distribution defaults smoothly to the background prior $q_a$.

#### 2. Pseudocount Regularization (Dirichlet Prior / Laplace Smoothing)
Small sample sizes ($M \ll 20$) lead to zero counts ($c_{i, a} = 0$), causing undefined log-odds ($\log(0) = -\infty$). Regularized frequencies are computed via:
$$P_{i, a} = \frac{c_{i, a} + \alpha \cdot q_a}{N_i + \alpha}$$
where:
- $\alpha \ge 0$ is the pseudocount weight (strength of the prior).
- $q_a > 0$ is the background proteome frequency ($\sum_{a \in \Sigma} q_a = 1.0$), by default Robinson & Robinson (1991).

Notice that for all $i$:
$$\sum_{a \in \Sigma} P_{i, a} = \frac{\sum_a c_{i, a} + \alpha \sum_a q_a}{N_i + \alpha} = \frac{N_i + \alpha}{N_i + \alpha} = 1.0$$

#### 3. Log-Odds Scoring Formula
The position-specific log-odds score $S_{i, a}$ evaluating residue $a$ at column $i$ is:
$$S_{i, a} = \lambda \cdot \log_b \left( \frac{P_{i, a}}{q_a} \right)$$
where:
- $b$ is the logarithm base ($b = 2.0$ for bits; $b = e$ for natural units).
- $\lambda$ is a scaling multiplier ($\lambda = 2.0$ yields half-bits; $\lambda = 3.0$ yields third-bits).

#### 4. Ambiguous Residue Scoring
For non-standard or ambiguous residues:
- `'B'` (Asx): $S_{i, B} = 0.5 \cdot S_{i, N} + 0.5 \cdot S_{i, D}$
- `'Z'` (Glx): $S_{i, Z} = 0.5 \cdot S_{i, Q} + 0.5 \cdot S_{i, E}$
- `'X'` (Unknown): Evaluated as the expected log-odds under the background distribution:
  $$S_{i, X} = \sum_{a \in \Sigma} q_a \cdot S_{i, a}$$

---

### 8.3 Worked Numerical Example (Viva Preparation)

Consider a toy alignment of $M = 3$ sequences at column $i = 0$:
- Observed residues: `H`, `H`, `H` ($c_{0, H} = 3$, all other $c_{0, a} = 0$).
- Background: Uniform distribution $q_a = \frac{1}{20} = 0.05$ for all $a \in \Sigma$.
- Pseudocount weight: $\alpha = 1.0$. Scaling: $\lambda = 2.0$, base $b = 2.0$ (half-bits).

#### Step 1: Effective Depth
$$N_0 = \sum_{a \in \Sigma} c_{0, a} = 3$$

#### Step 2: Regularized Probabilities
- For Histidine ($H$):
  $$P_{0, H} = \frac{c_{0, H} + \alpha \cdot q_H}{N_0 + \alpha} = \frac{3 + 1.0 \times 0.05}{3 + 1.0} = \frac{3.05}{4.0} = 0.7625$$
- For Alanine ($A$) or any unobserved residue:
  $$P_{0, A} = \frac{0 + 1.0 \times 0.05}{3 + 1.0} = \frac{0.05}{4.0} = 0.0125$$
- Verification of sum:
  $$\sum_{a \in \Sigma} P_{0, a} = 0.7625 + 19 \times 0.0125 = 0.7625 + 0.2375 = 1.000 \quad \checkmark$$

#### Step 3: Log-Odds Scores
- For Histidine ($H$):
  $$\frac{P_{0, H}}{q_H} = \frac{0.7625}{0.05} = 15.25$$
  $$S_{0, H} = 2.0 \cdot \log_2(15.25) = 2.0 \times 3.9307 = +7.861 \text{ half-bits}$$
- For Alanine ($A$):
  $$\frac{P_{0, A}}{q_A} = \frac{0.0125}{0.05} = 0.25$$
  $$S_{0, A} = 2.0 \cdot \log_2(0.25) = 2.0 \times (-2.0) = -4.000 \text{ half-bits}$$

*Explanation for Viva:* Conserved Histidine is rewarded with a strong positive score of $+7.86$ half-bits, while unobserved residues are penalized with $-4.00$ half-bits.

---

### 8.4 Profile-to-Sequence Dynamic Programming
When aligning an unaligned query sequence $Y = (y_1, \dots, y_n)$ against a PSSM profile of length $m$:
- Row dimension: Profile positions $i \in \{1, \dots, m\}$.
- Column dimension: Query residues $j \in \{1, \dots, n\}$.
- Column substitution score:
  $$s(i, y_j) = S_{i-1, y_j}$$

#### Global Profile Alignment (Gotoh Affine Gaps):
$$M(i, j) = S_{i-1, y_j} + \max\left( M(i-1, j-1), \, I_x(i-1, j-1), \, I_y(i-1, j-1) \right)$$
$$I_x(i, j) = \max\left( M(i-1, j) - g_{\text{open}}, \, I_x(i-1, j) - g_{\text{extend}} \right) \quad \text{[Profile deletion / gap in query]}$$
$$I_y(i, j) = \max\left( M(i, j-1) - g_{\text{open}}, \, I_y(i, j-1) - g_{\text{extend}} \right) \quad \text{[Query insertion / gap in profile]}$$

#### Local Profile Alignment (Smith–Waterman):
Floor all states at 0.0. The optimal local alignment begins at $(i^*, j^*) = \arg\max_{i, j} M(i, j)$ and terminates upon reaching 0.

#### Score Verification Invariant:
$$\text{Recomputed Score} = \sum_{\text{aligned columns}} S_{i, y_j} - \sum_{\text{gap runs}} \left( g_{\text{open}} + (k - 1) g_{\text{extend}} \right)$$
Must identically match dynamic programming terminal score within $10^{-4}$.

---

### 8.5 Comparison: Substitution Matrix vs. PSSM vs. Profile HMM

| Feature | Pairwise Matrix (BLOSUM62) | PSSM (ProteinAlign Lab) | Profile HMM (HMMER) |
| :--- | :--- | :--- | :--- |
| **Position Dependency** | Uniform (same for all positions) | Position-specific substitution | Position-specific substitution & position-specific gap rates |
| **Model Type** | Log-odds score lookup | Additive scoring profile matrix | Full probabilistic graphical model (state transitions) |
| **Gap Penalties** | Constant $d, e$ everywhere | Constant $d, e$ (or profile-guided) | Position-specific insertion/deletion probabilities |
| **Mathematical Engine** | Dynamic Programming (NW / SW) | Dynamic Programming (NW / SW) | Viterbi / Forward-Backward algorithms |
| **Sensitivity** | Good for $>30\%$ identity | Excellent for remote homologs ($15\text{--}30\%$) | Gold standard for remote homologs ($<20\%$) |
| **Complexity** | $\mathcal{O}(m \cdot n)$ | $\mathcal{O}(m \cdot n)$ | $\mathcal{O}(m \cdot n)$ with larger state machine |

---

## 9. Formal Assumptions, Conventions, and Mathematical Limitations

1. **Collinearity Assumption**: Sequences are assumed to be collinear; dynamic programming cannot detect inversions, circular permutations, or translocations.
2. **Context Independence**: The substitution score at position $k$ is assumed to be statistically independent of neighboring positions $k-1$ and $k+1$.
3. **Additive Scoring**: Total alignment score is strictly additive across columns.
4. **Symmetry of Substitution Matrices**: $S(a, b) = S(b, a)$ for standard pairwise matrices. Note that PSSM profile scoring is asymmetric because $P(i, b)$ is indexed by profile column position, not residue $a$.
5. **No Gap-to-Gap Scoring**: In Gotoh's model, transitions directly between $I_x$ and $I_y$ are disallowed ($I_x \leftrightarrow I_y$ forbidden) to prevent degenerate representations of substitution as insertion-deletion pairs.
