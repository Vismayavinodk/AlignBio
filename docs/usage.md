# ProteinAlign Lab — Usage Guide and Code Examples

This guide provides simple, step-by-step usage examples demonstrating how to use the `proteinalign` package for global and local protein sequence alignment.

---

## 1. Global Alignment with Basic Match/Mismatch Scoring

Use `SubstitutionMatrix.create_identity()` to define a simple reward and penalty scheme (e.g. +2 for matches, -1 for mismatches):

```python
from proteinalign import needleman_wunsch, SubstitutionMatrix

# Define two short protein sequences
seq1 = "HEAG"
seq2 = "HEAP"

# Create a simple scoring matrix: +2 for match, -1 for mismatch
simple_matrix = SubstitutionMatrix.create_identity(match=2.0, mismatch=-1.0)

# Run Needleman-Wunsch global alignment with a linear gap penalty of 3.0
result = needleman_wunsch(seq1, seq2, matrix=simple_matrix, gap_penalty=3.0)

# Inspect the results
print("Alignment Type:", result.alignment_type)
print("Optimal Score :", result.score)
print("Aligned Seq 1 :", result.seq1_aligned)
print("Aligned Seq 2 :", result.seq2_aligned)
print(f"Sequence Identity: {result.identity:.1f}%")
```

**Expected Output:**
```text
Alignment Type: global
Optimal Score : 5.0
Aligned Seq 1 : HEAG
Aligned Seq 2 : HEAP
Sequence Identity: 75.0%
```

*Explanation for Viva:* 3 matches (H, E, A) give $3 \times 2 = 6$. 1 mismatch (G vs P) gives $-1$. Total score is $6 - 1 = 5.0$.

---

## 2. Local Alignment with Basic Scoring

Use `smith_waterman` to locate a conserved motif embedded inside divergent or unrelated flanking sequences:

```python
from proteinalign import smith_waterman, SubstitutionMatrix

# Sequence with conserved motif "HEAG" flanked by dissimilar residues
seq1 = "DDHEAGRR"
seq2 = "WWHEAGWW"

# Simple matrix: match = +3, mismatch = -2
simple_matrix = SubstitutionMatrix.create_identity(match=3.0, mismatch=-2.0)

# Run Smith-Waterman local alignment
result = smith_waterman(seq1, seq2, matrix=simple_matrix, gap_penalty=4.0)

print("Optimal Local Score:", result.score)
print("Aligned Motif 1    :", result.seq1_aligned)
print("Aligned Motif 2    :", result.seq2_aligned)
print(f"Seq 1 Coordinates  : [{result.start1}:{result.end1}] -> '{result.seq1[result.start1:result.end1]}'")
print(f"Seq 2 Coordinates  : [{result.start2}:{result.end2}] -> '{result.seq2[result.start2:result.end2]}'")
```

**Expected Output:**
```text
Optimal Local Score: 12.0
Aligned Motif 1    : HEAG
Aligned Motif 2    : HEAG
Seq 1 Coordinates  : [2:6] -> 'HEAG'
Seq 2 Coordinates  : [2:6] -> 'HEAG'
```

*Explanation for Viva:* The algorithm zeroes out the negative flanking residues (`DD` vs `WW` and `RR` vs `WW`) and isolates only the positive segment `HEAG` with score $4 \times 3 = 12.0$.

---

## 3. Global Alignment with BLOSUM62

Standard NCBI BLOSUM62 is the default matrix for protein alignment:

```python
from proteinalign import needleman_wunsch, BLOSUM62

seq1 = "HEAGAWGHEE"
seq2 = "PAWHEAE"

# Align using BLOSUM62 and linear gap penalty of 8.0
result = needleman_wunsch(seq1, seq2, matrix=BLOSUM62, gap_penalty=8.0)

print(f"Score: {result.score:.1f}")
print("Query :", result.seq1_aligned)
print("Target:", result.seq2_aligned)
print(f"Identity: {result.identity:.1f}%")
```

---

## 4. Local Alignment with PAM250

Use Dayhoff's PAM250 matrix when comparing divergent sequences separated by deep evolutionary time:

```python
from proteinalign import smith_waterman, PAM250

seq1 = "MKFLILLFNILCLFPVLAAD"
seq2 = "MKFLILLVLFPAAD"

result = smith_waterman(seq1, seq2, matrix=PAM250, gap_penalty=8.0)

print(f"PAM250 Local Score: {result.score:.1f}")
print("Seq 1:", result.seq1_aligned)
print("Seq 2:", result.seq2_aligned)
print(f"Coordinates: Seq1 [{result.start1}:{result.end1}], Seq2 [{result.start2}:{result.end2}]")
```

---

## 5. Affine Gap Alignment with Gotoh (1982) 3-State DP

Under an affine gap model, opening an insertion or deletion costs `gap_open` ($d$), while extending an established gap costs `gap_extend` ($e$). A contiguous gap of length $k \ge 1$ incurs:
$$\text{Cost}(k) = \text{gap\_open} + (k - 1) \cdot \text{gap\_extend}$$

### 5.1 Global Affine Alignment (`needleman_wunsch_affine`)

```python
from proteinalign import needleman_wunsch, needleman_wunsch_affine, BLOSUM62

seq1 = "MKFLILLVLFPAAD"
seq2 = "MKFLILLFNILCLFPVLAAD"

# Linear model (gap_penalty = 8.0)
res_linear = needleman_wunsch(seq1, seq2, matrix=BLOSUM62, gap_penalty=8.0)

# Gotoh Affine model (gap_open = 10.0, gap_extend = 1.0)
res_affine = needleman_wunsch_affine(seq1, seq2, matrix=BLOSUM62, gap_open=10.0, gap_extend=1.0)

print(f"Linear Score : {res_linear.score:.1f}")
print("Linear Seq1  :", res_linear.seq1_aligned)
print("Linear Seq2  :", res_linear.seq2_aligned)
# Note: Linear model fragments into multiple scattered gaps: MKFLILL--V--LFP--AAD

print(f"\nAffine Score : {res_affine.score:.1f}")
print("Affine Seq1  :", res_affine.seq1_aligned)
print("Affine Seq2  :", res_affine.seq2_aligned)
# Note: Affine model groups indels into a contiguous block: MKFLIL----LVLFP--AAD
```

### 5.2 Local Affine Alignment (`smith_waterman_affine`)

```python
from proteinalign import smith_waterman_affine, BLOSUM62

query = "XXXXHEAGWWWWWWXXXX"
target = "YYYYHEAWAGWWWWWWYYYY"

# Local alignment with affine gap penalty
result = smith_waterman_affine(query, target, matrix=BLOSUM62, gap_open=4.0, gap_extend=1.0)

print(f"Local Affine Score: {result.score:.1f}")
print("Extracted Query   :", result.seq1_aligned)
print("Extracted Target  :", result.seq2_aligned)
print(f"Query Slice [{result.start1}:{result.end1}] : '{query[result.start1:result.end1]}'")
print(f"Target Slice [{result.start2}:{result.end2}]: '{target[result.start2:result.end2]}'")
```

---

## 7. Position-Specific Scoring Matrices (PSSMs) and Profile Alignment

A PSSM captures position-specific evolutionary conservation across a homologous protein family:

### 7.1 Constructing a Profile and PSSM from an MSA

```python
from proteinalign import build_profile_from_msa, PSSM

# Define a multiple sequence alignment (MSA) of related proteins
msa = {
    "prot_A": "HEAGAWG",
    "prot_B": "HEAPAWG",
    "prot_C": "HEAGAWA",
}

# 1. Build a validated frequency Profile with pseudocount regularization
profile = build_profile_from_msa(msa, name="test_family", pseudocount_weight=1.0)
print(f"Profile Length : {profile.alignment_length}")
print(f"Num Sequences  : {profile.num_sequences}")
print(f"Consensus      : {profile.consensus}")

# 2. Build a PSSM with half-bit scaling (scale = 2.0, base = 2.0)
pssm = PSSM.from_profile(profile, scale=2.0)
print(f"Score for 'H' at position 0: {pssm.score(0, 'H'):.2f} half-bits")
print(f"Score for 'D' at position 0: {pssm.score(0, 'D'):.2f} half-bits")
```

### 7.2 Profile-to-Sequence Global Alignment

```python
# Align a separate held-out query sequence against the profile
query_seq = "HEAPAWG"
result = pssm.align_global(query_seq, gap_open=10.0, gap_extend=1.0)

print(f"Global Alignment Score: {result.score:.1f}")
print("Profile Consensus     :", result.seq1_aligned)
print("Aligned Query         :", result.seq2_aligned)

# Verify score by independent column-wise summation
assert pssm.recompute_score(result) == result.score
```

### 7.3 Data Leakage Safeguards (Held-Out Evaluation)

```python
from proteinalign import split_msa_train_test, DataLeakageError

# Perform a valid leave-one-out train/test split on an MSA
train_msa, test_query = split_msa_train_test(msa, test_id="prot_A")

# Train profile strictly on remaining training sequences
train_pssm = PSSM.from_msa(train_msa, name="trained_profile")

# Attempting to align a training sequence triggers an error unless allow_leakage=True
try:
    train_pssm.align_global("HEAGAWG", query_id="prot_B")
except DataLeakageError as e:
    print("Caught expected leakage error:", e)

# Aligning the true held-out query succeeds cleanly
res_held_out = train_pssm.align_global(test_query, query_id="prot_A")
print(f"Held-out score: {res_held_out.score:.1f}")
```

---

## 8. Visualizing Alignments

You can format alignment outputs for clear terminal display using `format_alignment`:

```python
from proteinalign import needleman_wunsch, format_alignment, BLOSUM62

result = needleman_wunsch("HEAGAWGHEE", "PAWHEAE", matrix=BLOSUM62, gap_penalty=8.0)
print(format_alignment(result, line_width=60))
```

---

## 9. How to Run the Unit Tests

The test suite uses `pytest` and covers normal cases, edge cases, invariants, and independent validation against Biopython (127 automated tests):

```bash
# Run all 127 tests
pytest -v

# Run only profile tests
pytest tests/test_profiles.py -v

# Run only PSSM tests
pytest tests/test_pssm.py -v

# Run only affine gap tests
pytest tests/test_affine_gaps.py -v

# Run Biopython cross-validation tests
pytest tests/test_biopython_compat.py -v
```
