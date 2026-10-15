"""ProteinAlign Lab — Interactive Demonstration Script.
Run this script directly:
    python demo.py
"""

from proteinalign import (
    needleman_wunsch,
    smith_waterman,
    needleman_wunsch_affine,
    BLOSUM62,
    PAM250,
    SubstitutionMatrix,
    format_alignment,
    build_profile_from_msa,
    split_msa_train_test,
    PSSM,
    DataLeakageError,
)


def run_demo():
    print("=" * 70)
    print("      ProteinAlign Lab — Core Sequence Alignment Demo")
    print("=" * 70)

    # -------------------------------------------------------------
    # 1. Needleman–Wunsch (Global Alignment) with Simple Matrix
    # -------------------------------------------------------------
    print("\n[1] Global Alignment (Needleman–Wunsch) — Hand-Calculated Toy Example")
    print("-" * 70)
    seq_a = "HEAG"
    seq_b = "HEAP"
    simple_mat = SubstitutionMatrix.create_identity(match=2.0, mismatch=-1.0)
    res_simple = needleman_wunsch(seq_a, seq_b, matrix=simple_mat, gap_penalty=3.0)

    print(f"Query Sequence : {seq_a}")
    print(f"Target Sequence: {seq_b}")
    print(f"Matrix         : {simple_mat.name} (Match=+2, Mismatch=-1, Gap=3)")
    print(f"Optimal Score  : {res_simple.score}  (Expected: 3 matches * 2 - 1 mismatch * 1 = 5.0)")
    print(f"Aligned Seq 1  : {res_simple.seq1_aligned}")
    print(f"Aligned Seq 2  : {res_simple.seq2_aligned}")
    print(f"Percent Identity: {res_simple.identity:.1f}%")

    # -------------------------------------------------------------
    # 2. Needleman–Wunsch with BLOSUM62 vs PAM250
    # -------------------------------------------------------------
    print("\n[2] Global Alignment — Comparing BLOSUM62 vs PAM250")
    print("-" * 70)
    query = "HEAGAWGHEE"
    target = "PAWHEAE"

    res_blosum = needleman_wunsch(query, target, matrix=BLOSUM62, gap_penalty=8.0)
    res_pam = needleman_wunsch(query, target, matrix=PAM250, gap_penalty=8.0)

    print(f"Sequences: '{query}' vs '{target}'\n")
    print(f"BLOSUM62 Score: {res_blosum.score:.1f} | Identity: {res_blosum.identity:.1f}%")
    print(format_alignment(res_blosum))

    print(f"\nPAM250 Score  : {res_pam.score:.1f} | Identity: {res_pam.identity:.1f}%")
    print(format_alignment(res_pam))

    # -------------------------------------------------------------
    # 3. Smith–Waterman (Local Alignment)
    # -------------------------------------------------------------
    print("\n[3] Local Alignment (Smith–Waterman) — Extracting Conserved Motif")
    print("-" * 70)
    s1 = "DDDDPAWHEAERRRR"
    s2 = "WWWWPAWHEAEWWWW"

    res_local = smith_waterman(s1, s2, matrix=BLOSUM62, gap_penalty=8.0)
    print(f"Sequence 1: {s1}")
    print(f"Sequence 2: {s2}")
    print(f"Extracted Motif: '{res_local.seq1_aligned}'")
    print(f"Local Score    : {res_local.score:.1f}")
    print(f"Seq 1 Slices   : [{res_local.start1}:{res_local.end1}] -> '{s1[res_local.start1:res_local.end1]}'")
    print(f"Seq 2 Slices   : [{res_local.start2}:{res_local.end2}] -> '{s2[res_local.start2:res_local.end2]}'")

    # -------------------------------------------------------------
    # 4. Affine Gap Penalties (Gotoh 1982) vs Linear Gap Model
    # -------------------------------------------------------------
    print("\n[4] Affine Gap Alignment (Gotoh 1982) — Contiguous Indel Grouping")
    print("-" * 70)
    prot1 = "MKFLILLVLFPAAD"
    prot2 = "MKFLILLFNILCLFPVLAAD"

    res_lin = needleman_wunsch(prot1, prot2, matrix=BLOSUM62, gap_penalty=8.0)
    res_aff = needleman_wunsch_affine(prot1, prot2, matrix=BLOSUM62, gap_open=10.0, gap_extend=1.0)

    print(f"Protein 1: {prot1}")
    print(f"Protein 2: {prot2}\n")
    print(f"Linear Gap (d=8.0) — Score: {res_lin.score:.1f}")
    print(f"  Seq 1: {res_lin.seq1_aligned}")
    print(f"  Seq 2: {res_lin.seq2_aligned}")
    print("  Note: Linear model fragments indels into multiple scattered 2-residue gaps (--).")

    print(f"\nAffine Gap (d=10.0, e=1.0) — Score: {res_aff.score:.1f}")
    print(f"  Seq 1: {res_aff.seq1_aligned}")
    print(f"  Seq 2: {res_aff.seq2_aligned}")
    print("  Note: Affine model groups indels into a contiguous block (----) reflecting biological loop deletion.")

    # -------------------------------------------------------------
    # 5. Position-Specific Scoring Matrices (PSSM) & Profile Alignment
    # -------------------------------------------------------------
    print("\n[5] Position-Specific Scoring Matrix (PSSM) & Profile Alignment")
    print("-" * 70)
    print("  [Note: The toy MSA below is illustrative and demonstrates mathematical mechanics]")

    toy_msa = {
        "homolog_1": "HEAGAWG",
        "homolog_2": "HEAPAWG",
        "homolog_3": "HEAGAWA",
        "target_X":  "HEAPAWG",
    }

    # 1. Scientific Rigor: Split out held-out test sequence to prevent data leakage
    train_msa, held_out_seq = split_msa_train_test(toy_msa, test_id="target_X")
    print(f"Training MSA Homologs : {list(train_msa.keys())}")
    print(f"Held-Out Query Target : 'target_X' -> '{held_out_seq}'")

    # 2. Build Profile and PSSM with pseudocount regularization
    profile = build_profile_from_msa(train_msa, name="demo_family", pseudocount_weight=1.0)
    pssm = PSSM.from_profile(profile, scale=2.0)
    print(f"Profile Consensus     : {profile.consensus} (Length: {profile.alignment_length})")
    print(f"PSSM Score for 'H' at pos 0 (conserved)  : +{pssm.score(0, 'H'):.2f} half-bits")
    print(f"PSSM Score for 'D' at pos 0 (unobserved) : {pssm.score(0, 'D'):.2f} half-bits")

    # 3. Align held-out query against profile
    res_pssm = pssm.align_global(held_out_seq, gap_open=10.0, gap_extend=1.0, query_id="target_X")
    recomputed = pssm.recompute_score(res_pssm)

    print(f"\nProfile-to-Sequence Global Score: {res_pssm.score:.2f}")
    print(f"  Profile : {res_pssm.seq1_aligned}")
    print(f"  Query   : {res_pssm.seq2_aligned}")
    print(f"  Independent Score Recomputation : {recomputed:.2f} (Verified match: {abs(recomputed - res_pssm.score) < 1e-4})")

    # 4. Demonstrate Data Leakage Safeguard
    try:
        pssm.align_global(held_out_seq, query_id="homolog_1")
    except DataLeakageError:
        print("  Data Leakage Safeguard: Successfully prevented evaluating training ID 'homolog_1'!")

    print("\n" + "=" * 70)
    print("Demo completed successfully! To run automated tests, use: pytest")
    print("=" * 70)


if __name__ == "__main__":
    run_demo()
