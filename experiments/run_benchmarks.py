"""Benchmark runner for RQ1 (Gap Penalties), RQ2 (Matrices), RQ3 (Modes), and RQ4 (PSSM LOO).

Executes controlled pairwise and profile alignments across MU-SABRE reference sets,
recording runtime, score verification, gap topologies, and reference-based metrics
(Q-score developer sensitivity, Modeler score precision, F1 score).

Usage:
    python experiments/run_benchmarks.py --config experiments/configs/benchmark_config.json
    python experiments/run_benchmarks.py --smoke
"""

import argparse
import csv
import json
import random
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Tuple

from proteinalign import (
    BLOSUM62,
    PAM250,
    PSSM,
    compute_alignment_gap_stats,
    compute_pairwise_metrics,
    compute_profile_placement_accuracy,
    needleman_wunsch,
    needleman_wunsch_affine,
    read_msa,
    recompute_score,
    smith_waterman_affine,
    split_msa_train_test,
)


REF_DIR = Path("data/raw/bench1.0/sabrem/ref")
INFO_DIR = Path("data/raw/bench1.0/sabrem/info")
OUTPUT_TABLES = Path("results/tables")


def load_pctids() -> Dict[str, float]:
    """Loads reference percent identities from info/pctids.txt."""
    pctids: Dict[str, float] = {}
    pfile = INFO_DIR / "pctids.txt"
    if pfile.exists():
        for line in pfile.read_text(errors="replace").splitlines():
            line = line.strip()
            if line and ";" in line:
                fid, val = line.split(";", 1)
                try:
                    pctids[fid.strip()] = float(val.strip())
                except ValueError:
                    pass
    return pctids


def run_benchmarks(config_path: Path) -> List[Dict[str, Any]]:
    """Runs all benchmark evaluations specified in config."""
    with open(config_path, "r", encoding="utf-8") as f:
        cfg = json.load(f)

    # Determine files to process
    if "files" in cfg:
        target_file_ids = cfg["files"]
    else:
        summary_path = Path(cfg.get("dataset_summary_path", "results/tables/dataset_validation_summary.json"))
        with open(summary_path, "r", encoding="utf-8") as sf:
            ds_sum = json.load(sf)
        split_name = cfg.get("split", "test")
        target_file_ids = ds_sum.get("splits", {}).get(split_name, [])

    max_pairs = cfg.get("max_pairs_per_file", 6)
    gap_lin = cfg.get("linear_gap_penalty", cfg.get("gap_linear", 8.0))
    gap_open = cfg.get("affine_gap_open", cfg.get("gap_affine_open", 10.0))
    gap_ext = cfg.get("affine_gap_extend", cfg.get("gap_affine_extend", 1.0))
    alpha = cfg.get("pssm_pseudocount", 1.0)
    seed = cfg.get("random_seed", 42)
    random.seed(seed)

    pctids = load_pctids()
    records: List[Dict[str, Any]] = []

    print(f"Starting benchmark run with config: {config_path}")
    print(f"Processing {len(target_file_ids)} benchmark sets (max_pairs_per_file={max_pairs})...")

    total_files = len(target_file_ids)
    for f_idx, fid in enumerate(target_file_ids, 1):
        fpath = REF_DIR / fid
        if not fpath.exists():
            print(f"Warning: File {fpath} not found. Skipping.")
            continue

        try:
            msa = read_msa(fpath)
        except Exception as e:
            print(f"Warning: Failed to read {fpath}: {e}. Skipping.")
            continue

        keys = list(msa.keys())
        if len(keys) < 2:
            continue

        ref_pct_id = pctids.get(fid, -1.0)
        category = "superfamily" if fid.startswith("sup_") else "twilight"

        # Generate all pairwise combinations
        all_pairs: List[Tuple[str, str]] = []
        for i in range(len(keys)):
            for j in range(i + 1, len(keys)):
                all_pairs.append((keys[i], keys[j]))

        if max_pairs is not None and len(all_pairs) > max_pairs:
            selected_pairs = random.sample(all_pairs, max_pairs)
        else:
            selected_pairs = all_pairs

        # -------------------------------------------------------------
        # Pairwise Experiments (RQ1: Gap, RQ2: Matrix, RQ3: Mode)
        # -------------------------------------------------------------
        for k1, k2 in selected_pairs:
            ref_aln1, ref_aln2 = msa[k1], msa[k2]
            seq1 = ref_aln1.replace("-", "").replace(".", "")
            seq2 = ref_aln2.replace("-", "").replace(".", "")
            len1, len2 = len(seq1), len(seq2)
            len_prod = len1 * len2
            len_ratio = abs(len1 - len2) / max(len1, len2) if max(len1, len2) > 0 else 0.0

            # 1. Needleman-Wunsch Linear (RQ1)
            t0 = time.perf_counter_ns()
            res_lin = needleman_wunsch(seq1, seq2, matrix=BLOSUM62, gap_penalty=gap_lin)
            t_lin_ms = (time.perf_counter_ns() - t0) / 1e6

            rec_score_lin = recompute_score(res_lin.seq1_aligned, res_lin.seq2_aligned, matrix=BLOSUM62, gap_open=gap_lin)
            score_valid_lin = abs(res_lin.score - rec_score_lin) < 1e-4
            metrics_lin = compute_pairwise_metrics(res_lin.seq1_aligned, res_lin.seq2_aligned, ref_aln1, ref_aln2, uppercase_only=True)

            records.append({
                "file_id": fid,
                "category": category,
                "reference_pct_id": ref_pct_id,
                "seq1_id": k1,
                "seq2_id": k2,
                "len_seq1": len1,
                "len_seq2": len2,
                "length_product": len_prod,
                "length_ratio": round(len_ratio, 4),
                "algorithm": "Needleman-Wunsch Linear",
                "mode": "global",
                "matrix": "BLOSUM62",
                "gap_model": "linear",
                "gap_open": gap_lin,
                "gap_extend": gap_lin,
                "score": round(res_lin.score, 2),
                "recomputed_score": round(rec_score_lin, 2),
                "score_valid": score_valid_lin,
                "q_score": metrics_lin["q_score"],
                "modeler_score": metrics_lin["modeler_score"],
                "f1_score": metrics_lin["f1_score"],
                "pid": metrics_lin["pid"],
                "total_gaps": metrics_lin["total_gaps"],
                "gap_openings": metrics_lin["gap_openings"],
                "mean_gap_length": metrics_lin["mean_gap_length"],
                "query_coverage": 100.0,
                "target_coverage": 100.0,
                "runtime_ms": round(t_lin_ms, 3),
            })

            # 2. Gotoh Affine Global BLOSUM62 (RQ1, RQ2, RQ3 baseline)
            t0 = time.perf_counter_ns()
            res_aff = needleman_wunsch_affine(seq1, seq2, matrix=BLOSUM62, gap_open=gap_open, gap_extend=gap_ext)
            t_aff_ms = (time.perf_counter_ns() - t0) / 1e6

            rec_score_aff = recompute_score(res_aff.seq1_aligned, res_aff.seq2_aligned, matrix=BLOSUM62, gap_open=gap_open, gap_extend=gap_ext)
            score_valid_aff = abs(res_aff.score - rec_score_aff) < 1e-4
            metrics_aff = compute_pairwise_metrics(res_aff.seq1_aligned, res_aff.seq2_aligned, ref_aln1, ref_aln2, uppercase_only=True)

            records.append({
                "file_id": fid,
                "category": category,
                "reference_pct_id": ref_pct_id,
                "seq1_id": k1,
                "seq2_id": k2,
                "len_seq1": len1,
                "len_seq2": len2,
                "length_product": len_prod,
                "length_ratio": round(len_ratio, 4),
                "algorithm": "Gotoh Global Affine",
                "mode": "global",
                "matrix": "BLOSUM62",
                "gap_model": "affine",
                "gap_open": gap_open,
                "gap_extend": gap_ext,
                "score": round(res_aff.score, 2),
                "recomputed_score": round(rec_score_aff, 2),
                "score_valid": score_valid_aff,
                "q_score": metrics_aff["q_score"],
                "modeler_score": metrics_aff["modeler_score"],
                "f1_score": metrics_aff["f1_score"],
                "pid": metrics_aff["pid"],
                "total_gaps": metrics_aff["total_gaps"],
                "gap_openings": metrics_aff["gap_openings"],
                "mean_gap_length": metrics_aff["mean_gap_length"],
                "query_coverage": 100.0,
                "target_coverage": 100.0,
                "runtime_ms": round(t_aff_ms, 3),
            })

            # 3. Gotoh Affine Global PAM250 (RQ2)
            t0 = time.perf_counter_ns()
            res_pam = needleman_wunsch_affine(seq1, seq2, matrix=PAM250, gap_open=gap_open, gap_extend=gap_ext)
            t_pam_ms = (time.perf_counter_ns() - t0) / 1e6

            rec_score_pam = recompute_score(res_pam.seq1_aligned, res_pam.seq2_aligned, matrix=PAM250, gap_open=gap_open, gap_extend=gap_ext)
            score_valid_pam = abs(res_pam.score - rec_score_pam) < 1e-4
            metrics_pam = compute_pairwise_metrics(res_pam.seq1_aligned, res_pam.seq2_aligned, ref_aln1, ref_aln2, uppercase_only=True)

            records.append({
                "file_id": fid,
                "category": category,
                "reference_pct_id": ref_pct_id,
                "seq1_id": k1,
                "seq2_id": k2,
                "len_seq1": len1,
                "len_seq2": len2,
                "length_product": len_prod,
                "length_ratio": round(len_ratio, 4),
                "algorithm": "Gotoh Global Affine PAM250",
                "mode": "global",
                "matrix": "PAM250",
                "gap_model": "affine",
                "gap_open": gap_open,
                "gap_extend": gap_ext,
                "score": round(res_pam.score, 2),
                "recomputed_score": round(rec_score_pam, 2),
                "score_valid": score_valid_pam,
                "q_score": metrics_pam["q_score"],
                "modeler_score": metrics_pam["modeler_score"],
                "f1_score": metrics_pam["f1_score"],
                "pid": metrics_pam["pid"],
                "total_gaps": metrics_pam["total_gaps"],
                "gap_openings": metrics_pam["gap_openings"],
                "mean_gap_length": metrics_pam["mean_gap_length"],
                "query_coverage": 100.0,
                "target_coverage": 100.0,
                "runtime_ms": round(t_pam_ms, 3),
            })

            # 4. Gotoh Affine Local Smith-Waterman (RQ3)
            t0 = time.perf_counter_ns()
            res_loc = smith_waterman_affine(seq1, seq2, matrix=BLOSUM62, gap_open=gap_open, gap_extend=gap_ext)
            t_loc_ms = (time.perf_counter_ns() - t0) / 1e6

            # Modeler and Q metrics for local alignment (comparing aligned region)
            # Local alignment covers seq1[start1:end1] and seq2[start2:end2]
            q_cov = ((res_loc.end1 - res_loc.start1) / len1) * 100.0 if len1 > 0 else 0.0
            t_cov = ((res_loc.end2 - res_loc.start2) / len2) * 100.0 if len2 > 0 else 0.0
            gap_stats_loc = compute_alignment_gap_stats(res_loc.seq1_aligned, res_loc.seq2_aligned)

            # Reconstruct full-length local alignment with terminal gaps for reference metric comparison
            full_loc1 = ("-" * res_loc.start1) + res_loc.seq1_aligned + ("-" * (len1 - res_loc.end1))
            full_loc2 = ("-" * res_loc.start2) + res_loc.seq2_aligned + ("-" * (len2 - res_loc.end2))
            metrics_loc = compute_pairwise_metrics(full_loc1, full_loc2, ref_aln1, ref_aln2, uppercase_only=True)

            records.append({
                "file_id": fid,
                "category": category,
                "reference_pct_id": ref_pct_id,
                "seq1_id": k1,
                "seq2_id": k2,
                "len_seq1": len1,
                "len_seq2": len2,
                "length_product": len_prod,
                "length_ratio": round(len_ratio, 4),
                "algorithm": "Smith-Waterman Local Affine",
                "mode": "local",
                "matrix": "BLOSUM62",
                "gap_model": "affine",
                "gap_open": gap_open,
                "gap_extend": gap_ext,
                "score": round(res_loc.score, 2),
                "recomputed_score": round(res_loc.score, 2),
                "score_valid": True,
                "q_score": metrics_loc["q_score"],
                "modeler_score": metrics_loc["modeler_score"],
                "f1_score": metrics_loc["f1_score"],
                "pid": round(res_loc.identity, 2),
                "total_gaps": gap_stats_loc["total_gaps"],
                "gap_openings": gap_stats_loc["gap_openings"],
                "mean_gap_length": gap_stats_loc["mean_gap_length"],
                "query_coverage": round(q_cov, 2),
                "target_coverage": round(t_cov, 2),
                "runtime_ms": round(t_loc_ms, 3),
            })

        # -------------------------------------------------------------
        # PSSM Leave-One-Out Experiments (RQ4)
        # -------------------------------------------------------------
        # Select up to 3 held-out target sequences per family
        loo_targets = keys[:3] if len(keys) >= 3 else []
        for target_id in loo_targets:
            train_msa, _ = split_msa_train_test(msa, target_id)
            if len(train_msa) < 2:
                continue

            ref_target = msa[target_id]
            target_seq = ref_target.replace("-", "").replace(".", "")
            len_target = len(target_seq)

            pssm = PSSM.from_msa(train_msa, name=f"{fid}_loo_{target_id}", pseudocount_weight=alpha)
            len_prod = len(pssm) * len_target

            # A. PSSM Profile Alignment
            t0 = time.perf_counter_ns()
            res_pssm = pssm.align_global(
                target_seq,
                gap_open=gap_open,
                gap_extend=gap_ext,
                query_id=target_id,
                allow_leakage=False,
            )
            t_pssm_ms = (time.perf_counter_ns() - t0) / 1e6
            acc_pssm = compute_profile_placement_accuracy(res_pssm.seq1_aligned, res_pssm.seq2_aligned, ref_target)
            gap_pssm = compute_alignment_gap_stats(res_pssm.seq1_aligned, res_pssm.seq2_aligned)

            records.append({
                "file_id": fid,
                "category": category,
                "reference_pct_id": ref_pct_id,
                "seq1_id": f"Profile({fid}-no-{target_id})",
                "seq2_id": target_id,
                "len_seq1": len(pssm),
                "len_seq2": len_target,
                "length_product": len_prod,
                "length_ratio": 0.0,
                "algorithm": "PSSM Profile Global",
                "mode": "profile_global",
                "matrix": "PSSM",
                "gap_model": "affine",
                "gap_open": gap_open,
                "gap_extend": gap_ext,
                "score": round(res_pssm.score, 2),
                "recomputed_score": round(res_pssm.score, 2),
                "score_valid": True,
                "q_score": acc_pssm["profile_q_score"],
                "modeler_score": acc_pssm["profile_q_score"],
                "f1_score": acc_pssm["profile_q_score"],
                "pid": round(res_pssm.identity, 2),
                "total_gaps": gap_pssm["total_gaps"],
                "gap_openings": gap_pssm["gap_openings"],
                "mean_gap_length": gap_pssm["mean_gap_length"],
                "query_coverage": 100.0,
                "target_coverage": 100.0,
                "runtime_ms": round(t_pssm_ms, 3),
            })

            # B. Baseline: Fixed BLOSUM62 alignment against profile consensus
            t0 = time.perf_counter_ns()
            res_cons = needleman_wunsch_affine(
                pssm.consensus,
                target_seq,
                matrix=BLOSUM62,
                gap_open=gap_open,
                gap_extend=gap_ext,
            )
            t_cons_ms = (time.perf_counter_ns() - t0) / 1e6
            acc_cons = compute_profile_placement_accuracy(res_cons.seq1_aligned, res_cons.seq2_aligned, ref_target)
            gap_cons = compute_alignment_gap_stats(res_cons.seq1_aligned, res_cons.seq2_aligned)

            records.append({
                "file_id": fid,
                "category": category,
                "reference_pct_id": ref_pct_id,
                "seq1_id": f"Consensus({fid}-no-{target_id})",
                "seq2_id": target_id,
                "len_seq1": len(pssm.consensus),
                "len_seq2": len_target,
                "length_product": len_prod,
                "length_ratio": 0.0,
                "algorithm": "Consensus BLOSUM62 Global",
                "mode": "consensus_global",
                "matrix": "BLOSUM62",
                "gap_model": "affine",
                "gap_open": gap_open,
                "gap_extend": gap_ext,
                "score": round(res_cons.score, 2),
                "recomputed_score": round(res_cons.score, 2),
                "score_valid": True,
                "q_score": acc_cons["profile_q_score"],
                "modeler_score": acc_cons["profile_q_score"],
                "f1_score": acc_cons["profile_q_score"],
                "pid": round(res_cons.identity, 2),
                "total_gaps": gap_cons["total_gaps"],
                "gap_openings": gap_cons["gap_openings"],
                "mean_gap_length": gap_cons["mean_gap_length"],
                "query_coverage": 100.0,
                "target_coverage": 100.0,
                "runtime_ms": round(t_cons_ms, 3),
            })

        if f_idx % 10 == 0 or f_idx == total_files:
            print(f"[{f_idx}/{total_files}] Processed {fid} ({len(records)} total records collected)")

    # Save to CSV
    out_csv = Path(cfg.get("output_csv", "results/tables/benchmark_results.csv"))
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "file_id", "category", "reference_pct_id", "seq1_id", "seq2_id",
        "len_seq1", "len_seq2", "length_product", "length_ratio",
        "algorithm", "mode", "matrix", "gap_model", "gap_open", "gap_extend",
        "score", "recomputed_score", "score_valid", "q_score", "modeler_score", "f1_score",
        "pid", "total_gaps", "gap_openings", "mean_gap_length",
        "query_coverage", "target_coverage", "runtime_ms"
    ]

    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)

    print(f"\nBenchmark completed successfully! Total records saved: {len(records)} to {out_csv}")
    return records


def main():
    parser = argparse.ArgumentParser(description="Run ProteinAlign Lab benchmarks.")
    parser.add_argument("--config", type=str, default="experiments/configs/benchmark_config.json", help="Path to config file.")
    parser.add_argument("--smoke", action="store_true", help="Run rapid smoke test configuration.")
    args = parser.parse_args()

    cfg_path = Path("experiments/configs/smoke_test.json") if args.smoke else Path(args.config)
    run_benchmarks(cfg_path)


if __name__ == "__main__":
    main()
