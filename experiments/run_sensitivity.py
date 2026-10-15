"""Parameter sensitivity analysis runner for RQ5.

Explores dynamic programming behavior and alignment fidelity across:
1. 2D grid of gap-opening (6.0 - 16.0) and gap-extension (0.5 - 3.0) penalties under Gotoh DP.
2. 1D sweep of pseudocount smoothing weight alpha (0.1 - 10.0) in PSSM construction.

Usage:
    python experiments/run_sensitivity.py --config experiments/configs/sensitivity_config.json
"""

import argparse
import csv
import json
import time
from pathlib import Path
from typing import Any, Dict, List

from proteinalign import (
    BLOSUM62,
    PSSM,
    compute_alignment_gap_stats,
    compute_pairwise_metrics,
    compute_profile_placement_accuracy,
    needleman_wunsch_affine,
    read_msa,
    recompute_score,
    split_msa_train_test,
)


REF_DIR = Path("data/raw/bench1.0/sabrem/ref")
OUTPUT_TABLES = Path("results/tables")


def run_sensitivity(config_path: Path) -> List[Dict[str, Any]]:
    """Executes sensitivity grid search over gap penalties and pseudocount weights."""
    with open(config_path, "r", encoding="utf-8") as f:
        cfg = json.load(f)

    target_files = cfg.get("target_files", ["sup_003", "sup_005", "sup_008", "twi_028", "twi_050"])
    max_pairs = cfg.get("max_pairs_per_file", 3)
    gap_open_grid = cfg.get("gap_open_grid", [6.0, 8.0, 10.0, 12.0, 14.0, 16.0])
    gap_extend_grid = cfg.get("gap_extend_grid", [0.5, 1.0, 1.5, 2.0, 3.0])
    alpha_grid = cfg.get("pssm_alpha_grid", [0.1, 0.5, 1.0, 2.0, 5.0, 10.0])
    def_open = cfg.get("default_gap_open", 10.0)
    def_ext = cfg.get("default_gap_extend", 1.0)

    records: List[Dict[str, Any]] = []

    print(f"Starting parameter sensitivity analysis on {len(target_files)} benchmark families...")
    print(f"Gap penalty grid: {len(gap_open_grid)} open values x {len(gap_extend_grid)} extend values")
    print(f"PSSM alpha sweep: {len(alpha_grid)} values")

    # Load sequence pairs and MSAs
    loaded_data = []
    for fid in target_files:
        fpath = REF_DIR / fid
        if not fpath.exists():
            continue
        msa = read_msa(fpath)
        keys = list(msa.keys())
        if len(keys) < 2:
            continue

        pairs = []
        for i in range(len(keys)):
            for j in range(i + 1, len(keys)):
                pairs.append((keys[i], keys[j]))
        if len(pairs) > max_pairs:
            pairs = pairs[:max_pairs]

        loaded_data.append((fid, msa, pairs))

    # -------------------------------------------------------------
    # 1. 2D Gap Penalty Grid Search (Gotoh Affine)
    # -------------------------------------------------------------
    total_grid_points = len(gap_open_grid) * len(gap_extend_grid)
    grid_idx = 0

    for g_open in gap_open_grid:
        for g_ext in gap_extend_grid:
            grid_idx += 1
            if g_open < g_ext:
                continue

            for fid, msa, pairs in loaded_data:
                category = "superfamily" if fid.startswith("sup_") else "twilight"
                for k1, k2 in pairs:
                    ref1, ref2 = msa[k1], msa[k2]
                    s1 = ref1.replace("-", "").replace(".", "")
                    s2 = ref2.replace("-", "").replace(".", "")

                    t0 = time.perf_counter_ns()
                    res = needleman_wunsch_affine(s1, s2, matrix=BLOSUM62, gap_open=g_open, gap_extend=g_ext)
                    runtime_ms = (time.perf_counter_ns() - t0) / 1e6

                    rec_score = recompute_score(res.seq1_aligned, res.seq2_aligned, matrix=BLOSUM62, gap_open=g_open, gap_extend=g_ext)
                    metrics = compute_pairwise_metrics(res.seq1_aligned, res.seq2_aligned, ref1, ref2, uppercase_only=True)

                    records.append({
                        "experiment": "gap_penalty_grid",
                        "file_id": fid,
                        "category": category,
                        "pair_id": f"{k1}_vs_{k2}",
                        "param_gap_open": g_open,
                        "param_gap_extend": g_ext,
                        "param_alpha": "N/A",
                        "score": round(res.score, 2),
                        "recomputed_score": round(rec_score, 2),
                        "q_score": metrics["q_score"],
                        "modeler_score": metrics["modeler_score"],
                        "f1_score": metrics["f1_score"],
                        "pid": metrics["pid"],
                        "total_gaps": metrics["total_gaps"],
                        "gap_openings": metrics["gap_openings"],
                        "mean_gap_length": metrics["mean_gap_length"],
                        "runtime_ms": round(runtime_ms, 3),
                    })

    print(f"Gap penalty grid search completed. Records: {len(records)}")

    # -------------------------------------------------------------
    # 2. 1D PSSM Pseudocount Alpha Sweep
    # -------------------------------------------------------------
    for alpha in alpha_grid:
        for fid, msa, _ in loaded_data:
            keys = list(msa.keys())
            if len(keys) < 3:
                continue
            category = "superfamily" if fid.startswith("sup_") else "twilight"
            # LOO on first 2 sequences
            for target_id in keys[:2]:
                train_msa, _ = split_msa_train_test(msa, target_id)
                ref_target = msa[target_id]
                target_seq = ref_target.replace("-", "").replace(".", "")

                pssm = PSSM.from_msa(train_msa, name=f"{fid}_loo_{target_id}", pseudocount_weight=alpha)
                t0 = time.perf_counter_ns()
                res = pssm.align_global(target_seq, gap_open=def_open, gap_extend=def_ext, query_id=target_id)
                runtime_ms = (time.perf_counter_ns() - t0) / 1e6

                acc = compute_profile_placement_accuracy(res.seq1_aligned, res.seq2_aligned, ref_target)
                gap_stats = compute_alignment_gap_stats(res.seq1_aligned, res.seq2_aligned)

                records.append({
                    "experiment": "pssm_alpha_sweep",
                    "file_id": fid,
                    "category": category,
                    "pair_id": f"Profile_vs_{target_id}",
                    "param_gap_open": def_open,
                    "param_gap_extend": def_ext,
                    "param_alpha": alpha,
                    "score": round(res.score, 2),
                    "recomputed_score": round(res.score, 2),
                    "q_score": acc["profile_q_score"],
                    "modeler_score": acc["profile_q_score"],
                    "f1_score": acc["profile_q_score"],
                    "pid": round(res.identity, 2),
                    "total_gaps": gap_stats["total_gaps"],
                    "gap_openings": gap_stats["gap_openings"],
                    "mean_gap_length": gap_stats["mean_gap_length"],
                    "runtime_ms": round(runtime_ms, 3),
                })

    out_csv = Path(cfg.get("output_csv", "results/tables/sensitivity_results.csv"))
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "experiment", "file_id", "category", "pair_id",
        "param_gap_open", "param_gap_extend", "param_alpha",
        "score", "recomputed_score", "q_score", "modeler_score", "f1_score",
        "pid", "total_gaps", "gap_openings", "mean_gap_length", "runtime_ms"
    ]

    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)

    print(f"Sensitivity experiments completed successfully! Saved {len(records)} records to {out_csv}")
    return records


def main():
    parser = argparse.ArgumentParser(description="Run parameter sensitivity experiments.")
    parser.add_argument("--config", type=str, default="experiments/configs/sensitivity_config.json")
    args = parser.parse_args()
    run_sensitivity(Path(args.config))


if __name__ == "__main__":
    main()
