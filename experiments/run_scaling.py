"""Computational runtime and memory scaling experiment runner for RQ6.

Empirically benchmarks execution time (via time.perf_counter_ns) and heap memory
allocation (via tracemalloc) across sequence lengths L in [50, 100, 200, 300, 500, 750, 1000]
(length products m * n from 2,500 to 1,000,000) for all five core dynamic programming engines:
1. Needleman-Wunsch Global Linear
2. Smith-Waterman Local Linear
3. Gotoh Global Affine
4. Gotoh Local Affine
5. PSSM Profile-to-Sequence Global Affine

Usage:
    python experiments/run_scaling.py --config experiments/configs/scaling_config.json
"""

import argparse
import csv
import json
import time
import tracemalloc
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

from proteinalign import (
    BLOSUM62,
    PSSM,
    ROBINSON_FREQUENCIES,
    needleman_wunsch,
    needleman_wunsch_affine,
    smith_waterman,
    smith_waterman_affine,
)


OUTPUT_TABLES = Path("results/tables")


def generate_protein_sequence(length: int, seed: int) -> str:
    """Generates a synthetic protein sequence drawn from Robinson background frequencies."""
    rng = np.random.default_rng(seed)
    aas = list(ROBINSON_FREQUENCIES.keys())
    probs = [ROBINSON_FREQUENCIES[aa] for aa in aas]
    probs = [p / sum(probs) for p in probs]
    chosen = rng.choice(aas, size=length, p=probs)
    return "".join(chosen)


def generate_synthetic_profile(length: int, depth: int, seed: int) -> PSSM:
    """Constructs a synthetic PSSM of specified column length from synthetic homologous sequences."""
    aligned_seqs = [generate_protein_sequence(length, seed + i) for i in range(depth)]
    return PSSM.from_msa(aligned_seqs, name=f"synthetic_{length}", pseudocount_weight=1.0)


def measure_memory_peak(func, *args, **kwargs) -> float:
    """Measures peak heap allocation (KB) for a single execution using tracemalloc."""
    tracemalloc.start()
    try:
        func(*args, **kwargs)
        _, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()
    return round(peak / 1024.0, 2)


def run_scaling(config_path: Path) -> List[Dict[str, Any]]:
    """Runs scaling benchmarks across sequence lengths and algorithms."""
    with open(config_path, "r", encoding="utf-8") as f:
        cfg = json.load(f)

    lengths = cfg.get("sequence_lengths", [50, 100, 200, 300, 500, 750, 1000])
    repetitions = cfg.get("repetitions", 3)
    base_seed = cfg.get("random_seed", 42)

    records: List[Dict[str, Any]] = []

    print(f"Starting computational scaling experiments across lengths: {lengths}", flush=True)
    print(f"Repetitions per length: {repetitions} (memory measured separately to avoid timing distortion)", flush=True)

    for L in lengths:
        m, n = L, L
        len_prod = m * n
        print(f"\nBenchmarking length {L}x{L} (m*n = {len_prod:,})...", flush=True)

        # Baseline sequences for memory tracing
        seq1_mem = generate_protein_sequence(m, base_seed + L)
        seq2_mem = generate_protein_sequence(n, base_seed + L + 1)
        pssm_mem = generate_synthetic_profile(m, depth=5, seed=base_seed + L)

        # 1. Measure Peak Memory (isolated from timing)
        mem_nw = measure_memory_peak(needleman_wunsch, seq1_mem, seq2_mem, matrix=BLOSUM62, gap_penalty=8.0)
        mem_sw = measure_memory_peak(smith_waterman, seq1_mem, seq2_mem, matrix=BLOSUM62, gap_penalty=8.0)
        mem_ga = measure_memory_peak(needleman_wunsch_affine, seq1_mem, seq2_mem, matrix=BLOSUM62, gap_open=10.0, gap_extend=1.0)
        mem_la = measure_memory_peak(smith_waterman_affine, seq1_mem, seq2_mem, matrix=BLOSUM62, gap_open=10.0, gap_extend=1.0)
        mem_pa = measure_memory_peak(pssm_mem.align_global, seq2_mem, gap_open=10.0, gap_extend=1.0)

        # 2. Measure Pure Wall-Clock Runtime (no tracemalloc overhead)
        for rep in range(repetitions):
            seed1 = base_seed + L * 100 + rep * 2
            seed2 = base_seed + L * 100 + rep * 2 + 1
            seq1 = generate_protein_sequence(m, seed1)
            seq2 = generate_protein_sequence(n, seed2)
            pssm = generate_synthetic_profile(m, depth=5, seed=seed1)

            # NW Linear
            t0 = time.perf_counter_ns()
            res_nw = needleman_wunsch(seq1, seq2, matrix=BLOSUM62, gap_penalty=8.0)
            elapsed_nw = (time.perf_counter_ns() - t0) / 1e6

            records.append({
                "algorithm": "Needleman-Wunsch Linear",
                "length_m": m,
                "length_n": n,
                "length_product": len_prod,
                "repetition": rep,
                "runtime_ms": round(elapsed_nw, 3),
                "peak_memory_kb": mem_nw,
                "score": round(res_nw.score, 2),
            })

            # SW Linear
            t0 = time.perf_counter_ns()
            res_sw = smith_waterman(seq1, seq2, matrix=BLOSUM62, gap_penalty=8.0)
            elapsed_sw = (time.perf_counter_ns() - t0) / 1e6

            records.append({
                "algorithm": "Smith-Waterman Linear",
                "length_m": m,
                "length_n": n,
                "length_product": len_prod,
                "repetition": rep,
                "runtime_ms": round(elapsed_sw, 3),
                "peak_memory_kb": mem_sw,
                "score": round(res_sw.score, 2),
            })

            # Gotoh Global Affine
            t0 = time.perf_counter_ns()
            res_ga = needleman_wunsch_affine(seq1, seq2, matrix=BLOSUM62, gap_open=10.0, gap_extend=1.0)
            elapsed_ga = (time.perf_counter_ns() - t0) / 1e6

            records.append({
                "algorithm": "Gotoh Global Affine",
                "length_m": m,
                "length_n": n,
                "length_product": len_prod,
                "repetition": rep,
                "runtime_ms": round(elapsed_ga, 3),
                "peak_memory_kb": mem_ga,
                "score": round(res_ga.score, 2),
            })

            # Gotoh Local Affine
            t0 = time.perf_counter_ns()
            res_la = smith_waterman_affine(seq1, seq2, matrix=BLOSUM62, gap_open=10.0, gap_extend=1.0)
            elapsed_la = (time.perf_counter_ns() - t0) / 1e6

            records.append({
                "algorithm": "Gotoh Local Affine",
                "length_m": m,
                "length_n": n,
                "length_product": len_prod,
                "repetition": rep,
                "runtime_ms": round(elapsed_la, 3),
                "peak_memory_kb": mem_la,
                "score": round(res_la.score, 2),
            })

            # PSSM Global Affine
            t0 = time.perf_counter_ns()
            res_pa = pssm.align_global(seq2, gap_open=10.0, gap_extend=1.0)
            elapsed_pa = (time.perf_counter_ns() - t0) / 1e6

            records.append({
                "algorithm": "PSSM Global Affine",
                "length_m": m,
                "length_n": n,
                "length_product": len_prod,
                "repetition": rep,
                "runtime_ms": round(elapsed_pa, 3),
                "peak_memory_kb": mem_pa,
                "score": round(res_pa.score, 2),
            })

        print(f"L={L}: NW={elapsed_nw:.1f}ms, SW={elapsed_sw:.1f}ms, GotohAffine={elapsed_ga:.1f}ms, PSSM={elapsed_pa:.1f}ms | Peak Mem: NW={mem_nw}KB, GotohAffine={mem_ga}KB", flush=True)

    out_csv = Path(cfg.get("output_csv", "results/tables/scaling_results.csv"))
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "algorithm", "length_m", "length_n", "length_product",
        "repetition", "runtime_ms", "peak_memory_kb", "score"
    ]

    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)

    print(f"\nScaling benchmark completed! Saved {len(records)} records to {out_csv}", flush=True)
    return records


def main():
    parser = argparse.ArgumentParser(description="Run computational scaling benchmarks.")
    parser.add_argument("--config", type=str, default="experiments/configs/scaling_config.json")
    args = parser.parse_args()
    run_scaling(Path(args.config))


if __name__ == "__main__":
    main()
