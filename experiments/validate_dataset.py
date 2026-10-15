"""Dataset validation, quantitative auditing, and train/test split generation for MU-SABRE.

This script scans all reference alignments in data/raw/bench1.0/sabrem/ref, parses their
contents, verifies sequence lengths, counts valid amino acids, identifies core uppercase
residues, checks for malformed records, and establishes a stratified family-level
Train (60%), Validation (20%), and Test (20%) partitioning.

Usage:
    python experiments/validate_dataset.py
"""

from pathlib import Path
import json
import csv
from typing import Dict, List, Any
import numpy as np

from proteinalign.io import read_msa
from proteinalign.scoring import AMINO_ACIDS


REF_DIR = Path("data/raw/bench1.0/sabrem/ref")
INFO_DIR = Path("data/raw/bench1.0/sabrem/info")
OUTPUT_TABLES = Path("results/tables")


def validate_dataset() -> Dict[str, Any]:
    """Inspects and validates all files in the MU-SABRE benchmark dataset."""
    if not REF_DIR.exists():
        raise FileNotFoundError(f"Reference alignment directory not found: {REF_DIR}")

    files = sorted(list(REF_DIR.glob("*")))
    total_files = len(files)

    # Load percent identity metadata if present
    pctids: Dict[str, float] = {}
    pctids_file = INFO_DIR / "pctids.txt"
    if pctids_file.exists():
        for line in pctids_file.read_text(errors="replace").splitlines():
            line = line.strip()
            if line and ";" in line:
                fid, val = line.split(";", 1)
                try:
                    pctids[fid.strip()] = float(val.strip())
                except ValueError:
                    pass

    valid_files = 0
    excluded_files = []
    file_stats = []

    total_sequences = 0
    all_alignment_lengths = []
    all_seq_lengths = []
    all_core_fractions = []
    all_gap_fractions = []

    for file_path in files:
        if not file_path.is_file():
            continue

        fid = file_path.name
        try:
            msa = read_msa(file_path)
        except Exception as e:
            excluded_files.append({"file_id": fid, "reason": f"Parse error: {e}"})
            continue

        if not msa:
            excluded_files.append({"file_id": fid, "reason": "Empty MSA file"})
            continue

        num_seqs = len(msa)
        lengths = [len(seq) for seq in msa.values()]
        aln_len = lengths[0]

        # Verify equal length across all sequences in file
        if len(set(lengths)) > 1:
            excluded_files.append({"file_id": fid, "reason": f"Unequal lengths: {set(lengths)}"})
            continue

        # Inspect residue composition
        total_chars = 0
        core_chars = 0  # uppercase core
        gap_chars = 0
        invalid_symbols = set()

        for seq in msa.values():
            for ch in seq:
                total_chars += 1
                if ch in ("-", "."):
                    gap_chars += 1
                elif ch.isupper() and ch.upper() in AMINO_ACIDS:
                    core_chars += 1
                elif ch.islower() and ch.upper() in AMINO_ACIDS:
                    pass  # variable loop residue
                elif ch.upper() in ("B", "Z", "X"):
                    pass  # ambiguous residue
                else:
                    invalid_symbols.add(ch)

        if invalid_symbols:
            excluded_files.append({
                "file_id": fid,
                "reason": f"Invalid symbols detected: {invalid_symbols}"
            })
            continue

        valid_files += 1
        total_sequences += num_seqs
        all_alignment_lengths.append(aln_len)
        core_frac = core_chars / total_chars if total_chars > 0 else 0.0
        gap_frac = gap_chars / total_chars if total_chars > 0 else 0.0
        all_core_fractions.append(core_frac)
        all_gap_fractions.append(gap_frac)

        # Record pairwise count: K * (K - 1) / 2
        num_pairs = (num_seqs * (num_seqs - 1)) // 2

        file_stat = {
            "file_id": fid,
            "category": "superfamily" if fid.startswith("sup_") else "twilight",
            "num_sequences": num_seqs,
            "alignment_length": aln_len,
            "pairwise_combinations": num_pairs,
            "core_fraction": round(core_frac, 4),
            "gap_fraction": round(gap_frac, 4),
            "reference_pct_id": pctids.get(fid, -1.0),
        }
        file_stats.append(file_stat)

    # Stratified Train/Val/Test Split (60% Train, 20% Val, 20% Test)
    # Stratified by reference_pct_id tiers (<15%, 15-25%, >25%)
    np.random.seed(42)
    sorted_files = sorted(file_stats, key=lambda x: (x["reference_pct_id"], x["file_id"]))
    
    tier_low = [f for f in sorted_files if f["reference_pct_id"] < 15.0]
    tier_mid = [f for f in sorted_files if 15.0 <= f["reference_pct_id"] < 25.0]
    tier_high = [f for f in sorted_files if f["reference_pct_id"] >= 25.0]

    splits = {"train": [], "validation": [], "test": []}

    for tier in [tier_low, tier_mid, tier_high]:
        indices = np.random.permutation(len(tier))
        n_train = int(round(0.6 * len(tier)))
        n_val = int(round(0.2 * len(tier)))

        for i, idx in enumerate(indices):
            f_id = tier[idx]["file_id"]
            if i < n_train:
                splits["train"].append(f_id)
            elif i < n_train + n_val:
                splits["validation"].append(f_id)
            else:
                splits["test"].append(f_id)

    summary = {
        "dataset_name": "MU-SABRE (SABmark v1.65 re-aligned by MUSTANG)",
        "source_directory": str(REF_DIR),
        "compiler": "Dr. Robert C. Edgar (bench1.0/sabrem)",
        "total_files_examined": total_files,
        "valid_files_parsed": valid_files,
        "excluded_files_count": len(excluded_files),
        "excluded_files": excluded_files,
        "total_sequences": total_sequences,
        "total_pairwise_combinations": sum(f["pairwise_combinations"] for f in file_stats),
        "alignment_length_stats": {
            "min": int(np.min(all_alignment_lengths)),
            "median": float(np.median(all_alignment_lengths)),
            "mean": float(np.mean(all_alignment_lengths)),
            "max": int(np.max(all_alignment_lengths)),
        },
        "core_fraction_stats": {
            "mean": float(np.mean(all_core_fractions)),
            "median": float(np.median(all_core_fractions)),
        },
        "gap_fraction_stats": {
            "mean": float(np.mean(all_gap_fractions)),
            "median": float(np.median(all_gap_fractions)),
        },
        "split_counts": {
            "train": len(splits["train"]),
            "validation": len(splits["validation"]),
            "test": len(splits["test"]),
        },
        "splits": splits,
    }

    # Save outputs
    OUTPUT_TABLES.mkdir(parents=True, exist_ok=True)
    summary_json_path = OUTPUT_TABLES / "dataset_validation_summary.json"
    summary_json_path.write_text(json.dumps(summary, indent=2))

    summary_csv_path = OUTPUT_TABLES / "dataset_files_catalog.csv"
    with open(summary_csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "file_id", "category", "num_sequences", "alignment_length",
            "pairwise_combinations", "core_fraction", "gap_fraction", "reference_pct_id"
        ])
        writer.writeheader()
        writer.writerows(file_stats)

    print(f"Dataset validation complete. Parsed {valid_files}/{total_files} files.")
    print(f"Total sequences: {total_sequences} | Total pairwise combinations: {summary['total_pairwise_combinations']}")
    print(f"Stratified split: {len(splits['train'])} train, {len(splits['validation'])} val, {len(splits['test'])} test.")
    print(f"Saved: {summary_json_path} and {summary_csv_path}")

    return summary


if __name__ == "__main__":
    validate_dataset()
