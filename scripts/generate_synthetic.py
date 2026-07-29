#!/usr/bin/env python3
"""
generate_synthetic.py
=====================
Command-line entry point for generating a synthetic longitudinal MS
dataset using the ms_percolation_framework synthetic data module.

Usage examples
--------------
# Minimal (uses configs/synthetic.yaml):
python scripts/generate_synthetic.py

# Custom output path and patient count:
python scripts/generate_synthetic.py \\
    --config configs/synthetic.yaml \\
    --output data/synthetic_ms_3000.parquet \\
    --n-patients 3000

# Force CSV output:
python scripts/generate_synthetic.py \\
    --config configs/synthetic.yaml \\
    --output data/synthetic_ms.csv \\
    --format csv

# Generate a small validation subset (fast):
python scripts/generate_synthetic.py \\
    --config configs/synthetic.yaml \\
    --output data/synthetic_small.parquet \\
    --n-patients 200 --seed 99

# Validate against generating Q and print statistics:
python scripts/generate_synthetic.py \\
    --config configs/synthetic.yaml \\
    --output data/synth.parquet \\
    --validate
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Make src/ importable when run from project root
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.utils import setup_logging, get_logger, load_config, ensure_dir
from src.synthetic_data import (
    generate_synthetic_ms_dataset,
    build_Q_from_config,
    compare_empirical_to_Q,
)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate a synthetic longitudinal MS dataset.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--config", "-c",
        default="configs/synthetic.yaml",
        help="Path to synthetic data YAML / JSON config.",
    )
    parser.add_argument(
        "--output", "-o",
        default="data/synthetic_ms_longitudinal.parquet",
        help="Output file path (.parquet or .csv).",
    )
    parser.add_argument(
        "--n-patients", "-n",
        type=int,
        default=None,
        help="Number of patients (overrides config value).",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed (overrides config value).",
    )
    parser.add_argument(
        "--format",
        choices=["parquet", "csv"],
        default=None,
        help="Output format.  Inferred from extension if not set.",
    )
    parser.add_argument(
        "--validate",
        action="store_true",
        help="After generation, compare empirical Q to generating Q and print results.",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
    )
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    args = parse_args()
    setup_logging(level=args.log_level)
    logger = get_logger(__name__)

    # ── Load config ───────────────────────────────────────────────────────
    cfg = load_config(args.config)

    # ── Generate dataset ──────────────────────────────────────────────────
    df = generate_synthetic_ms_dataset(
        cfg=cfg,
        n_patients=args.n_patients,
        seed=args.seed,
    )

    # ── Save output ───────────────────────────────────────────────────────
    out_path = Path(args.output)
    ensure_dir(out_path.parent)

    # Determine format
    fmt = args.format
    if fmt is None:
        ext = out_path.suffix.lower()
        fmt = "csv" if ext == ".csv" else "parquet"

    if fmt == "parquet":
        df.to_parquet(out_path, index=False)
    else:
        df.to_csv(out_path, index=False)

    logger.info(
        "Dataset saved: %s  (%d rows, %d patients, %.1f MB)",
        out_path,
        len(df),
        df["patient_id"].nunique(),
        out_path.stat().st_size / 1e6,
    )

    # ── Optional: empirical Q validation ─────────────────────────────────
    if args.validate:
        logger.info("Running Q comparison (empirical vs generating Q)…")
        n_states = int(cfg["states"]["n_states"])
        Q_true = build_Q_from_config(cfg)
        comparison_df = compare_empirical_to_Q(df, Q_true)
        comparison_df = comparison_df.sort_values("relative_error", ascending=False)

        print("\n── Empirical Q vs Generating Q (top 10 by error) ──────────────")
        print(comparison_df.head(10).to_string(index=False))
        print()
        mean_rel_err = float(comparison_df["relative_error"].mean())
        logger.info("Mean relative error across non-zero Q entries: %.3f", mean_rel_err)

        val_path = out_path.parent / (out_path.stem + "_Q_validation.csv")
        comparison_df.to_csv(val_path, index=False)
        logger.info("Q validation saved: %s", val_path)


if __name__ == "__main__":
    main()
