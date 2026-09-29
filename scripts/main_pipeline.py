#!/usr/bin/env python3
"""
main_pipeline.py
================
End-to-end pipeline for the MS Percolation Framework.

Steps
-----
1.  Load / generate data.
2.  Fit CT-MSM  →  Q matrix.
3.  Build directed weighted graph G(Q).
4.  Run threshold sweep  →  PercolationCurve.
5.  Compute patient-level δ(t) series.
6.  Fit percolation-augmented MSM (PA-MSM).
7.  Compare baseline vs PA-MSM (AIC/BIC/LR).
8.  Cluster patient trajectories by percolation features.
9.  Save results (CSV / HDF5).
10. Generate figures.

Usage
-----
python scripts/main_pipeline.py --config configs/default.yaml [--data path/to/data.csv]

"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

# Make src/ importable when running from project root

import numpy as np
import pandas as pd

from ms_synth.utils import setup_logging, get_logger, set_seed, load_config, ensure_dir, save_dict_json
from ms_synth.data import generate_synthetic_data, validate_data, load_data
from ms_synth.msm_fit import fit_msm, Q_to_transition_matrix
from ms_synth.graph_utils import build_graph, build_theta_grid, normalise_weights
from ms_synth.percolation import (
    threshold_sweep,
    compute_patient_delta_series,
    compute_summary_statistics,
    percolation_vulnerability_index,
)
from ms_synth.augmentation import (
    augment_Q,
    compare_models,
    cluster_trajectories_by_percolation,
    fit_augmented_msm,
)
from ms_synth.visualization import (
    plot_percolation_curve,
    plot_component_sizes,
    plot_susceptibility,
    plot_network_snapshot,
    plot_patient_delta,
    plot_model_comparison,
    plot_summary_dashboard,
)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="MS Percolation Framework – End-to-End Pipeline",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--config", "-c",
        default="configs/default.yaml",
        help="Path to YAML / JSON configuration file.",
    )
    parser.add_argument(
        "--data", "-d",
        default=None,
        help="Path to longitudinal MS data CSV (if None, synthetic data is used).",
    )
    parser.add_argument(
        "--output", "-o",
        default="results",
        help="Output directory for results and figures.",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Override random seed from config.",
    )
    parser.add_argument(
        "--n-patients",
        type=int,
        default=None,
        help="Override number of synthetic patients from config.",
    )
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    args = parse_args()

    # ── Logging ──────────────────────────────────────────────────────────────
    out_dir = Path(args.output)
    ensure_dir(out_dir)
    log_path = out_dir / "pipeline.log"
    setup_logging(level=args.log_level, log_file=log_path)
    logger = get_logger(__name__)
    logger.info("=" * 60)
    logger.info("MS Percolation Framework – pipeline start")
    logger.info("=" * 60)

    # ── Config ────────────────────────────────────────────────────────────────
    cfg = load_config(args.config)
    seed = args.seed if args.seed is not None else cfg.get("reproducibility", {}).get("seed", 42)
    set_seed(seed)
    logger.info("Using seed=%d", seed)

    n_states = cfg["states"]["n_states"]
    state_labels = {int(k): v for k, v in cfg["states"]["labels"].items()}
    starting_states = cfg["states"]["starting_states"]

    figs_dir = ensure_dir(out_dir / "figures")

    # ── Step 1: Data ──────────────────────────────────────────────────────────
    if args.data is not None:
        logger.info("Loading data from %s", args.data)
        df = load_data(args.data, cfg=cfg)
    else:
        syn_cfg = cfg["data"]["synthetic"]
        n_patients = args.n_patients or syn_cfg["n_patients"]
        logger.info("Generating synthetic data: n_patients=%d", n_patients)
        df = generate_synthetic_data(
            n_patients=n_patients,
            n_states=n_states,
            max_time=syn_cfg["max_time"],
            time_step=syn_cfg["time_step"],
            seed=seed,
        )
        df.to_csv(out_dir / "synthetic_data.csv", index=False)
        logger.info("Synthetic data saved to %s/synthetic_data.csv", out_dir)

    logger.info("Data: %d rows, %d patients, %d states",
                len(df), df["patient_id"].nunique(), n_states)

    # ── Step 2: Fit CT-MSM ────────────────────────────────────────────────────
    logger.info("Fitting CT-MSM (method=%s)…", cfg["msm"]["method"])
    Q = fit_msm(df, n_states=n_states, cfg=cfg)

    logger.info("Q matrix (diagonal):\n%s", np.diag(Q))
    np.savetxt(out_dir / "Q_matrix.csv", Q, delimiter=",")
    logger.info("Q matrix saved.")

    # ── Step 3: Build directed graph ──────────────────────────────────────────
    logger.info("Building directed weighted graph from Q…")
    G = build_graph(Q, state_labels=state_labels)
    logger.info("Graph: %d nodes, %d edges", G.number_of_nodes(), G.number_of_edges())

    # ── Step 4: Threshold sweep (percolation analysis) ────────────────────────
    logger.info("Running threshold sweep…")
    curve = threshold_sweep(G, cfg=cfg, starting_states=starting_states)

    # Save curve to CSV
    curve_df = curve.to_dataframe()
    curve_df.to_csv(out_dir / "percolation_curve.csv", index=False)
    logger.info("Percolation curve saved.")

    # Summary statistics
    stats = compute_summary_statistics(curve)
    save_dict_json(stats, out_dir / "percolation_summary.json")
    logger.info("Percolation summary: %s", stats)

    # ── Step 5: Patient-level δ(t) series ────────────────────────────────────
    logger.info("Computing patient-level δ(t) series…")
    # Use a smaller subset for speed in the demo (all patients for production)
    demo_patients = df["patient_id"].unique()[:50]
    df_demo = df[df["patient_id"].isin(demo_patients)].copy()
    patient_delta_df = compute_patient_delta_series(
        df_demo, Q_pop=Q, curve_pop=curve, cfg=cfg
    )
    patient_delta_df.to_csv(out_dir / "patient_delta_series.csv", index=False)
    logger.info("Patient δ(t) series saved.")

    # ── Step 6: Trajectory clustering ────────────────────────────────────────
    logger.info("Clustering patient trajectories by percolation features…")
    n_clusters = 3
    patient_delta_df = cluster_trajectories_by_percolation(
        patient_delta_df, n_clusters=n_clusters, seed=seed
    )
    logger.info(
        "Cluster distribution:\n%s",
        patient_delta_df.groupby("patient_id")["trajectory_cluster"]
        .first()
        .value_counts()
        .to_string(),
    )

    # ── Step 7: Augmented MSM ─────────────────────────────────────────────────
    logger.info("Fitting percolation-augmented MSM (PA-MSM)…")
    Q_aug, aug_metrics = fit_augmented_msm(Q, df_demo, curve, cfg=cfg)
    np.savetxt(out_dir / "Q_augmented.csv", Q_aug, delimiter=",")
    save_dict_json(aug_metrics, out_dir / "pa_msm_metrics.json")
    logger.info("PA-MSM metrics: %s", aug_metrics)

    # ── Step 8: Model comparison ──────────────────────────────────────────────
    logger.info("Comparing baseline vs PA-MSM…")
    comparison = compare_models(Q, Q_aug, df_demo, n_extra_params=1)
    save_dict_json(comparison, out_dir / "model_comparison.json")
    logger.info(
        "Model comparison: preferred=%s  ΔAIC=%.2f  p=%.4f",
        comparison["preferred"],
        comparison["delta_aic"],
        comparison["lr_pvalue"],
    )

    # ── Step 9: Figures ───────────────────────────────────────────────────────
    logger.info("Generating figures…")

    dpi = cfg.get("output", {}).get("dpi", 150)

    plot_percolation_curve(curve, save_path=figs_dir / "percolation_curve.png", dpi=dpi)
    plot_component_sizes(curve, save_path=figs_dir / "component_sizes.png", dpi=dpi)
    plot_susceptibility(curve, save_path=figs_dir / "susceptibility.png", dpi=dpi)

    plot_network_snapshot(
        G, theta=0.5 * curve.theta_c, state_labels=state_labels,
        save_path=figs_dir / "network_subcritical.png", dpi=dpi,
    )
    plot_network_snapshot(
        G, theta=1.5 * curve.theta_c, state_labels=state_labels,
        save_path=figs_dir / "network_supercritical.png", dpi=dpi,
    )

    plot_patient_delta(
        patient_delta_df,
        theta_c=curve.theta_c,
        save_path=figs_dir / "patient_delta.png",
        dpi=dpi,
    )
    plot_model_comparison(comparison, save_path=figs_dir / "model_comparison.png", dpi=dpi)

    plot_summary_dashboard(
        curve, G,
        patient_delta_df=patient_delta_df,
        comparison_result=comparison,
        save_path=figs_dir / "summary_dashboard.png",
        dpi=dpi,
    )

    logger.info("All figures saved to %s", figs_dir)

    # ── Final summary ─────────────────────────────────────────────────────────
    logger.info("=" * 60)
    logger.info("Pipeline complete.")
    logger.info("  θ_c            = %.4f", curve.theta_c)
    logger.info("  I_perc         = %.4f", curve.percolation_integral)
    logger.info("  ΔAIC           = %.2f  (preferred: %s)", comparison["delta_aic"], comparison["preferred"])
    logger.info("  Results in:    %s", out_dir.resolve())
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
