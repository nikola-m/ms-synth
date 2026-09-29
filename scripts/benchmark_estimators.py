#!/usr/bin/env python3
"""
Known-truth benchmark of CT-MSM generator-matrix estimators (ADEMP design).

Aims
    Quantify how well three estimators recover the generating intensities
    when the latent disease process is observed only at irregular visits.
Data-generating mechanisms
    A  cohort size N in {250, 500, 1000, 2000, 4000} x {homogeneous,
       heterogeneous} population, default visit schedule;
    B  homogeneous population, N = 1000, visit intervals scaled by
       {0.5, 1, 2} (mean 3.5/4.5, 7/9, 14/18 months in RRMS/SPMS).
    "Homogeneous" switches off latent-class, age, sex and DMT scaling so
    that every patient's generator equals Q_base exactly.
Estimands
    The 26 non-zero off-diagonal intensities of Q_base.
Methods
    oracle  N_ij/T_i from the exact latent paths (upper bound);
    crude   N_ij/T_i from consecutive visits (common shortcut);
    panel   panel-likelihood MLE, P(dt) = expm(Q dt), structural zeros known.
Performance measures
    Relative bias and relative absolute error per intensity, summarised
    over intensities and replicates (Monte Carlo SE reported downstream).

Output: results/benchmark_estimators.csv (one row per estimate).
Usage:  python scripts/benchmark_estimators.py [--reps 5] [--jobs 4] [--quick]
"""
from __future__ import annotations

import argparse
import copy
import logging
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from ms_synth.msm_fit import estimate_Q_crude, estimate_Q_from_paths, estimate_Q_panel
from ms_synth.synthetic_data import build_Q_from_config, generate_synthetic_ms_dataset
from ms_synth.utils import load_config

ROOT = Path(__file__).resolve().parent.parent
CONFIG = ROOT / "configs" / "synthetic.yaml"
OUT = ROOT / "results" / "benchmark_estimators.csv"
BASE_SEED = 1000


def make_config(homogeneous: bool, interval_mult: float) -> dict:
    cfg = copy.deepcopy(load_config(str(CONFIG)))
    if homogeneous:
        cfg["cohort"]["class_fractions"] = {"stable": 0.0, "moderate": 1.0, "aggressive": 0.0}
        ase = cfg["age_sex_effects"]
        ase["age_progression_multiplier_per_decade"] = 1.0
        ase["male_progression_multiplier"] = 1.0
        cfg["demographics"]["dmt_status_probs"] = {"none": 1.0, "moderate_dmt": 0.0, "high_dmt": 0.0}
    v = cfg["visits"]
    v["mean_interval_rrms_months"] *= interval_mult
    v["mean_interval_spms_months"] *= interval_mult
    return cfg


def run_task(task: dict) -> list[dict]:
    logging.disable(logging.CRITICAL)
    cfg = make_config(task["homogeneous"], task["interval_mult"])
    Q = build_Q_from_config(cfg)
    nz = [(i, j) for i in range(12) for j in range(12) if i != j and Q[i, j] > 0]
    t0 = time.time()
    df, latent = generate_synthetic_ms_dataset(
        cfg=cfg, n_patients=task["n"], seed=task["seed"], return_latent=True)
    t_gen = time.time() - t0
    est = {}
    t0 = time.time(); est["oracle"] = estimate_Q_from_paths(latent, 12); t_or = time.time() - t0
    t0 = time.time(); est["crude"] = estimate_Q_crude(df, 12); t_cr = time.time() - t0
    t0 = time.time(); est["panel"] = estimate_Q_panel(df, 12, allowed=nz); t_pa = time.time() - t0
    times = {"oracle": t_or, "crude": t_cr, "panel": t_pa}
    rows = []
    for name, Qh in est.items():
        for (i, j) in nz:
            rows.append({
                "scenario": task["scenario"], "population":
                "homogeneous" if task["homogeneous"] else "heterogeneous",
                "n_patients": task["n"], "interval_mult": task["interval_mult"],
                "rep": task["rep"], "seed": task["seed"], "estimator": name,
                "from_state": i, "to_state": j, "q_true": float(Q[i, j]),
                "q_hat": float(Qh[i, j]), "n_visits": len(df),
                "fit_seconds": round(times[name], 3), "gen_seconds": round(t_gen, 3),
            })
    return rows


def build_tasks(reps: int, quick: bool) -> list[dict]:
    sizes = [250, 500, 1000] if quick else [250, 500, 1000, 2000, 4000]
    tasks = []
    for r in range(reps):
        seed = BASE_SEED + r
        for hom in (True, False):
            for n in sizes:
                tasks.append(dict(scenario="A_size", homogeneous=hom, n=n,
                                  interval_mult=1.0, rep=r, seed=seed))
        for m in (0.5, 2.0):
            tasks.append(dict(scenario="B_visits", homogeneous=True, n=1000,
                              interval_mult=m, rep=r, seed=seed))
    # biggest first for better load balancing
    return sorted(tasks, key=lambda t: -t["n"])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reps", type=int, default=5)
    ap.add_argument("--jobs", type=int, default=1)
    ap.add_argument("--quick", action="store_true", help="small sizes for CI")
    ap.add_argument("--out", default=str(OUT))
    a = ap.parse_args()
    tasks = build_tasks(a.reps, a.quick)
    rows: list[dict] = []
    t0 = time.time()
    if a.jobs > 1:
        with ProcessPoolExecutor(max_workers=a.jobs) as ex:
            for k, res in enumerate(ex.map(run_task, tasks), 1):
                rows.extend(res)
                print(f"[{k}/{len(tasks)}] done ({time.time()-t0:.0f}s)", flush=True)
    else:
        for k, t in enumerate(tasks, 1):
            rows.extend(run_task(t))
            print(f"[{k}/{len(tasks)}] {t['scenario']} n={t['n']} rep={t['rep']} "
                  f"({time.time()-t0:.0f}s)", flush=True)
    out = pd.DataFrame(rows).sort_values(
        ["scenario", "population", "n_patients", "interval_mult", "rep",
         "estimator", "from_state", "to_state"]).reset_index(drop=True)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    out.drop(columns=["fit_seconds", "gen_seconds"]).to_csv(a.out, index=False)
    out.groupby(["estimator"])[["fit_seconds"]].describe().to_csv(
        Path(a.out).with_name("benchmark_timings.csv"))
    print(f"wrote {a.out} ({len(out)} rows) in {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
