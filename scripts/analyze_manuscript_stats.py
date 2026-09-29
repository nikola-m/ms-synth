#!/usr/bin/env python3
"""
Reproduce every statistic reported in the manuscript for the reference
cohort (N = 500, seed 42) and write results/results.json.

All definitions live in ms_synth.summaries (shared with make_figures.py).
Usage: python scripts/analyze_manuscript_stats.py [--out results/results.json]
"""
from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import numpy as np
from lifelines import CoxPHFitter, KaplanMeierFitter

from ms_synth import __version__
from ms_synth import summaries as S
from ms_synth.msm_fit import estimate_Q_crude, estimate_Q_from_paths, estimate_Q_panel
from ms_synth.synthetic_data import (
    build_Q_from_config, compare_empirical_to_Q, generate_synthetic_ms_dataset)
from ms_synth.utils import load_config

ROOT = Path(__file__).resolve().parent.parent
SEED, N = 42, 500


def km_summary(times, events, at=(5, 10, 15, 20)) -> dict:
    k = KaplanMeierFitter().fit(times, events)
    out = {str(t): round(100 * float(1 - k.survival_function_at_times(t).iloc[0]), 1) for t in at}
    med = k.median_survival_time_
    out["median"] = None if not np.isfinite(med) else round(float(med), 1)
    return out


def cox(pat, cols) -> dict:
    d = pat[cols + ["spms_time", "spms_event"]]
    s = CoxPHFitter().fit(d, "spms_time", "spms_event").summary
    return {c: {"HR": round(float(np.exp(s.loc[c, "coef"])), 2),
                "lo": round(float(np.exp(s.loc[c, "coef lower 95%"])), 2),
                "hi": round(float(np.exp(s.loc[c, "coef upper 95%"])), 2),
                "p": round(float(s.loc[c, "p"]), 4)} for c in cols}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(ROOT / "results" / "results.json"))
    a = ap.parse_args()
    logging.disable(logging.CRITICAL)

    cfg = load_config(str(ROOT / "configs" / "synthetic.yaml"))
    df, latent = generate_synthetic_ms_dataset(cfg=cfg, n_patients=N, seed=SEED, return_latent=True)
    pat = S.patient_table(df)
    iv = S.visit_intervals_months(df)
    R: dict = {"ms_synth_version": __version__, "seed": SEED, "n_patients": N}

    R.update(total_visits=int(len(df)),
             visits_per_patient_mean=round(float(pat.n_visits.mean()), 1),
             visits_per_patient_median=round(float(pat.n_visits.median()), 1),
             female_pct=round(100 * float((pat.sex == "F").mean()), 1),
             onset_age_mean=round(float(pat.age_at_onset.mean()), 1),
             onset_age_sd=round(float(pat.age_at_onset.std()), 1),
             fu_median=round(float(pat.fu.median()), 1), fu_mean=round(float(pat.fu.mean()), 1),
             fu_min=round(float(pat.fu.min()), 1), fu_max=round(float(pat.fu.max()), 1),
             interval_mean_m=round(float(iv.mean()), 2),
             interval_median_m=round(float(np.median(iv)), 2),
             interval_sd_m=round(float(iv.std(ddof=1)), 2),
             edss_base_mean=round(float(pat.edss_base.mean()), 2),
             edss_base_sd=round(float(pat.edss_base.std()), 2),
             edss_last_mean=round(float(pat.edss_last.mean()), 2),
             slope_mean=round(float(pat.slope.mean()), 3),
             slope_median=round(float(pat.slope.median()), 3),
             pos_slope_pct=round(100 * S.positive_slope_fraction(pat), 1),
             fast_prog_pct=round(100 * S.fast_progressor_fraction(pat), 1))

    R["spms_overall"] = km_summary(pat.spms_time, pat.spms_event)
    R["spms_overall_pct"] = round(100 * float(pat.ever_spms.mean()), 1)
    R["by_class"] = {}
    for c in ("stable", "moderate", "aggressive"):
        m = pat.cls == c
        km = km_summary(pat.loc[m, "spms_time"], pat.loc[m, "spms_event"], at=(10, 20))
        crude = {str(t): round(100 * float((pat.loc[m, "t_spms"] <= t).mean()), 1) for t in (10, 20)}
        R["by_class"][c] = {"n": int(m.sum()), "pct": round(100 * float(m.mean()), 1),
                            "spms_km": km, "spms_crude": crude,
                            "edss_last": round(float(pat.loc[m, "edss_last"].mean()), 2)}
    R["milestones"] = {}
    for name, col in (("EDSS3", "t_e3"), ("EDSS6", "t_e6"), ("EDSS8", "t_e8"), ("SPMS", "t_spms")):
        t, e = S.milestone(pat, col)
        R["milestones"][name] = km_summary(t, e, at=(10, 20))

    arr_rr, arr_sp = S.pooled_arr_by_phase(df, pat)
    R["arr_rrms"], R["arr_spms"] = round(float(arr_rr), 3), round(float(arr_sp), 3)
    arr = S.pooled_arr_by_group(pat, "dmt")
    R["dmt"] = {d: {"n": int((pat.dmt == d).sum()), "pct": round(100 * float((pat.dmt == d).mean()), 1),
                    "arr": round(float(arr[d]), 3),
                    "rel_reduction": None if d == "none" else round(100 * (1 - arr[d] / arr["none"]), 1)}
                for d in ("none", "moderate_dmt", "high_dmt")}

    # Generator recovery on the reference cohort: three estimators
    Q = build_Q_from_config(cfg)
    nz = [(i, j) for i in range(12) for j in range(12) if i != j and Q[i, j] > 0]
    comp = compare_empirical_to_Q(df, Q)
    R["Q_mean_relerr"] = round(float(comp.relative_error.mean()), 3)
    R["Q_median_relerr"] = round(float(comp.relative_error.median()), 3)
    R["Q_n_pairs"] = int(len(comp))
    est = {"oracle": estimate_Q_from_paths(latent, 12), "crude": estimate_Q_crude(df, 12),
           "panel": estimate_Q_panel(df, 12, allowed=nz)}
    R["Q_estimators"] = {}
    R["Q_table"] = [{"from_state": i, "to_state": j, "Q_true": float(Q[i, j]),
                     **{f"Q_{k}": round(float(v[i, j]), 4) for k, v in est.items()}} for i, j in nz]
    for k, v in est.items():
        rel = np.array([(v[i, j] - Q[i, j]) / Q[i, j] for i, j in nz])
        R["Q_estimators"][k] = {"mean_abs_relerr": round(float(np.abs(rel).mean()), 3),
                                "median_abs_relerr": round(float(np.median(np.abs(rel))), 3),
                                "mean_rel_bias": round(float(rel.mean()), 3)}

    R["cox"] = {"early_bin": cox(pat, ["early_relapse_bin"]),
                "early_cnt": cox(pat, ["early_relapse_cnt"]),
                "adjusted": cox(pat, ["age_z", "male", "edssbase_z"])}

    out = Path(a.out); out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(R, indent=2) + "\n")
    print(f"ms-synth {__version__}: N={N}, seed={SEED}, visits={R['total_visits']}, "
          f"SPMS KM20={R['spms_overall']['20']}%, ARR RRMS/SPMS={R['arr_rrms']}/{R['arr_spms']}")
    for k, v in R["Q_estimators"].items():
        print(f"  {k:6s} |rel err| mean {v['mean_abs_relerr']}  bias {v['mean_rel_bias']:+}")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
