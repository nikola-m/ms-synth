"""
ms_synth.summaries
==================
Single source of truth for the cohort-level quantities reported in the
manuscript (tables *and* figures use these functions).

Definitions
-----------
* Follow-up: time of the last visit (years from onset).
* Time to SPMS: first visit whose latent state lies in the SPMS phase;
  patients without such a visit are censored at their last visit.
* Time to EDSS >= k: first visit with true (noise-free) EDSS >= k; censored
  at the last visit otherwise.
* EDSS slope: ordinary least-squares slope of observed EDSS on time.
  A slope is counted as positive only if it exceeds ``SLOPE_TOL`` so that
  numerically flat trajectories (|slope| ~ 1e-16, whose sign depends on the
  BLAS build) are classified identically on every platform.
* Annualised relapse rate (ARR): pooled person-time, i.e. total observed
  relapse events / total observed follow-up in the stratum. Phase-specific
  person-time is partitioned at the first SPMS visit.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

SLOPE_TOL = 1e-9          # yr^-1; see module docstring
FAST_PROGRESSION = 0.5    # EDSS points / yr
EARLY_WINDOW_YR = 2.0     # window for "early relapses"


def _ols_slope(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 2 or x.std() == 0:
        return 0.0
    return float(np.polyfit(x, y, 1)[0])


def patient_table(df: pd.DataFrame) -> pd.DataFrame:
    """One row per patient with covariates, follow-up and event times."""
    d = df.sort_values(["patient_id", "disease_duration_yr"])
    rows = []
    for pid, g in d.groupby("patient_id", sort=False):
        t = g["disease_duration_yr"].to_numpy(dtype=float)
        et = g["EDSS_true"].to_numpy(dtype=float)
        spms = t[(g["phase"] == "SPMS").to_numpy()]
        rec = {
            "patient_id": pid,
            "age_at_onset": float(g["age_at_onset"].iloc[0]),
            "sex": g["sex"].iloc[0],
            "dmt": g["DMT_status"].iloc[0],
            "cls": g["latent_class"].iloc[0],
            "fu": float(t.max()),
            "n_visits": int(len(g)),
            "edss_base": float(g["EDSS"].iloc[0]),
            "edss_last": float(g["EDSS"].iloc[-1]),
            "n_relapse": int(g["relapse"].sum()),
            "t_spms": float(spms.min()) if len(spms) else np.nan,
            "slope": _ols_slope(t, g["EDSS"].to_numpy(dtype=float)),
            "early_relapse_cnt": int(g.loc[g["disease_duration_yr"] <= EARLY_WINDOW_YR,
                                           "relapse"].sum()),
        }
        for k in (3, 6, 8):
            hit = t[et >= k]
            rec[f"t_e{k}"] = float(hit.min()) if len(hit) else np.nan
        rows.append(rec)
    pat = pd.DataFrame(rows).set_index("patient_id")
    pat["ever_spms"] = pat["t_spms"].notna()
    pat["spms_event"] = pat["ever_spms"].astype(int)
    pat["spms_time"] = pat["t_spms"].where(pat["ever_spms"], pat["fu"])
    pat["male"] = (pat["sex"] == "M").astype(int)
    pat["early_relapse_bin"] = (pat["early_relapse_cnt"] >= 1).astype(int)
    pat["age_z"] = (pat["age_at_onset"] - pat["age_at_onset"].mean()) / pat["age_at_onset"].std()
    pat["edssbase_z"] = (pat["edss_base"] - pat["edss_base"].mean()) / pat["edss_base"].std()
    return pat


def milestone(pat: pd.DataFrame, col: str) -> tuple[pd.Series, pd.Series]:
    """(time, event) for a first-reaching milestone column such as 't_e6'."""
    ev = pat[col].notna().astype(int)
    return pat[col].where(ev == 1, pat["fu"]), ev


def visit_intervals_months(df: pd.DataFrame) -> np.ndarray:
    d = df.sort_values(["patient_id", "disease_duration_yr"])
    out = []
    for _, g in d.groupby("patient_id", sort=False):
        out.append(np.diff(g["disease_duration_yr"].to_numpy(dtype=float)) * 12.0)
    return np.concatenate(out) if out else np.array([])


def pooled_arr_by_phase(df: pd.DataFrame, pat: pd.DataFrame) -> tuple[float, float]:
    """Pooled person-time ARR in the RRMS and SPMS phases."""
    rr_t = rr_e = sp_t = sp_e = 0.0
    for pid, g in df.groupby("patient_id", sort=False):
        conv, last = pat.at[pid, "t_spms"], pat.at[pid, "fu"]
        rr_t += conv if not np.isnan(conv) else last
        rr_e += g.loc[g["phase"] == "RRMS", "relapse"].sum()
        if not np.isnan(conv):
            sp_t += max(last - conv, 0.0)
            sp_e += g.loc[g["phase"] == "SPMS", "relapse"].sum()
    return rr_e / max(rr_t, 1e-9), sp_e / max(sp_t, 1e-9)


def pooled_arr_by_group(pat: pd.DataFrame, col: str) -> pd.Series:
    """Pooled person-time ARR per level of a patient-level column."""
    g = pat.groupby(col)
    return g["n_relapse"].sum() / g["fu"].sum().clip(lower=1e-9)


def positive_slope_fraction(pat: pd.DataFrame) -> float:
    return float((pat["slope"] > SLOPE_TOL).mean())


def fast_progressor_fraction(pat: pd.DataFrame) -> float:
    return float((pat["slope"] > FAST_PROGRESSION).mean())
