"""
synthetic_data.py
=================
Known-truth synthetic longitudinal MS cohort generator (ms-synth).

Architectural overview
----------------------
The generator uses a **continuous-time Markov jump process** (Gillespie
algorithm) to simulate patient trajectories through a configurable
state space.  Five design pillars:

1. **Clinically grounded Q matrix**
   Transition intensities (yr⁻¹) are NOT extracted from or fitted to any
   registry microdata. They are hand-specified and then calibrated so that
   AGGREGATE cohort outputs (SPMS conversion, relapse rates, time-to-
   milestone) approximate published natural-history targets. They reflect
   key clinical dynamics: high reversibility in early RRMS, progressive
   one-way intensification, PIRA pathways, and relapse cycling.

2. **Patient heterogeneity via latent classes**
   Each patient is assigned to one of three latent progression classes
   (stable / moderate / aggressive) that scale different sub-blocks of Q.
   This reproduces the well-known bimodal EDSS slope distribution.

3. **Covariate effects**
   Age-at-onset, sex, and DMT status multiplicatively modulate transition
   rates consistent with published hazard ratios.

4. **Realistic visit scheduling**
   Irregular observations are generated per patient: Gamma-distributed
   (shape=4) phase-dependent inter-visit intervals, extra visits post-
   relapse, and random censoring (dropout).

5. **EDSS measurement noise**
   True EDSS (from state centre) is perturbed with Gaussian noise and
   rounded to the clinical 0.5-step grid.

Data schema (output DataFrame columns)
---------------------------------------
patient_id          int       Unique patient identifier
visit_num           int       Sequential visit index (0-based)
visit_date          datetime  Absolute visit date (onset = random calendar date)
days_from_onset     float     Days since MS onset
EDSS                float     Observed EDSS (noisy, rounded to 0.5)
EDSS_true           float     True underlying EDSS (unrounded, for validation)
state               int       Macrostate index [0, n_states)
state_label         str       Human-readable macrostate label
relapse             int       1 if relapse occurred in window around visit
phase               str       "RRMS" | "SPMS" | "ABS"
age_at_visit        float     Patient age at this visit (years)
sex                 str       "F" | "M"
age_at_onset        float     Age at MS onset (constant per patient)
disease_duration_yr float     Years from onset at this visit
DMT_status          str       "none" | "moderate_dmt" | "high_dmt"
latent_class        str       "stable" | "moderate" | "aggressive"

Notes on public datasets
------------------------
Fully open, large (≥1000 patients) longitudinal EDSS+relapse datasets are
rare due to data-sharing restrictions.  Known candidates:

  • BRAINTEASER MS clinical dataset (Zenodo 10.5281/zenodo.7554921)
    – 1792 patients, EDSS + relapse visits, anonymised Italian/Spanish MS
    registry; however missing DMT granularity and limited relapses detail.
  • MSBase subset releases on Zenodo (doi:10.5281/zenodo.8373740)
    – aggregated statistics only; individual-level not public.
  • NARCOMS registry (Patient-determined EDSS, US): summary only.

Given these limitations this module generates fully synthetic data.
The ``load_or_generate`` function checks for a user-supplied real-data path
before falling back to synthesis.

Calibration targets (aggregate outputs are tuned toward these; the Q
entries are not fitted to their microdata)
-------------------------------------------------------------------
Weinshenker BG et al. (1989)  Brain 112(1):133-146.       [natural history]
Scalfari A et al. (2010)      Brain 133(7):1914-1929.     [London Ontario]
Confavreux C, Vukusic S (2006) Brain 129(3):606-616.     [Lyon]
Koch M et al. (2010)          JNNP 81:1039-1043.          [BCMS, SPMS]
Leray E et al. (2010)         Brain 133:1900-1913.        [Rennes, 2-stage]
Signori A et al. (2023)       JNNP 94:23-30.              [Big MS Data]
Polman CH et al. (2006)       NEJM 354:899-910.           [AFFIRM/natalizumab]
Kappos L et al. (2018)        Lancet 391:1263-1273.       [EXPAND/siponimod]

"""

from __future__ import annotations

import datetime
import logging
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.linalg import expm

from .utils import get_logger, set_seed, load_config, ensure_dir

logger = get_logger(__name__)


# ---------------------------------------------------------------------------
# State metadata helpers
# ---------------------------------------------------------------------------

# Default 12-state definitions if config not provided
_DEFAULT_STATE_DEFS: dict[int, dict] = {
    0:  {"label": "mild_RRMS_inactive",   "edss_centre": 0.75, "relapse": False, "phase": "RRMS"},
    1:  {"label": "mild_RRMS_active",     "edss_centre": 0.75, "relapse": True,  "phase": "RRMS"},
    2:  {"label": "mod_RRMS_inactive",    "edss_centre": 2.75, "relapse": False, "phase": "RRMS"},
    3:  {"label": "mod_RRMS_active",      "edss_centre": 2.75, "relapse": True,  "phase": "RRMS"},
    4:  {"label": "adv_RRMS_inactive",    "edss_centre": 4.75, "relapse": False, "phase": "RRMS"},
    5:  {"label": "adv_RRMS_active",      "edss_centre": 4.75, "relapse": True,  "phase": "RRMS"},
    6:  {"label": "early_prog_relapsing", "edss_centre": 6.50, "relapse": True,  "phase": "SPMS"},
    7:  {"label": "early_prog_inactive",  "edss_centre": 6.50, "relapse": False, "phase": "SPMS"},
    8:  {"label": "late_prog_relapsing",  "edss_centre": 8.00, "relapse": True,  "phase": "SPMS"},
    9:  {"label": "late_prog_inactive",   "edss_centre": 8.00, "relapse": False, "phase": "SPMS"},
    10: {"label": "severe_SPMS",          "edss_centre": 9.25, "relapse": False, "phase": "SPMS"},
    11: {"label": "absorbed",             "edss_centre": 9.50, "relapse": False, "phase": "ABS"},
}

# Default sparse Q entries (from_state, to_state, rate_per_year)
_DEFAULT_Q_ENTRIES: list[tuple[int, int, float]] = [
    (0, 1, 0.60), (1, 0, 0.85),
    (0, 2, 0.08), (1, 2, 0.05), (1, 3, 0.10),
    (2, 3, 0.45), (3, 2, 0.60),
    (2, 0, 0.04), (3, 1, 0.06),
    (2, 4, 0.07), (3, 4, 0.05), (3, 5, 0.08),
    (4, 5, 0.30), (5, 4, 0.40),
    (4, 7, 0.06), (5, 6, 0.08), (5, 7, 0.04),
    (6, 7, 0.50), (7, 6, 0.15),
    (6, 8, 0.10), (7, 8, 0.12),
    (8, 9, 0.35), (9, 8, 0.08),
    (8, 10, 0.08), (9, 10, 0.10),
    (10, 11, 0.05),
]


# ---------------------------------------------------------------------------
# Q matrix construction
# ---------------------------------------------------------------------------


def build_Q_from_config(cfg: dict[str, Any]) -> np.ndarray:
    """Assemble the base infinitesimal generator Q from sparse config entries.

    Parameters
    ----------
    cfg : dict
        Full configuration dictionary.  Reads ``cfg["states"]["n_states"]``
        and ``cfg["Q_base"]`` (list of [from, to, rate] triples).

    Returns
    -------
    Q : np.ndarray, shape (n, n)
        Valid generator matrix (rows sum to ≤ 0, off-diag ≥ 0).
    """
    n = int(cfg["states"]["n_states"])
    Q = np.zeros((n, n), dtype=float)

    entries = cfg.get("Q_base", [])
    if not entries:
        logger.warning("No Q_base in config; using built-in defaults.")
        entries_parsed = _DEFAULT_Q_ENTRIES
    else:
        entries_parsed = [(int(e[0]), int(e[1]), float(e[2])) for e in entries]

    for i, j, rate in entries_parsed:
        if 0 <= i < n and 0 <= j < n and i != j:
            Q[i, j] = rate
        else:
            logger.warning("Q entry (%d, %d, %.3f) out of range; skipped.", i, j, rate)

    # Fix diagonal
    np.fill_diagonal(Q, 0.0)
    np.fill_diagonal(Q, -Q.sum(axis=1))
    logger.debug("Built base Q: max off-diag = %.4f", float(np.max(Q - np.diag(np.diag(Q)))))
    return Q


def scale_Q_for_patient(
    Q_base: np.ndarray,
    latent_class: str,
    age_at_onset: float,
    sex: str,
    dmt_status: str,
    cfg: dict[str, Any],
    mean_onset_age: float = 32.0,
) -> np.ndarray:
    """Return a patient-specific Q matrix scaled by class, age, sex, and DMT.

    The scaling is applied multiplicatively to sub-blocks of Q:
    - **Progression edges** (i → j where j > i+1 or phase changes):
      scaled by ``progression_scale × age_factor × sex_factor × dmt_progression``.
    - **Relapse edges** (transitions to/from relapse-active states within band):
      scaled by ``relapse_scale × dmt_relapse``.
    - **Recovery edges** (j < i, returning towards remission):
      scaled by ``recovery_scale``.

    Parameters
    ----------
    Q_base : np.ndarray, shape (n, n)
    latent_class : str   "stable" | "moderate" | "aggressive"
    age_at_onset : float
    sex : str            "F" | "M"
    dmt_status : str     "none" | "moderate_dmt" | "high_dmt"
    cfg : dict
    mean_onset_age : float

    Returns
    -------
    Q_pat : np.ndarray, shape (n, n)
    """
    n = Q_base.shape[0]

    # ── Class factors ──────────────────────────────────────────────────────
    class_factors = cfg.get("class_Q_factors", {}).get(latent_class, {})
    prog_scale = float(class_factors.get("progression_scale", 1.0))
    relapse_scale = float(class_factors.get("relapse_scale", 1.0))
    recovery_scale = float(class_factors.get("recovery_scale", 1.0))

    # ── Age factor (per decade above mean) ────────────────────────────────
    age_eff = cfg.get("age_sex_effects", {})
    age_mult_per_decade = float(age_eff.get("age_progression_multiplier_per_decade", 1.25))
    decades_above = (age_at_onset - mean_onset_age) / 10.0
    age_factor = age_mult_per_decade ** decades_above  # >1 for older onset

    # ── Sex factor ─────────────────────────────────────────────────────────
    sex_mult = float(age_eff.get("male_progression_multiplier", 1.15))
    sex_factor = sex_mult if sex == "M" else 1.0

    # ── DMT factors ────────────────────────────────────────────────────────
    dmt_effs = age_eff.get("dmt_effects", {}).get(dmt_status, {})
    dmt_relapse_mult = float(dmt_effs.get("relapse_mult", 1.0))
    dmt_prog_mult = float(dmt_effs.get("progression_mult", 1.0))

    # ── State metadata for edge classification ─────────────────────────────
    state_defs = _parse_state_defs(cfg)

    Q_pat = Q_base.copy()

    for i in range(n):
        for j in range(n):
            if i == j or Q_base[i, j] == 0:
                continue

            phase_i = state_defs[i]["phase"]
            phase_j = state_defs[j]["phase"]
            edss_i = state_defs[i]["edss_centre"]
            edss_j = state_defs[j]["edss_centre"]
            relapse_i = state_defs[i]["relapse"]
            relapse_j = state_defs[j]["relapse"]

            # Classify edge
            is_phase_progression = (phase_j in {"SPMS", "ABS"} and phase_i == "RRMS")
            is_progression = (edss_j > edss_i + 0.5) or is_phase_progression
            is_recovery = (edss_j < edss_i - 0.5) and (phase_j != "ABS")
            is_relapse_edge = (relapse_j and not relapse_i) or (relapse_i and not relapse_j)

            if is_progression:
                Q_pat[i, j] *= prog_scale * age_factor * sex_factor * dmt_prog_mult
            elif is_recovery:
                Q_pat[i, j] *= recovery_scale
            elif is_relapse_edge:
                Q_pat[i, j] *= relapse_scale * dmt_relapse_mult

    # Re-fix diagonal
    np.fill_diagonal(Q_pat, 0.0)
    np.fill_diagonal(Q_pat, -Q_pat.sum(axis=1))
    return Q_pat


def _parse_state_defs(cfg: dict[str, Any]) -> dict[int, dict]:
    """Extract state definitions dict from config, falling back to defaults."""
    raw = cfg.get("states", {}).get("definitions", {})
    if not raw:
        return _DEFAULT_STATE_DEFS.copy()

    out = {}
    for k_raw, v in raw.items():
        k = int(k_raw)
        edss_low = float(v.get("edss_low", 0.0))
        edss_high = float(v.get("edss_high", 0.0))
        out[k] = {
            "label": str(v.get("label", str(k))),
            "edss_centre": (edss_low + edss_high) / 2.0,
            "relapse": bool(v.get("relapse", False)),
            "phase": str(v.get("phase", "RRMS")),
        }
    # Fill any missing states from defaults
    n = cfg["states"]["n_states"]
    for k in range(n):
        if k not in out:
            out[k] = _DEFAULT_STATE_DEFS.get(k, {
                "label": str(k), "edss_centre": 0.0, "relapse": False, "phase": "RRMS"
            })
    return out


# ---------------------------------------------------------------------------
# Gillespie / CTMC trajectory simulator
# ---------------------------------------------------------------------------


def simulate_patient_trajectory(
    Q_pat: np.ndarray,
    state_defs: dict[int, dict],
    rng: np.random.Generator,
    initial_state: int,
    max_time_years: float,
    min_time_years: float,
    dropout_rate: float,
    mean_interval_rrms: float,
    mean_interval_spms: float,
    post_relapse_interval: float,
    jitter_months: float,
    edss_noise_sigma: float,
    edss_round_to: float = 0.5,
    latent_sink: list | None = None,
) -> list[dict]:
    """Simulate a single patient's MS trajectory using the Gillespie algorithm.

    If ``latent_sink`` is a list, a dict ``{"events": [(t, state), ...],
    "follow_up": float}`` describing the exact latent path is appended to
    it. This consumes no random numbers, so enabling it never changes the
    simulated cohort.

    The Gillespie algorithm generates *exact* realisations of the CTMC:
    at each state the sojourn time is exponentially distributed with rate
    |Q[state, state]| = Σ_{j≠s} Q[s, j], and the next state is chosen
    proportionally to the individual rates.

    After generating the continuous event times, the trajectory is
    *observed* at irregular clinic visits (non-informative sampling
    conditional on the disease state) to produce the panel dataset.

    Parameters
    ----------
    Q_pat : np.ndarray, shape (n, n)
        Patient-specific generator matrix.
    state_defs : dict[int, dict]
        State metadata (label, edss_centre, relapse, phase).
    rng : np.random.Generator
    initial_state : int
    max_time_years : float
    min_time_years : float
    dropout_rate : float   Per-year probability of censoring.
    mean_interval_rrms : float   Mean visit interval (months) in RRMS phase.
    mean_interval_spms : float   Mean visit interval (months) in SPMS phase.
    post_relapse_interval : float  Extra visit (months) after a relapse.
    jitter_months : float   ± uniform jitter on visit intervals.
    edss_noise_sigma : float
    edss_round_to : float

    Returns
    -------
    list[dict]
        One dict per visit with all schema fields.
    """
    n = Q_pat.shape[0]

    # ── 1. Gillespie simulation ───────────────────────────────────────────
    events: list[tuple[float, int]] = []  # (time, state) at each jump
    current_state = initial_state
    t = 0.0
    events.append((0.0, current_state))

    while t < max_time_years:
        rate_out = max(-float(Q_pat[current_state, current_state]), 1e-12)
        sojourn = rng.exponential(1.0 / rate_out)
        t_jump = t + sojourn

        if t_jump >= max_time_years:
            break

        # Choose destination
        row = Q_pat[current_state].copy()
        row[current_state] = 0.0
        row = np.clip(row, 0.0, None)
        total = row.sum()
        if total <= 0:
            break  # absorbing
        probs = row / total
        next_state = int(rng.choice(n, p=probs))
        current_state = next_state
        t = t_jump
        events.append((t, current_state))

    if not events:
        events.append((0.0, initial_state))

    # ── 2. Observation schedule (irregular visits) ────────────────────────
    # Sample follow-up duration (subject to dropout)
    min_fu = min(min_time_years, max_time_years)
    follow_up = _sample_followup(
        rng, min_fu, max_time_years, dropout_rate
    )

    visit_times = _generate_visit_schedule(
        rng=rng,
        events=events,
        follow_up=follow_up,
        state_defs=state_defs,
        mean_interval_rrms=mean_interval_rrms,
        mean_interval_spms=mean_interval_spms,
        post_relapse_interval=post_relapse_interval,
        jitter_months=jitter_months,
    )

    if not visit_times:
        visit_times = [0.0]

    if latent_sink is not None:
        latent_sink.append({"events": list(events), "follow_up": float(follow_up)})

    # ── 3. Build visit records ────────────────────────────────────────────
    records = []
    for v_num, t_visit in enumerate(visit_times):
        # Interpolate true state at t_visit from Gillespie events
        true_state = _state_at_time(events, t_visit)
        meta = state_defs[true_state]

        # True EDSS: centre ± within-band noise
        edss_centre = float(meta["edss_centre"])
        edss_true = float(np.clip(
            edss_centre + rng.normal(0.0, 0.15),
            0.0, 9.5
        ))

        # Observed EDSS: add measurement noise, snap to grid
        edss_obs_raw = edss_true + rng.normal(0.0, edss_noise_sigma)
        edss_obs = _snap_to_grid(float(np.clip(edss_obs_raw, 0.0, 9.5)), edss_round_to)

        # Relapse: window ±0.25 years around any event that involved active state
        relapse_obs = _relapse_in_window(events, t_visit, state_defs, window=0.25)

        records.append({
            "visit_num": v_num,
            "time_years": round(t_visit, 6),
            "state": true_state,
            "state_label": meta["label"],
            "phase": meta["phase"],
            "EDSS_true": round(edss_true, 2),
            "EDSS": edss_obs,
            "relapse": int(relapse_obs),
        })

    return records


# ---------------------------------------------------------------------------
# Visit scheduling helpers
# ---------------------------------------------------------------------------


def _sample_followup(
    rng: np.random.Generator,
    min_fu: float,
    max_fu: float,
    dropout_rate: float,
) -> float:
    """Sample follow-up duration with random censoring.

    Censoring modelled as Exponential(dropout_rate); follow-up is
    min(censoring_time, max_fu), floored at min_fu.
    """
    if dropout_rate <= 0:
        return max_fu
    censor_time = rng.exponential(1.0 / dropout_rate)
    return float(np.clip(censor_time, min_fu, max_fu))


def _generate_visit_schedule(
    rng: np.random.Generator,
    events: list[tuple[float, int]],
    follow_up: float,
    state_defs: dict[int, dict],
    mean_interval_rrms: float,
    mean_interval_spms: float,
    post_relapse_interval: float,
    jitter_months: float,
) -> list[float]:
    """Generate irregular clinic visit times for one patient.

    Visits use Gamma-distributed (shape=4) phase-dependent inter-visit
    intervals (not a Poisson/exponential process).
    A post-relapse extra visit is inserted after any relapse event.

    Parameters
    ----------
    events : list[(time, state)]
        Gillespie event list (continuous-time jumps).
    follow_up : float   Total follow-up in years.

    Returns
    -------
    list[float]   Sorted visit times in years.
    """
    visits: list[float] = [0.0]  # always have baseline visit
    t = 0.0

    # Determine phase at each time for interval choice
    event_idx = 0

    while True:
        # Current state
        cur_state = _state_at_time(events, t)
        cur_phase = state_defs[cur_state]["phase"]

        if cur_phase == "RRMS":
            mean_months = mean_interval_rrms
        else:
            mean_months = mean_interval_spms

        # Sample interval (Gamma to avoid too-short intervals)
        interval_yr = _sample_interval(rng, mean_months, jitter_months)
        t_next = t + interval_yr

        if t_next > follow_up:
            break

        # Check if a relapse event happened in (t, t_next) → add extra visit
        # (outcome-dependent: the extra visit depends on the latent jump
        # time; set post_relapse_visit_months <= 0 to disable)
        for ev_t, ev_s in (events if post_relapse_interval > 0 else ()):
            if t < ev_t < t_next and state_defs[ev_s]["relapse"]:
                extra = ev_t + post_relapse_interval / 12.0
                if extra < follow_up and extra not in visits:
                    visits.append(round(extra, 6))
                break

        visits.append(round(t_next, 6))
        t = t_next

    visits = sorted(set(visits))
    return visits


def _sample_interval(
    rng: np.random.Generator,
    mean_months: float,
    jitter_months: float,
) -> float:
    """Sample one inter-visit interval in years."""
    # Gamma with mean=mean_months, shape=4 (moderate variance)
    shape = 4.0
    scale = mean_months / shape
    raw = rng.gamma(shape=shape, scale=scale)
    raw = max(raw + rng.uniform(-jitter_months, jitter_months), 1.0)
    return raw / 12.0  # convert months → years


def _state_at_time(events: list[tuple[float, int]], t: float) -> int:
    """Return the state occupied at time t from a Gillespie event list."""
    state = events[0][1]
    for ev_t, ev_s in events:
        if ev_t <= t:
            state = ev_s
        else:
            break
    return state


def _relapse_in_window(
    events: list[tuple[float, int]],
    t_visit: float,
    state_defs: dict[int, dict],
    window: float = 0.25,
) -> bool:
    """Check if any relapse event occurred within ±window years of t_visit."""
    for ev_t, ev_s in events:
        if abs(ev_t - t_visit) <= window and state_defs[ev_s]["relapse"]:
            return True
    return False


def _snap_to_grid(value: float, step: float) -> float:
    """Round a float to the nearest multiple of ``step``."""
    if step <= 0:
        return value
    return round(round(value / step) * step, 1)


# ---------------------------------------------------------------------------
# Patient demographics sampler
# ---------------------------------------------------------------------------


def sample_demographics(
    rng: np.random.Generator,
    cfg: dict[str, Any],
) -> tuple[str, float, str, str]:
    """Sample one patient's demographics from config distributions.

    Parameters
    ----------
    rng : np.random.Generator
    cfg : dict

    Returns
    -------
    (latent_class, age_at_onset, sex, dmt_status) : tuple
    """
    demo = cfg.get("demographics", {})
    cohort = cfg.get("cohort", {})

    # Latent class
    class_fracs = cohort.get("class_fractions", {"stable": 0.3, "moderate": 0.45, "aggressive": 0.25})
    classes = list(class_fracs.keys())
    class_probs = np.array([class_fracs[c] for c in classes], dtype=float)
    class_probs /= class_probs.sum()
    latent_class = str(classes[rng.choice(len(classes), p=class_probs)])

    # Age at onset
    age_cfg = demo.get("onset_age", {})
    dist = age_cfg.get("distribution", "normal")
    age: float
    if dist == "normal":
        age = float(rng.normal(
            float(age_cfg.get("mean", 32.0)),
            float(age_cfg.get("std", 9.0)),
        ))
        age = float(np.clip(age, float(age_cfg.get("min", 16.0)), float(age_cfg.get("max", 60.0))))
    else:
        age = float(age_cfg.get("mean", 32.0))

    # Sex
    female_prob = float(demo.get("sex_female_prob", 0.68))
    sex = "F" if rng.random() < female_prob else "M"

    # DMT status
    dmt_probs_raw = demo.get("dmt_status_probs", {"none": 0.2, "moderate_dmt": 0.35, "high_dmt": 0.45})
    dmt_keys = list(dmt_probs_raw.keys())
    dmt_probs = np.array([dmt_probs_raw[k] for k in dmt_keys], dtype=float)
    dmt_probs /= dmt_probs.sum()
    dmt_status = str(dmt_keys[rng.choice(len(dmt_keys), p=dmt_probs)])

    return latent_class, age, sex, dmt_status


# ---------------------------------------------------------------------------
# Main dataset generator
# ---------------------------------------------------------------------------


def generate_synthetic_ms_dataset(
    cfg: dict[str, Any] | None = None,
    cfg_path: str | Path | None = None,
    n_patients: int | None = None,
    seed: int | None = None,
    return_latent: bool = False,
) -> pd.DataFrame | tuple[pd.DataFrame, dict]:
    """Generate a synthetic longitudinal MS dataset.

    With ``return_latent=True`` the function returns ``(df, latent)`` where
    ``latent[patient_id] = {"events", "follow_up", "Q_patient"}`` holds the
    exact continuous-time path and the patient-specific generator (the
    known truth). The observed ``df`` is identical in both modes.

    Combines all sub-components (Q construction, per-patient scaling,
    Gillespie simulation, visit scheduling, EDSS noise) into a single
    cohort DataFrame.

    Parameters
    ----------
    cfg : dict, optional
        Pre-loaded configuration.  If None, ``cfg_path`` must be provided.
    cfg_path : str or Path, optional
        Path to YAML / JSON config file.
    n_patients : int, optional
        Override ``cfg["cohort"]["n_patients"]``.
    seed : int, optional
        Override ``cfg["reproducibility"]["seed"]``.

    Returns
    -------
    pd.DataFrame
        Full longitudinal dataset with all schema columns.

    Raises
    ------
    ValueError
        If neither cfg nor cfg_path is provided.
    """
    if cfg is None:
        if cfg_path is None:
            raise ValueError("Provide either cfg dict or cfg_path.")
        cfg = load_config(cfg_path)

    # ── Setup ──────────────────────────────────────────────────────────────
    _seed = seed if seed is not None else int(cfg.get("reproducibility", {}).get("seed", 42))
    set_seed(_seed)
    rng = np.random.default_rng(_seed)

    _n = n_patients if n_patients is not None else int(cfg["cohort"]["n_patients"])
    n_states = int(cfg["states"]["n_states"])
    state_defs = _parse_state_defs(cfg)

    onset_probs = cfg["states"].get("onset_state_probs")
    if onset_probs is None:
        onset_probs = [1.0 / n_states] * n_states
    onset_probs = np.array(onset_probs, dtype=float)
    onset_probs = onset_probs[:n_states]
    onset_probs /= onset_probs.sum()

    # Visit config
    vis_cfg = cfg.get("visits", {})
    max_fu = float(vis_cfg.get("max_followup_years", 20.0))
    min_fu = float(vis_cfg.get("min_followup_years", 3.0))
    dropout_rate = float(vis_cfg.get("dropout_rate_per_year", 0.03))
    mean_rrms = float(vis_cfg.get("mean_interval_rrms_months", 7.0))
    mean_spms = float(vis_cfg.get("mean_interval_spms_months", 9.0))
    post_relapse = float(vis_cfg.get("post_relapse_visit_months", 1.5))
    jitter = float(vis_cfg.get("jitter_months", 2.0))

    edss_sigma = float(cfg.get("edss_noise", {}).get("sigma", 0.30))
    edss_round = float(cfg.get("edss_noise", {}).get("round_to", 0.5))

    # Build base Q
    Q_base = build_Q_from_config(cfg)

    logger.info(
        "Starting synthetic dataset generation: n_patients=%d, n_states=%d, seed=%d",
        _n, n_states, _seed,
    )

    # ── Per-patient simulation ─────────────────────────────────────────────
    all_records: list[dict] = []
    latent: dict = {}
    mean_onset_age = float(cfg.get("demographics", {}).get("onset_age", {}).get("mean", 32.0))

    # Random calendar onset dates spanning a ~20-year accrual window
    # (simulates a real registry where patients enrol at different times)
    ref_date = datetime.date(2000, 1, 1)
    accrual_days = 365 * 20

    stats = {"n_spms": 0, "total_visits": 0}

    for pid in range(_n):
        if pid % 500 == 0 and pid > 0:
            logger.info("  Generated %d / %d patients…", pid, _n)

        # Demographics
        latent_class, age_onset, sex, dmt_status = sample_demographics(rng, cfg)

        # Patient-specific Q
        Q_pat = scale_Q_for_patient(
            Q_base, latent_class, age_onset, sex, dmt_status, cfg, mean_onset_age
        )

        # Initial state
        init_state = int(rng.choice(n_states, p=onset_probs))

        # Gillespie trajectory + visits
        _sink: list | None = [] if return_latent else None
        visit_records = simulate_patient_trajectory(
            Q_pat=Q_pat,
            state_defs=state_defs,
            rng=rng,
            initial_state=init_state,
            max_time_years=max_fu,
            min_time_years=min_fu,
            dropout_rate=dropout_rate,
            mean_interval_rrms=mean_rrms,
            mean_interval_spms=mean_spms,
            post_relapse_interval=post_relapse,
            jitter_months=jitter,
            edss_noise_sigma=edss_sigma,
            edss_round_to=edss_round,
            latent_sink=_sink,
        )
        if return_latent:
            latent[pid] = {**_sink[0], "Q_patient": Q_pat.copy()}

        # Calendar onset date (random within accrual window)
        onset_day_offset = int(rng.integers(0, accrual_days))
        onset_date = ref_date + datetime.timedelta(days=onset_day_offset)

        # Check if patient ever progressed to SPMS
        if any(r["phase"] == "SPMS" for r in visit_records):
            stats["n_spms"] += 1

        for rec in visit_records:
            t_yr = rec["time_years"]
            visit_date = onset_date + datetime.timedelta(days=int(t_yr * 365.25))
            all_records.append({
                "patient_id": pid,
                "visit_num": rec["visit_num"],
                "visit_date": visit_date,
                "days_from_onset": round(t_yr * 365.25, 1),
                "EDSS": rec["EDSS"],
                "EDSS_true": rec["EDSS_true"],
                "state": rec["state"],
                "state_label": rec["state_label"],
                "relapse": rec["relapse"],
                "phase": rec["phase"],
                "age_at_visit": round(age_onset + t_yr, 2),
                "sex": sex,
                "age_at_onset": round(age_onset, 2),
                "disease_duration_yr": round(t_yr, 4),
                "DMT_status": dmt_status,
                "latent_class": latent_class,
            })

        stats["total_visits"] += len(visit_records)

    # ── Assemble DataFrame ─────────────────────────────────────────────────
    df = pd.DataFrame(all_records)
    df["visit_date"] = pd.to_datetime(df["visit_date"])

    pct_spms = 100.0 * stats["n_spms"] / max(_n, 1)
    mean_visits = stats["total_visits"] / max(_n, 1)
    logger.info(
        "Generation complete: %d patients, %d visits, %.1f%% SPMS conversion, "
        "mean %.1f visits/patient",
        _n, len(df), pct_spms, mean_visits,
    )

    # ── Validation logging ─────────────────────────────────────────────────
    _log_validation_stats(df)

    if return_latent:
        return df, latent
    return df


# ---------------------------------------------------------------------------
# Validation helpers
# ---------------------------------------------------------------------------


def _log_validation_stats(df: pd.DataFrame) -> None:
    """Log key epidemiological statistics for quick plausibility check."""
    if df.empty:
        return
    n_pat = df["patient_id"].nunique()
    visits_per_pat = df.groupby("patient_id").size()
    spms_pats = df[df["phase"] == "SPMS"]["patient_id"].nunique()
    _g = df.groupby("patient_id")
    relapse_rate_yr = (
        _g["relapse"].sum() / _g["disease_duration_yr"].max().clip(lower=1e-6)
    ).mean()
    median_fu = df.groupby("patient_id")["disease_duration_yr"].max().median()
    class_dist = df.drop_duplicates("patient_id")["latent_class"].value_counts().to_dict()
    sex_dist = df.drop_duplicates("patient_id")["sex"].value_counts().to_dict()

    logger.info("── Validation statistics ──────────────────────────────────")
    logger.info("  Patients:              %d", n_pat)
    logger.info("  Total visits:          %d", len(df))
    logger.info("  Median visits/patient: %.1f", float(visits_per_pat.median()))
    logger.info("  Median follow-up (yr): %.1f", float(median_fu))
    logger.info("  SPMS conversion:       %d (%.1f%%)", spms_pats, 100.0 * spms_pats / n_pat)
    logger.info("  Mean relapse rate/yr:  %.3f", float(relapse_rate_yr))
    logger.info("  Latent class dist:     %s", class_dist)
    logger.info("  Sex distribution:      %s", sex_dist)
    logger.info("───────────────────────────────────────────────────────────")


def compute_empirical_transition_counts(
    df: pd.DataFrame,
    n_states: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute empirical inter-visit transition counts from the dataset.

    Thin wrapper around ``data.extract_transitions`` that maps the
    synthetic dataset's column names.

    Parameters
    ----------
    df : pd.DataFrame
    n_states : int

    Returns
    -------
    (count_matrix, total_time) : (np.ndarray shape (n,n), np.ndarray shape (n,))
    """
    _tmp = df.rename(columns={
        "patient_id": "patient_id",
        "disease_duration_yr": "time_years",
        "state": "edss_state",
    })
    from .data import extract_transitions
    return extract_transitions(_tmp, n_states)


def compare_empirical_to_Q(
    df: pd.DataFrame,
    Q_true: np.ndarray,
) -> pd.DataFrame:
    """Compare empirical transition rates to the generating Q matrix.

    Computes λ̂_{ij} = N_{ij}/T_i and returns a DataFrame of true vs
    estimated rates for each non-zero pair.

    Parameters
    ----------
    df : pd.DataFrame
    Q_true : np.ndarray

    Returns
    -------
    pd.DataFrame
        Columns: from_state, to_state, Q_true, Q_empirical, relative_error.
    """
    n = Q_true.shape[0]
    cnt, tot = compute_empirical_transition_counts(df, n)
    rows = []
    for i in range(n):
        for j in range(n):
            if i == j or Q_true[i, j] == 0:
                continue
            q_emp = float(cnt[i, j]) / max(float(tot[i]), 1e-9)
            rows.append({
                "from_state": i,
                "to_state": j,
                "Q_true": float(Q_true[i, j]),
                "Q_empirical": round(q_emp, 5),
                "relative_error": round(abs(q_emp - Q_true[i, j]) / max(Q_true[i, j], 1e-9), 4),
            })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Public entry point (compatible with data.py's load_or_generate pattern)
# ---------------------------------------------------------------------------


def load_or_generate(
    cfg: dict[str, Any],
    real_data_path: str | Path | None = None,
) -> pd.DataFrame:
    """Load real data (if available) or generate synthetic data.

    If ``real_data_path`` is provided and the file exists, load it;
    otherwise generate synthetically.  Ensures the output conforms
    to the standard pipeline schema (columns: patient_id, time_years,
    edss_state for compatibility with ``data.validate_data``).

    Parameters
    ----------
    cfg : dict
    real_data_path : str or Path, optional

    Returns
    -------
    pd.DataFrame
    """
    if real_data_path is not None:
        p = Path(real_data_path)
        if p.exists():
            logger.info("Loading real data from %s", p)
            from .data import load_data
            df = load_data(p, cfg=cfg)
            return df
        else:
            logger.warning("Real data path '%s' not found; generating synthetic data.", p)

    df = generate_synthetic_ms_dataset(cfg=cfg)

    # Add pipeline-compatible alias columns
    df["time_years"] = df["disease_duration_yr"]
    df["edss_state"] = df["state"]
    return df
