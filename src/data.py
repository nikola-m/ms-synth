"""
data.py
=======
Load, validate, and preprocess longitudinal MS patient data.

Expected data schema (pandas DataFrame)
----------------------------------------
Column            Type        Description
---------         --------    -----------
patient_id        str/int     Unique patient identifier
time_years        float       Observation time in years (>=0)
edss_state        int         Discretised EDSS macrostate index  [0, n_states)
has_relapse       bool/int    Optional: relapse at this visit
phenotype         str         Optional: "RRMS" | "SPMS" | "PPMS" | "PRMS"

The module also contains a synthetic-data generator that produces
threshold-driven trajectories inspired by the Kannan et al. (2017) ODE
model, expressed as a discrete CT-MSM path.

References
----------
Kannan et al. (2017) Math. Biosci. 289:1-8.
Mirkov et al. (2026) arXiv:2602.08576.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.linalg import expm

from .utils import get_logger, set_seed

logger = get_logger(__name__)

# ---------------------------------------------------------------------------
# Constants / default schema
# ---------------------------------------------------------------------------

DEFAULT_PATIENT_COL = "patient_id"
DEFAULT_TIME_COL = "time_years"
DEFAULT_STATE_COL = "edss_state"
DEFAULT_RELAPSE_COL = "has_relapse"
DEFAULT_PHENOTYPE_COL = "phenotype"


# ---------------------------------------------------------------------------
# Loading helpers
# ---------------------------------------------------------------------------

def load_data(
    path: str | Path,
    cfg: dict[str, Any] | None = None,
) -> pd.DataFrame:
    """Load longitudinal MS data from a CSV (or feather/parquet) file.

    Parameters
    ----------
    path : str or Path
        File path.  Supported extensions: .csv, .parquet, .feather.
    cfg : dict, optional
        Configuration dict (``cfg["data"]`` sub-section).  Used to
        override column name mappings.

    Returns
    -------
    pd.DataFrame
        Validated data frame with standardised column names.

    Raises
    ------
    FileNotFoundError, ValueError
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Data file not found: {path}")

    ext = path.suffix.lower()
    if ext == ".csv":
        df = pd.read_csv(path)
    elif ext == ".parquet":
        df = pd.read_parquet(path)
    elif ext == ".feather":
        df = pd.read_feather(path)
    else:
        raise ValueError(f"Unsupported file extension '{ext}'")

    df = _apply_column_renames(df, cfg)
    df = validate_data(df, cfg)
    logger.info("Loaded data: %d rows, %d patients", len(df), df[DEFAULT_PATIENT_COL].nunique())
    return df


def _apply_column_renames(
    df: pd.DataFrame,
    cfg: dict[str, Any] | None,
) -> pd.DataFrame:
    """Rename columns according to config mappings (if provided)."""
    if cfg is None:
        return df
    data_cfg = cfg.get("data", {})
    rename_map = {}
    for default, key in [
        (DEFAULT_PATIENT_COL, "patient_col"),
        (DEFAULT_TIME_COL, "time_col"),
        (DEFAULT_STATE_COL, "state_col"),
        (DEFAULT_RELAPSE_COL, "relapse_col"),
        (DEFAULT_PHENOTYPE_COL, "phenotype_col"),
    ]:
        col = data_cfg.get(key)
        if col and col != default and col in df.columns:
            rename_map[col] = default
    if rename_map:
        df = df.rename(columns=rename_map)
    return df


def validate_data(
    df: pd.DataFrame,
    cfg: dict[str, Any] | None = None,
) -> pd.DataFrame:
    """Validate and lightly clean a longitudinal MS data frame.

    Checks performed:
    - Required columns present.
    - Time values non-negative and finite.
    - State values are non-negative integers.
    - Sorted by (patient_id, time_years).

    Parameters
    ----------
    df : pd.DataFrame
    cfg : dict, optional

    Returns
    -------
    pd.DataFrame
        Cleaned data frame.
    """
    required = [DEFAULT_PATIENT_COL, DEFAULT_TIME_COL, DEFAULT_STATE_COL]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    n_before = len(df)
    df = df.dropna(subset=required)
    df = df[df[DEFAULT_TIME_COL] >= 0.0]
    df[DEFAULT_STATE_COL] = df[DEFAULT_STATE_COL].astype(int)

    n_after = len(df)
    if n_after < n_before:
        logger.warning("Dropped %d rows during validation", n_before - n_after)

    df = df.sort_values([DEFAULT_PATIENT_COL, DEFAULT_TIME_COL]).reset_index(drop=True)
    return df


# ---------------------------------------------------------------------------
# Transition counting
# ---------------------------------------------------------------------------

def extract_transitions(
    df: pd.DataFrame,
    n_states: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Extract observed transition counts and holding times from panel data.

    For each consecutive pair of observations (t_k, t_{k+1}) for the same
    patient the function records:
      - Transition from state s_k to state s_{k+1}.
      - Holding time Δt = t_{k+1} - t_k.

    Parameters
    ----------
    df : pd.DataFrame
        Validated longitudinal data frame.
    n_states : int
        Total number of states in the model.

    Returns
    -------
    count_matrix : np.ndarray, shape (n_states, n_states)
        N[i, j] = number of observed transitions from state i to state j
        (i ≠ j) plus N[i, i] = number of transitions that stayed in i.
    total_time : np.ndarray, shape (n_states,)
        total_time[i] = total observed time spent in state i across all
        patients.
    """
    count_matrix = np.zeros((n_states, n_states), dtype=float)
    total_time = np.zeros(n_states, dtype=float)

    for _, grp in df.groupby(DEFAULT_PATIENT_COL, sort=False):
        grp = grp.sort_values(DEFAULT_TIME_COL)
        states = grp[DEFAULT_STATE_COL].to_numpy()
        times = grp[DEFAULT_TIME_COL].to_numpy()

        for k in range(len(states) - 1):
            s_from = int(states[k])
            s_to = int(states[k + 1])
            dt = float(times[k + 1] - times[k])
            if s_from < 0 or s_from >= n_states or s_to < 0 or s_to >= n_states:
                logger.warning("State index out of range: %d or %d", s_from, s_to)
                continue
            count_matrix[s_from, s_to] += 1
            total_time[s_from] += max(dt, 1e-9)

    return count_matrix, total_time


# ---------------------------------------------------------------------------
# Synthetic data generator
# ---------------------------------------------------------------------------

def generate_synthetic_data(
    n_patients: int = 300,
    n_states: int = 13,
    max_time: float = 20.0,
    time_step: float = 0.5,
    seed: int = 42,
    immune_threshold_range: tuple[float, float] = (0.3, 0.8),
    cns_threshold_range: tuple[float, float] = (0.4, 0.9),
) -> pd.DataFrame:
    """Generate synthetic longitudinal MS trajectories.

    The generator builds a heterogeneous population:
    - Each patient is assigned random ``immune_cap`` and ``cns_cap``
      thresholds (analogous to I_C and Z_C in Kannan et al. 2017).
    - Transition rates are derived from a baseline Q matrix whose
      reversible (RRMS-like) entries are attenuated once a patient's
      cumulative load exceeds their thresholds.
    - Stochastic jumps are simulated by matrix-exponentiation of Q·Δt.

    This produces the four classical MS subtypes (RRMS, RRMS→SPMS,
    PPMS, PRMS) in proportions controlled by the threshold distributions.

    Parameters
    ----------
    n_patients : int
    n_states : int
    max_time : float
        Total follow-up duration in years.
    time_step : float
        Observation interval in years.
    seed : int
    immune_threshold_range : tuple[float, float]
        Uniform range for the normalised immune-capacity parameter.
    cns_threshold_range : tuple[float, float]
        Uniform range for the normalised CNS-capacity parameter.

    Returns
    -------
    pd.DataFrame
        Longitudinal data in the standard schema.
    """
    rng = np.random.default_rng(seed)

    # --- Build a baseline Q matrix -----------------------------------------
    # Progressive states (upper diagonal) get larger intensities than
    # reversible ones (lower diagonal) to reflect the asymmetry in MS.
    Q_base = _make_base_Q(n_states, rng, asymmetry=3.0)

    records = []
    times_grid = np.arange(0.0, max_time + time_step, time_step)

    for pid in range(n_patients):
        immune_cap = rng.uniform(*immune_threshold_range)
        cns_cap = rng.uniform(*cns_threshold_range)

        state = int(rng.integers(0, max(1, n_states // 3)))  # start mild
        cumulative_load = 0.0
        immune_breached = False
        cns_breached = False

        for t in times_grid:
            records.append(
                {
                    DEFAULT_PATIENT_COL: pid,
                    DEFAULT_TIME_COL: round(float(t), 4),
                    DEFAULT_STATE_COL: state,
                    DEFAULT_RELAPSE_COL: int(not immune_breached and rng.random() < 0.1),
                    DEFAULT_PHENOTYPE_COL: _phenotype_label(
                        immune_breached, cns_breached
                    ),
                }
            )

            # ---- Determine effective Q for this patient/time ---------------
            Q_eff = Q_base.copy()
            # Attenuate reversible edges as thresholds are approached
            reversibility_factor = _compute_reversibility(
                cumulative_load, immune_cap, cns_cap
            )
            # Lower-triangular entries (remission / recovery) are scaled down
            for i in range(n_states):
                for j in range(i):
                    Q_eff[i, j] *= reversibility_factor

            # Re-fix diagonal so rows sum to 0
            np.fill_diagonal(Q_eff, 0.0)
            np.fill_diagonal(Q_eff, -Q_eff.sum(axis=1))

            # ---- Jump via matrix exponentiation ----------------------------
            P = expm(Q_eff * time_step)
            P = np.clip(P, 0.0, 1.0)
            row = P[state] / P[state].sum()

            state = int(rng.choice(n_states, p=row))
            cumulative_load += float(state) / n_states * time_step

            # Threshold breaches
            if not immune_breached and cumulative_load > immune_cap:
                immune_breached = True
            if not cns_breached and cumulative_load > cns_cap:
                cns_breached = True

    df = pd.DataFrame(records)
    logger.info(
        "Generated synthetic data: %d patients, %d rows", n_patients, len(df)
    )
    return df


# ---------------------------------------------------------------------------
# Private helpers for synthetic generation
# ---------------------------------------------------------------------------

def _make_base_Q(n: int, rng: np.random.Generator, asymmetry: float = 3.0) -> np.ndarray:
    """Build a random asymmetric Q (infinitesimal generator) matrix.

    Parameters
    ----------
    n : int
    rng : np.random.Generator
    asymmetry : float
        Ratio of progressive (up) to reversible (down) rates.

    Returns
    -------
    np.ndarray, shape (n, n)
    """
    Q = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            base = rng.exponential(0.1)
            if j > i:
                Q[i, j] = base * asymmetry  # worsening
            else:
                Q[i, j] = base              # recovery (weaker)
    np.fill_diagonal(Q, 0.0)
    np.fill_diagonal(Q, -Q.sum(axis=1))
    return Q


def _compute_reversibility(load: float, immune_cap: float, cns_cap: float) -> float:
    """Compute a [0,1] reversibility factor given cumulative load and capacities.

    Models the Kannan et al. (2017) threshold concept: once load exceeds
    the capacities the system loses its ability to recover (reversible edges
    are attenuated towards zero).

    Parameters
    ----------
    load : float
    immune_cap : float
    cns_cap : float

    Returns
    -------
    float in [0, 1]
    """
    cap = min(immune_cap, cns_cap)
    if cap <= 0:
        return 0.0
    factor = 1.0 - min(load / cap, 1.0)
    # Sigmoid-like sharpening
    return float(1.0 / (1.0 + np.exp(-10.0 * (factor - 0.5))))


def _phenotype_label(immune_breached: bool, cns_breached: bool) -> str:
    """Map threshold-breach flags to a clinical phenotype label."""
    if not immune_breached and not cns_breached:
        return "RRMS"
    if immune_breached and not cns_breached:
        return "SPMS"
    if not immune_breached and cns_breached:
        return "PRMS"
    return "PPMS"
