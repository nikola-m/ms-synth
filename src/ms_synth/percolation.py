"""
percolation.py
==============
Core percolation analysis on the directed state-transition graph of a
CT-MSM fitted to MS progression data.

Scientific background
---------------------
Mirkov et al. (2026) hypothesise that the RRMS→SPMS transition mirrors a
percolation phase transition in the brain's connectivity network: the system
maintains a Giant Strongly Connected Component (GSCC) while lesion/
transition-intensity thresholds are below a critical value θ_c, then
fragments abruptly above it.  Here we operationalise this on the *abstract*
state-transition graph whose edges are CT-MSM intensities.

As the threshold θ increases (weak edges pruned), the directed graph loses:
  1. The GSCC (largest strongly connected component) – represents the pool
     of macrostates that can still cycle between relapse and remission.
  2. The out-component of mild RRMS states – how far progressive worsening
     can spread.
  3. The in-component of progressive states – which early states can feed
     into progressive trajectories.

A critical threshold θ_c is identified as the inflection / steepest drop
of the normalised GSCC size curve S(θ) = |GSCC(θ)| / N.

Derived percolation statistics (after Lebkuecher et al. 2024 / Lall 2024):
  - Percolation integral  I_perc  = ∫ S(θ) dθ  (area under the GSCC curve).
    Higher I_perc → more robust network (slower GSCC collapse).
  - Susceptibility  χ(θ) = Var(component sizes) near θ_c.
  - Proximity to criticality  δ(θ) = θ_c - θ  (positive = sub-critical).
  - Time-varying δ(t):  for each patient visit, the patient-specific Q is
    used to compute δ, yielding a trajectory of "distance from the tipping
    point" over time.

References
----------
Kannan et al. (2017) Math. Biosci. 289:1-8.
Mirkov et al. (2026) arXiv:2602.08576.
Lebkuecher et al. (2024) Neuropsychology 38:42-57.
Lall (2024) PhD thesis, Montclair State University.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter

from .graph_utils import (
    build_theta_grid,
    get_gscc,
    get_in_component,
    get_out_component,
    threshold_graph,
    normalise_weights,
)
from .utils import get_logger

logger = get_logger(__name__)


# ---------------------------------------------------------------------------
# Data containers
# ---------------------------------------------------------------------------


@dataclass
class PercolationCurve:
    """Results of a single threshold sweep on a directed graph.

    Attributes
    ----------
    theta_grid : np.ndarray, shape (T,)
        Threshold values used.
    gscc_sizes : np.ndarray, shape (T,)
        Absolute GSCC node count at each θ.
    gscc_fractions : np.ndarray, shape (T,)
        GSCC size normalised by total node count.
    in_comp_sizes : np.ndarray, shape (T,)
        In-component size (relative to starting_states) at each θ.
    out_comp_sizes : np.ndarray, shape (T,)
        Out-component size at each θ.
    n_edges : np.ndarray, shape (T,)
        Number of surviving edges at each θ.
    susceptibility : np.ndarray, shape (T,)
        Variance of non-trivial SCC sizes at each θ.
    theta_c : float
        Estimated critical threshold (inflection point of gscc_fractions).
    theta_c_method : str
        Criticality detection method used.
    percolation_integral : float
        ∫ gscc_fractions(θ) dθ  (area under the GSCC curve).
    susceptibility_integral : float
        ∫ susceptibility(θ) dθ  (integrated susceptibility).
    n_states : int
        Total number of states (nodes).
    starting_states : list[int]
        Source states for in/out component computation.
    """
    theta_grid: np.ndarray
    gscc_sizes: np.ndarray
    gscc_fractions: np.ndarray
    in_comp_sizes: np.ndarray
    out_comp_sizes: np.ndarray
    n_edges: np.ndarray
    susceptibility: np.ndarray
    theta_c: float
    theta_c_method: str
    percolation_integral: float
    susceptibility_integral: float
    n_states: int
    starting_states: list[int] = field(default_factory=list)

    def delta(self, theta: float) -> float:
        """Proximity to criticality: δ = θ_c - θ  (positive = sub-critical).

        Parameters
        ----------
        theta : float
            Current threshold value.

        Returns
        -------
        float
        """
        return float(self.theta_c - theta)

    def to_dataframe(self) -> pd.DataFrame:
        """Export sweep results as a tidy DataFrame.

        Returns
        -------
        pd.DataFrame
            Columns: theta, gscc_size, gscc_fraction, in_comp_size,
            out_comp_size, n_edges, susceptibility, delta.
        """
        return pd.DataFrame(
            {
                "theta": self.theta_grid,
                "gscc_size": self.gscc_sizes,
                "gscc_fraction": self.gscc_fractions,
                "in_comp_size": self.in_comp_sizes,
                "out_comp_size": self.out_comp_sizes,
                "n_edges": self.n_edges,
                "susceptibility": self.susceptibility,
                "delta": self.theta_c - self.theta_grid,
            }
        )


# ---------------------------------------------------------------------------
# Main threshold sweep
# ---------------------------------------------------------------------------


def threshold_sweep(
    G,  # nx.DiGraph
    cfg: dict[str, Any] | None = None,
    starting_states: list[int] | None = None,
) -> PercolationCurve:
    """Sweep the threshold θ and compute directed percolation metrics.

    At each value of θ edges with weight < θ are removed; the resulting
    directed graph is analysed for GSCC size, in/out component sizes,
    and susceptibility (variance of SCC sizes).

    Parameters
    ----------
    G : nx.DiGraph
        Full directed graph (e.g., from ``graph_utils.build_graph``).
    cfg : dict, optional
        Configuration dict.  Reads ``cfg["percolation"]`` sub-section.
    starting_states : list[int], optional
        RRMS-like source nodes for in/out component computation.
        Falls back to ``cfg["states"]["starting_states"]`` or [0, 1, 2].

    Returns
    -------
    PercolationCurve
    """
    perc_cfg = (cfg or {}).get("percolation", {})
    n_theta = int(perc_cfg.get("n_theta", 100))
    method = perc_cfg.get("threshold_method", "absolute")
    theta_min = float(perc_cfg.get("theta_min", 0.0))
    theta_max = perc_cfg.get("theta_max", None)
    if theta_max is not None:
        theta_max = float(theta_max)
    criticality_method = perc_cfg.get("criticality_method", "steepest_drop")

    if starting_states is None:
        starting_states = (cfg or {}).get("states", {}).get("starting_states", [0, 1, 2])

    # Optionally normalise for "normalized_prob" mode
    G_work = G
    if method == "normalized_prob":
        G_work = normalise_weights(G, method="row")

    theta_grid = build_theta_grid(G_work, n_theta=n_theta, theta_min=theta_min, theta_max=theta_max)
    n_states = G_work.number_of_nodes()

    gscc_sizes = np.zeros(n_theta, dtype=float)
    in_comp_sizes = np.zeros(n_theta, dtype=float)
    out_comp_sizes = np.zeros(n_theta, dtype=float)
    n_edges_arr = np.zeros(n_theta, dtype=float)
    susceptibility = np.zeros(n_theta, dtype=float)

    logger.info(
        "Starting threshold sweep: n_theta=%d, method=%s, θ ∈ [%.4f, %.4f]",
        n_theta, method, theta_grid[0], theta_grid[-1],
    )

    for k, theta in enumerate(theta_grid):
        G_k = threshold_graph(G_work, theta=theta, method=method)
        n_edges_arr[k] = G_k.number_of_edges()

        # GSCC
        gscc = get_gscc(G_k)
        gscc_sizes[k] = len(gscc)

        # In / out components relative to starting states
        in_comp = get_in_component(G_k, starting_states)
        out_comp = get_out_component(G_k, starting_states)
        in_comp_sizes[k] = len(in_comp)
        out_comp_sizes[k] = len(out_comp)

        # Susceptibility: variance of SCC sizes (excluding trivial size-1 SCCs)
        sccs = list(nx.strongly_connected_components(G_k))
        sizes = [len(s) for s in sccs if len(s) > 1]
        susceptibility[k] = float(np.var(sizes)) if len(sizes) > 1 else 0.0

    import networkx as nx  # ensure import for inner use

    gscc_fractions = gscc_sizes / max(n_states, 1)

    # Smooth the GSCC curve for reliable criticality detection
    if n_theta >= 11:
        win = min(11, n_theta if n_theta % 2 == 1 else n_theta - 1)
        gscc_smooth = savgol_filter(gscc_fractions, window_length=win, polyorder=2)
    else:
        gscc_smooth = gscc_fractions.copy()

    # Detect θ_c
    theta_c = _detect_criticality(
        theta_grid, gscc_smooth, method=criticality_method
    )

    # Percolation integral (area under GSCC curve) using trapezoidal rule
    perc_integral = float(np.trapz(gscc_fractions, theta_grid))
    susc_integral = float(np.trapz(susceptibility, theta_grid))

    curve = PercolationCurve(
        theta_grid=theta_grid,
        gscc_sizes=gscc_sizes,
        gscc_fractions=gscc_fractions,
        in_comp_sizes=in_comp_sizes,
        out_comp_sizes=out_comp_sizes,
        n_edges=n_edges_arr,
        susceptibility=susceptibility,
        theta_c=theta_c,
        theta_c_method=criticality_method,
        percolation_integral=perc_integral,
        susceptibility_integral=susc_integral,
        n_states=n_states,
        starting_states=list(starting_states),
    )

    logger.info(
        "Threshold sweep complete: θ_c=%.4f, I_perc=%.4f, I_susc=%.4f",
        theta_c, perc_integral, susc_integral,
    )
    return curve


# Lazy networkx import at module body (avoid ordering issue with type hints)
import networkx as nx  # noqa: E402


# ---------------------------------------------------------------------------
# Criticality detection
# ---------------------------------------------------------------------------


def _detect_criticality(
    theta_grid: np.ndarray,
    gscc_fractions: np.ndarray,
    method: str = "steepest_drop",
) -> float:
    """Identify the critical threshold θ_c from the GSCC(θ) curve.

    Parameters
    ----------
    theta_grid : np.ndarray
    gscc_fractions : np.ndarray
        Normalised GSCC size at each θ (may be smoothed).
    method : str
        "steepest_drop"  – θ_c = argmax |dS/dθ|  (finite difference).
        "half_max"       – θ_c = θ where S(θ) first drops below 0.5.
        "sigmoid_fit"    – fit a decreasing sigmoid and return inflection.

    Returns
    -------
    theta_c : float
    """
    if method == "steepest_drop":
        return _criticality_steepest_drop(theta_grid, gscc_fractions)
    elif method == "half_max":
        return _criticality_half_max(theta_grid, gscc_fractions)
    elif method == "sigmoid_fit":
        return _criticality_sigmoid_fit(theta_grid, gscc_fractions)
    else:
        logger.warning("Unknown criticality method '%s'; using steepest_drop", method)
        return _criticality_steepest_drop(theta_grid, gscc_fractions)


def _criticality_steepest_drop(
    theta: np.ndarray, S: np.ndarray
) -> float:
    """θ_c = θ at the steepest descent of S(θ).

    Uses finite differences of the smoothed curve.
    """
    if len(S) < 2:
        return float(theta[0])
    dS = np.diff(S)
    # Most negative gradient = steepest drop
    idx = int(np.argmin(dS))
    # Interpolate between idx and idx+1
    theta_c = 0.5 * (float(theta[idx]) + float(theta[idx + 1]))
    return theta_c


def _criticality_half_max(
    theta: np.ndarray, S: np.ndarray
) -> float:
    """θ_c = first θ where S(θ) drops below 0.5·S_max."""
    S_max = float(np.max(S))
    half = 0.5 * S_max
    below = np.where(S < half)[0]
    if len(below) == 0:
        return float(theta[-1])
    return float(theta[below[0]])


def _criticality_sigmoid_fit(
    theta: np.ndarray, S: np.ndarray
) -> float:
    """Fit a decreasing logistic S(θ) = L / (1 + exp(k(θ - θ_c))) and return θ_c."""

    def sigmoid(x, L, k, x0):
        return L / (1.0 + np.exp(k * (x - x0)))

    try:
        # Initial guess
        x0_guess = float(theta[len(theta) // 2])
        popt, _ = curve_fit(
            sigmoid,
            theta,
            S,
            p0=[float(np.max(S)), -5.0, x0_guess],
            maxfev=5000,
        )
        return float(popt[2])
    except RuntimeError:
        logger.warning("Sigmoid fit failed; falling back to steepest_drop")
        return _criticality_steepest_drop(theta, S)


# ---------------------------------------------------------------------------
# Patient-level (time-varying) percolation analysis
# ---------------------------------------------------------------------------


def compute_patient_delta_series(
    df,  # pd.DataFrame
    Q_pop: np.ndarray,
    curve_pop: PercolationCurve,
    cfg: dict[str, Any] | None = None,
    shrinkage_weight: float = 0.3,
) -> pd.DataFrame:
    """Compute time-varying proximity-to-criticality δ(t) for each patient.

    For each patient visit we estimate a patient-specific Q (empirical Bayes
    shrinkage towards Q_pop), build the corresponding graph, compute the GSCC
    at the *population-level* θ_c, and return δ = θ_c - θ_eff where θ_eff
    is the effective threshold at which that patient's GSCC first drops below
    the population GSCC fraction at θ_c.

    A simpler proxy (used when insufficient per-patient data is available):
    δ_proxy(t) = mean(λ_{ij}(t)) - θ_c_pop

    Parameters
    ----------
    df : pd.DataFrame
        Longitudinal data.
    Q_pop : np.ndarray
        Population-level generator matrix.
    curve_pop : PercolationCurve
        Population-level threshold sweep result.
    cfg : dict, optional
    shrinkage_weight : float

    Returns
    -------
    pd.DataFrame
        Columns: patient_id, time_years, delta, mean_intensity,
        gscc_fraction_at_theta_c, percolation_integral_proxy.
    """
    from .data import (
        DEFAULT_PATIENT_COL,
        DEFAULT_TIME_COL,
        DEFAULT_STATE_COL,
        extract_transitions,
    )
    from .msm_fit import empirical_bayes_shrinkage
    from .graph_utils import build_graph

    perc_cfg = (cfg or {}).get("percolation", {})
    method = perc_cfg.get("threshold_method", "absolute")
    starting_states = (cfg or {}).get("states", {}).get("starting_states", [0, 1, 2])
    n_states = Q_pop.shape[0]
    theta_c = curve_pop.theta_c

    records = []

    patients = df[DEFAULT_PATIENT_COL].unique()
    logger.info(
        "Computing patient-level δ(t) for %d patients…", len(patients)
    )

    for pid in patients:
        pat_df = df[df[DEFAULT_PATIENT_COL] == pid].sort_values(DEFAULT_TIME_COL)
        times = pat_df[DEFAULT_TIME_COL].to_numpy()

        # Per-patient Q via empirical Bayes shrinkage
        count_mat, total_t = extract_transitions(pat_df, n_states)
        Q_pat = empirical_bayes_shrinkage(Q_pop, count_mat, total_t, shrinkage_weight)

        G_pat = build_graph(Q_pat)
        # Mean off-diagonal intensity as a scalar summary
        off_diag = Q_pat[~np.eye(n_states, dtype=bool)]
        mean_intensity = float(np.mean(off_diag[off_diag > 0])) if np.any(off_diag > 0) else 0.0

        # GSCC fraction at population θ_c
        G_at_tc = threshold_graph(G_pat, theta=theta_c, method=method)
        gscc_at_tc = len(get_gscc(G_at_tc)) / n_states

        # Percolation integral proxy: AUC of GSCC over population theta grid
        gscc_fracs = []
        for theta in curve_pop.theta_grid:
            G_k = threshold_graph(G_pat, theta=theta, method=method)
            gscc_fracs.append(len(get_gscc(G_k)) / n_states)
        perc_integral_pat = float(np.trapz(gscc_fracs, curve_pop.theta_grid))

        # δ: signed distance from theta_c
        delta_proxy = mean_intensity - theta_c

        for t in times:
            records.append(
                {
                    DEFAULT_PATIENT_COL: pid,
                    DEFAULT_TIME_COL: t,
                    "delta": delta_proxy,
                    "mean_intensity": mean_intensity,
                    "gscc_fraction_at_theta_c": gscc_at_tc,
                    "percolation_integral_proxy": perc_integral_pat,
                }
            )

    result_df = pd.DataFrame(records)
    logger.info("Patient δ(t) series computed: %d rows", len(result_df))
    return result_df


# ---------------------------------------------------------------------------
# Percolation summary statistics
# ---------------------------------------------------------------------------


def compute_summary_statistics(curve: PercolationCurve) -> dict[str, float]:
    """Compute scalar summary statistics from a PercolationCurve.

    Returns a dictionary suitable for logging, CSV export, or model
    augmentation feature construction.

    Parameters
    ----------
    curve : PercolationCurve

    Returns
    -------
    dict with keys:
        theta_c, percolation_integral, susceptibility_integral,
        gscc_at_theta_c, slope_at_theta_c, theta_range,
        n_states, n_starting_states.
    """
    # GSCC fraction at θ_c
    idx_c = int(np.argmin(np.abs(curve.theta_grid - curve.theta_c)))
    gscc_at_tc = float(curve.gscc_fractions[idx_c])

    # Slope at θ_c (finite difference)
    if idx_c > 0 and idx_c < len(curve.theta_grid) - 1:
        slope = float(
            (curve.gscc_fractions[idx_c + 1] - curve.gscc_fractions[idx_c - 1])
            / (curve.theta_grid[idx_c + 1] - curve.theta_grid[idx_c - 1] + 1e-12)
        )
    else:
        slope = 0.0

    theta_range = float(curve.theta_grid[-1] - curve.theta_grid[0])

    stats = {
        "theta_c": float(curve.theta_c),
        "percolation_integral": float(curve.percolation_integral),
        "susceptibility_integral": float(curve.susceptibility_integral),
        "gscc_at_theta_c": gscc_at_tc,
        "slope_at_theta_c": slope,
        "theta_range": theta_range,
        "n_states": int(curve.n_states),
        "n_starting_states": len(curve.starting_states),
    }
    return stats


# ---------------------------------------------------------------------------
# Vulnerability index (Lebkuecher / Lall inspired)
# ---------------------------------------------------------------------------


def percolation_vulnerability_index(
    curves: list[PercolationCurve],
    ref_integral: float | None = None,
) -> np.ndarray:
    """Compute a vulnerability index for a list of PercolationCurves.

    Inspired by Lebkuecher et al. (2024): lower percolation integral →
    more vulnerable network.  Returns standardised (z-scored) vulnerability
    scores (positive = more vulnerable than mean).

    Parameters
    ----------
    curves : list[PercolationCurve]
    ref_integral : float, optional
        Reference (population) integral.  If provided, scores are relative
        to this reference rather than the sample mean.

    Returns
    -------
    np.ndarray, shape (len(curves),)
        Vulnerability indices.  Higher = more fragile.
    """
    integrals = np.array([c.percolation_integral for c in curves])
    # Vulnerability = negative of integral (less area = more vulnerable)
    vuln = -integrals
    if ref_integral is not None:
        vuln = vuln - (-ref_integral)
    else:
        vuln = vuln - float(np.mean(vuln))
    std = float(np.std(vuln))
    if std > 0:
        vuln = vuln / std
    return vuln
