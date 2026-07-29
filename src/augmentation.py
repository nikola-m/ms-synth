"""
augmentation.py
===============
Embed percolation-derived features back into the multi-state Markov model
(CT-MSM) and into trajectory analysis models.

Two augmentation strategies are provided:

1. **Percolation-Augmented MSM (PA-MSM)**
   Selected transition intensities in Q are made to depend on the
   percolation proximity-to-criticality covariate δ:

     λ_{ij}^{aug}(δ) = λ_{ij} · f(δ)

   where f is a monotone coupling function:
   - "logistic":  f(δ) = 1 / (1 + exp(-k · δ))
     (reversible edges → 0 as δ → -∞; progressive edges unaffected)
   - "exp_decay": f(δ) = exp(-k · max(-δ, 0))

2. **Latent-Class / Mixture Trajectory Model**
   Uses percolation integral and susceptibility integral as time-varying
   predictors of EDSS trajectory group membership (via multinomial logistic
   regression or k-means clustering on the percolation feature time series).

Model comparison (AIC / BIC) between baseline Q and augmented Q is also
implemented using the observed-data panel log-likelihood.

References
----------
Kannan et al. (2017) Math. Biosci. 289:1-8.
Mirkov et al. (2026) arXiv:2602.08576.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
from scipy.linalg import expm
from scipy.optimize import minimize

from .percolation import PercolationCurve
from .utils import get_logger

logger = get_logger(__name__)


# ---------------------------------------------------------------------------
# Percolation-Augmented MSM
# ---------------------------------------------------------------------------


def augment_Q(
    Q: np.ndarray,
    delta: float,
    cfg: dict[str, Any] | None = None,
) -> np.ndarray:
    """Return an augmented Q matrix where reversible intensities are scaled
    by a percolation-covariate coupling function evaluated at δ.

    The coupling attenuates reversible (remission) edges as the system
    approaches the percolation threshold (δ → 0) and beyond (δ < 0),
    while leaving progressive (worsening) edges unchanged.  This
    operationalises the Kannan et al. (2017) idea that crossing the
    immune/CNS threshold severs the negative-feedback loop.

    Parameters
    ----------
    Q : np.ndarray, shape (n, n)
        Population-level generator.
    delta : float
        Proximity to criticality (θ_c - θ_current).  Positive = sub-critical
        (reversible dynamics intact); negative = super-critical (fragmented).
    cfg : dict, optional
        Reads ``cfg["augmentation"]`` for coupling parameters.

    Returns
    -------
    Q_aug : np.ndarray, shape (n, n)
        Augmented generator with modified reversible intensities.
    """
    aug_cfg = (cfg or {}).get("augmentation", {})
    coupling_form = aug_cfg.get("coupling_form", "logistic")
    k = float(aug_cfg.get("coupling_strength", 2.0))
    modulate = aug_cfg.get("modulate_edges", "reversible")

    n = Q.shape[0]
    Q_aug = Q.copy()

    # Identify reversible edges: i→j such that Q[j, i] > 0
    for i in range(n):
        for j in range(n):
            if i == j or Q[i, j] <= 0:
                continue
            is_reversible = Q[j, i] > 0

            if modulate == "reversible" and not is_reversible:
                continue  # leave progressive edges untouched
            if modulate == "progressive" and is_reversible:
                continue

            # Apply coupling
            f = _coupling_function(delta, k, coupling_form)
            Q_aug[i, j] = Q[i, j] * f

    # Re-fix diagonal
    np.fill_diagonal(Q_aug, 0.0)
    np.fill_diagonal(Q_aug, -Q_aug.sum(axis=1))
    return Q_aug


def _coupling_function(delta: float, k: float, form: str) -> float:
    """Compute the coupling factor f(δ) ∈ [0, 1].

    Parameters
    ----------
    delta : float
        Proximity to criticality.
    k : float
        Steepness / rate.
    form : str
        "logistic" or "exp_decay".

    Returns
    -------
    float in [0, 1]
    """
    if form == "logistic":
        return float(1.0 / (1.0 + np.exp(-k * delta)))
    elif form == "exp_decay":
        return float(np.exp(-k * max(-delta, 0.0)))
    else:
        raise ValueError(f"Unknown coupling form: '{form}'")


# ---------------------------------------------------------------------------
# Panel log-likelihood (for AIC/BIC)
# ---------------------------------------------------------------------------


def panel_loglikelihood(
    Q: np.ndarray,
    df,  # pd.DataFrame
) -> float:
    """Compute the panel-data log-likelihood for a given generator Q.

    ℓ(Q) = Σ_transitions  log P(s_{t+1} | s_t, Δt; Q)
    where P(Δt) = expm(Q · Δt).

    Parameters
    ----------
    Q : np.ndarray
    df : pd.DataFrame

    Returns
    -------
    float
    """
    from .data import DEFAULT_PATIENT_COL, DEFAULT_TIME_COL, DEFAULT_STATE_COL

    ll = 0.0
    for _, grp in df.groupby(DEFAULT_PATIENT_COL, sort=False):
        grp = grp.sort_values(DEFAULT_TIME_COL)
        states = grp[DEFAULT_STATE_COL].to_numpy(dtype=int)
        times = grp[DEFAULT_TIME_COL].to_numpy(dtype=float)
        for k in range(len(states) - 1):
            dt = float(times[k + 1] - times[k])
            if dt <= 0:
                continue
            P = expm(Q * dt)
            P = np.clip(P, 1e-15, 1.0)
            s_from, s_to = int(states[k]), int(states[k + 1])
            ll += np.log(P[s_from, s_to])
    return float(ll)


def aic_bic(
    Q: np.ndarray,
    df,  # pd.DataFrame
    n_free_params: int | None = None,
) -> dict[str, float]:
    """Compute AIC and BIC for a fitted MSM.

    Parameters
    ----------
    Q : np.ndarray
    df : pd.DataFrame
    n_free_params : int, optional
        Number of free parameters.  Defaults to number of non-zero
        off-diagonal entries.

    Returns
    -------
    dict with keys: loglik, n_params, aic, bic.
    """
    n = Q.shape[0]
    if n_free_params is None:
        n_free_params = int(np.sum(Q > 1e-9)) - n  # off-diag non-zero

    from .data import DEFAULT_PATIENT_COL
    n_obs = len(df) - df[DEFAULT_PATIENT_COL].nunique()  # number of transitions

    ll = panel_loglikelihood(Q, df)
    aic = -2.0 * ll + 2.0 * n_free_params
    bic = -2.0 * ll + np.log(max(n_obs, 1)) * n_free_params

    metrics = {
        "loglik": ll,
        "n_params": n_free_params,
        "aic": aic,
        "bic": bic,
    }
    logger.debug("AIC=%.2f  BIC=%.2f  LL=%.2f  k=%d", aic, bic, ll, n_free_params)
    return metrics


def compare_models(
    Q_baseline: np.ndarray,
    Q_augmented: np.ndarray,
    df,  # pd.DataFrame
    n_extra_params: int = 1,
) -> dict[str, Any]:
    """Compare baseline and augmented MSM via AIC, BIC, and LR test.

    Parameters
    ----------
    Q_baseline : np.ndarray
    Q_augmented : np.ndarray
    df : pd.DataFrame
    n_extra_params : int
        Number of additional parameters in the augmented model (e.g., the
        coupling strength k).

    Returns
    -------
    dict with keys: baseline, augmented, delta_aic, delta_bic, lr_stat,
    lr_pvalue, preferred.
    """
    from scipy.stats import chi2

    base_metrics = aic_bic(Q_baseline, df)
    aug_metrics = aic_bic(Q_augmented, df, n_free_params=base_metrics["n_params"] + n_extra_params)

    delta_aic = float(base_metrics["aic"] - aug_metrics["aic"])
    delta_bic = float(base_metrics["bic"] - aug_metrics["bic"])

    lr_stat = float(2.0 * (aug_metrics["loglik"] - base_metrics["loglik"]))
    lr_pvalue = float(1.0 - chi2.cdf(max(lr_stat, 0.0), df=n_extra_params))

    preferred = "augmented" if delta_aic > 0 else "baseline"

    result = {
        "baseline": base_metrics,
        "augmented": aug_metrics,
        "delta_aic": delta_aic,
        "delta_bic": delta_bic,
        "lr_stat": lr_stat,
        "lr_pvalue": lr_pvalue,
        "preferred": preferred,
    }
    logger.info(
        "Model comparison: ΔAIC=%.2f  ΔBIC=%.2f  LR=%.2f (p=%.4f)  preferred=%s",
        delta_aic, delta_bic, lr_stat, lr_pvalue, preferred,
    )
    return result


# ---------------------------------------------------------------------------
# Trajectory clustering with percolation features
# ---------------------------------------------------------------------------


def cluster_trajectories_by_percolation(
    patient_delta_df: pd.DataFrame,
    n_clusters: int = 3,
    feature_cols: list[str] | None = None,
    seed: int = 42,
) -> pd.DataFrame:
    """Cluster patients by their percolation feature time series.

    Uses k-means on per-patient summary statistics of the δ(t) and
    percolation_integral_proxy series to identify distinct trajectory
    groups (analogous to latent-class growth analysis).

    Parameters
    ----------
    patient_delta_df : pd.DataFrame
        Output of ``percolation.compute_patient_delta_series``.
    n_clusters : int
    feature_cols : list[str], optional
        Columns to use as clustering features.  Defaults to
        ["delta", "mean_intensity", "percolation_integral_proxy"].
    seed : int

    Returns
    -------
    pd.DataFrame
        patient_delta_df with an added ``trajectory_cluster`` column.
    """
    from sklearn.cluster import KMeans
    from sklearn.preprocessing import StandardScaler
    from .data import DEFAULT_PATIENT_COL, DEFAULT_TIME_COL

    if feature_cols is None:
        feature_cols = ["delta", "mean_intensity", "percolation_integral_proxy"]

    # Summarise per patient: mean and std of each feature
    agg_dict = {}
    for col in feature_cols:
        if col in patient_delta_df.columns:
            agg_dict[col + "_mean"] = (col, "mean")
            agg_dict[col + "_std"] = (col, "std")

    patient_features = (
        patient_delta_df.groupby(DEFAULT_PATIENT_COL)
        .agg(**agg_dict)
        .reset_index()
        .fillna(0.0)
    )

    feat_cols = [c for c in patient_features.columns if c != DEFAULT_PATIENT_COL]
    X = patient_features[feat_cols].to_numpy()

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    km = KMeans(n_clusters=n_clusters, random_state=seed, n_init="auto")
    labels = km.fit_predict(X_scaled)
    patient_features["trajectory_cluster"] = labels

    result = patient_delta_df.merge(
        patient_features[[DEFAULT_PATIENT_COL, "trajectory_cluster"]],
        on=DEFAULT_PATIENT_COL,
        how="left",
    )

    logger.info(
        "Trajectory clustering: %d clusters, distribution %s",
        n_clusters,
        dict(pd.Series(labels).value_counts()),
    )
    return result


def fit_augmented_msm(
    Q_init: np.ndarray,
    df,
    curve_pop: PercolationCurve,
    cfg: dict[str, Any] | None = None,
) -> tuple[np.ndarray, dict[str, float]]:
    """Optimise the coupling strength k of the PA-MSM via panel MLE.

    Holding the shape of Q fixed, searches for the coupling parameter k
    that maximises the panel log-likelihood of the PA-MSM.

    Parameters
    ----------
    Q_init : np.ndarray
        Baseline (unadjusted) generator.
    df : pd.DataFrame
    curve_pop : PercolationCurve
        Population-level percolation results (provides θ_c and δ values).
    cfg : dict, optional

    Returns
    -------
    Q_aug_opt : np.ndarray
        Augmented Q at the optimal coupling strength.
    fit_info : dict
        Contains: optimal_k, loglik, aic, bic.
    """
    aug_cfg = (cfg or {}).get("augmentation", {})
    coupling_form = aug_cfg.get("coupling_form", "logistic")
    delta_pop = 0.0  # use population-level delta at theta_c (by definition 0)

    def _neg_ll(log_k: np.ndarray) -> float:
        k_val = float(np.exp(log_k[0]))
        Q_a = augment_Q(Q_init, delta=delta_pop, cfg={
            "augmentation": {
                "coupling_form": coupling_form,
                "coupling_strength": k_val,
                "modulate_edges": aug_cfg.get("modulate_edges", "reversible"),
            }
        })
        return -panel_loglikelihood(Q_a, df)

    result = minimize(
        _neg_ll,
        x0=np.array([np.log(2.0)]),
        method="L-BFGS-B",
        bounds=[(-5.0, 5.0)],
    )
    k_opt = float(np.exp(result.x[0]))
    aug_cfg_opt = dict(aug_cfg)
    aug_cfg_opt["coupling_strength"] = k_opt

    Q_aug_opt = augment_Q(Q_init, delta=delta_pop, cfg={"augmentation": aug_cfg_opt})
    metrics = aic_bic(Q_aug_opt, df, n_free_params=int(np.sum(Q_init > 1e-9)) + 1)
    metrics["optimal_k"] = k_opt

    logger.info("PA-MSM fitted: optimal_k=%.4f, AIC=%.2f, BIC=%.2f", k_opt, metrics["aic"], metrics["bic"])
    return Q_aug_opt, metrics
