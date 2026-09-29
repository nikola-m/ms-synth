"""
msm_fit.py
==========
Fit a continuous-time multi-state Markov model (CT-MSM) to longitudinal
MS panel data and return the estimated infinitesimal generator matrix Q.

Background
----------
A CT-MSM is characterised by the generator matrix Q (n_states × n_states)
where:
  - Q[i, j] >= 0  for i ≠ j  (transition intensities λ_{ij})
  - Q[i, i]  = -sum_{j≠i} Q[i, j]  (row-sum zero)

The transition probability matrix over interval Δt is:
  P(Δt) = expm(Q · Δt)

Maximum-likelihood estimation maximises the observed-data log-likelihood
  ℓ(Q) = Σ_{i≠j} [ N_{ij} · log λ_{ij}  -  T_i · λ_{ij} ]
where N_{ij} is the count of i→j transitions and T_i is total time spent
in state i.  This has the closed-form solution:
  λ̂_{ij} = N_{ij} / T_i
subject to a non-negativity constraint and optional regularisation.

For more complex panel data (interval-censored observations, arbitrary
time gaps) we also provide an EM / numerical optimisation fallback.

References
----------
Kalbfleisch & Lawless (1985) JASA 80:863-871.
Kannan et al. (2017) Math. Biosci. 289:1-8.
"""

from __future__ import annotations

import logging
from typing import Any

import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize

from .data import (
    extract_transitions,
    DEFAULT_PATIENT_COL,
    DEFAULT_TIME_COL,
    DEFAULT_STATE_COL,
)
from .utils import get_logger

logger = get_logger(__name__)

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def fit_msm(
    df,  # pd.DataFrame – avoid circular import at module level
    n_states: int,
    cfg: dict[str, Any] | None = None,
) -> np.ndarray:
    """Fit a CT-MSM to longitudinal data and return the Q matrix.

    Parameters
    ----------
    df : pd.DataFrame
        Validated longitudinal MS data frame (standard schema).
    n_states : int
        Number of macrostates.
    cfg : dict, optional
        Configuration dict.  Reads ``cfg["msm"]`` for method and hyper-
        parameters.

    Returns
    -------
    Q : np.ndarray, shape (n_states, n_states)
        Estimated infinitesimal generator with rows summing to ≤ 0.

    Raises
    ------
    ValueError
        If estimation fails or produces an invalid generator.
    """
    msm_cfg = (cfg or {}).get("msm", {})
    method = msm_cfg.get("method", "mle").lower()
    regularization = float(msm_cfg.get("regularization", 1e-4))
    allowed = msm_cfg.get("allowed_transitions")  # adjacency mask or None

    count_matrix, total_time = extract_transitions(df, n_states)

    if method == "mle":
        Q = _mle_estimate(count_matrix, total_time, regularization, allowed)
    elif method in ("em", "panel"):
        # "em" is kept as a backwards-compatible alias: the estimator is a
        # direct numerical MLE of the panel likelihood (Kalbfleisch & Lawless
        # 1985), not an EM algorithm.
        Q = _em_estimate(df, n_states, regularization, allowed)
    else:
        raise ValueError(f"Unknown MSM fitting method: '{method}'")

    _validate_Q(Q)
    logger.info(
        "Fitted Q via %s: max off-diag = %.4f, min off-diag = %.4f",
        method,
        float(np.max(Q - np.diag(np.diag(Q)))),
        float(np.min(Q - np.diag(np.diag(Q)))),
    )
    return Q


# ---------------------------------------------------------------------------
# MLE estimator (closed-form)
# ---------------------------------------------------------------------------


def _mle_estimate(
    count_matrix: np.ndarray,
    total_time: np.ndarray,
    regularization: float,
    allowed: list | None,
) -> np.ndarray:
    """Closed-form MLE for the generator matrix.

    λ̂_{ij} = (N_{ij} + ε) / (T_i + ε·n)
    where ε = regularization provides a Laplace-smoothing / ridge penalty
    that avoids division-by-zero for unvisited states.

    Parameters
    ----------
    count_matrix : np.ndarray, shape (n, n)
    total_time : np.ndarray, shape (n,)
    regularization : float
    allowed : list of [i, j] pairs or None

    Returns
    -------
    Q : np.ndarray, shape (n, n)
    """
    n = count_matrix.shape[0]
    eps = regularization
    Q = np.zeros((n, n))

    # Build adjacency mask
    mask = np.ones((n, n), dtype=bool)
    np.fill_diagonal(mask, False)
    if allowed is not None:
        mask[:] = False
        np.fill_diagonal(mask, False)
        for pair in allowed:
            i, j = int(pair[0]), int(pair[1])
            mask[i, j] = True

    for i in range(n):
        Ti = float(total_time[i]) + eps * n
        for j in range(n):
            if i == j or not mask[i, j]:
                continue
            Q[i, j] = (count_matrix[i, j] + eps) / Ti

    # Fix diagonal
    np.fill_diagonal(Q, 0.0)
    np.fill_diagonal(Q, -Q.sum(axis=1))
    return Q


# ---------------------------------------------------------------------------
# EM / numerical optimisation estimator (interval-censored panel data)
# ---------------------------------------------------------------------------


def _em_estimate(
    df,
    n_states: int,
    regularization: float,
    allowed: list | None,
    max_iter: int = 200,
    tol: float = 1e-5,
) -> np.ndarray:
    """Numerical MLE via L-BFGS-B optimising the panel log-likelihood.

    This handles arbitrary observation intervals by computing P(Δt)=expm(Q·Δt)
    at each observed gap.  The parameter vector is the upper-triangular
    (off-diagonal) entries of Q in log-space, keeping Q valid.

    Parameters
    ----------
    df : pd.DataFrame
    n_states : int
    regularization : float
        L2 penalty coefficient on log-intensity parameters.
    allowed : list or None
    max_iter : int
    tol : float

    Returns
    -------
    Q : np.ndarray, shape (n_states, n_states)
    """
    import pandas as pd  # local import to avoid circular
    from .data import DEFAULT_PATIENT_COL, DEFAULT_TIME_COL, DEFAULT_STATE_COL

    logger.info("EM/numerical MSM fitting (n_states=%d)…", n_states)

    # Collect (from_state, to_state, dt) pairs
    transitions = []
    for _, grp in df.groupby(DEFAULT_PATIENT_COL, sort=False):
        grp = grp.sort_values(DEFAULT_TIME_COL)
        states = grp[DEFAULT_STATE_COL].to_numpy(dtype=int)
        times = grp[DEFAULT_TIME_COL].to_numpy(dtype=float)
        for k in range(len(states) - 1):
            dt = float(times[k + 1] - times[k])
            if dt > 0:
                transitions.append((int(states[k]), int(states[k + 1]), dt))

    n = n_states
    # Free parameters: off-diagonal entries (log-space)
    idx_pairs = [(i, j) for i in range(n) for j in range(n) if i != j]
    if allowed is not None:
        allowed_set = {(int(p[0]), int(p[1])) for p in allowed}
        idx_pairs = [p for p in idx_pairs if p in allowed_set]
    n_params = len(idx_pairs)
    pair_to_idx = {p: k for k, p in enumerate(idx_pairs)}

    def _params_to_Q(log_lams: np.ndarray) -> np.ndarray:
        Q_ = np.zeros((n, n))
        for k, (i, j) in enumerate(idx_pairs):
            Q_[i, j] = np.exp(log_lams[k])
        np.fill_diagonal(Q_, 0.0)
        np.fill_diagonal(Q_, -Q_.sum(axis=1))
        return Q_

    a_arr = np.array([t[0] for t in transitions], dtype=int)
    b_arr = np.array([t[1] for t in transitions], dtype=int)
    dt_arr = np.array([t[2] for t in transitions], dtype=float)

    def _neg_loglik(log_lams: np.ndarray) -> float:
        Q_ = _params_to_Q(log_lams)
        probs = panel_transition_probs(Q_, a_arr, b_arr, dt_arr)
        nll = -float(np.sum(np.log(np.clip(probs, 1e-15, None))))
        # L2 regularisation
        nll += regularization * float(np.sum(log_lams ** 2))
        return nll

    # Initialise at the crude count/time estimate (faster convergence)
    cnt = np.zeros((n, n)); tot = np.zeros(n)
    np.add.at(cnt, (a_arr, b_arr), 1.0); np.add.at(tot, a_arr, dt_arr)
    x0 = np.array([np.log(max(cnt[i, j] / max(tot[i], 1e-9), 1e-3))
                   for (i, j) in idx_pairs])
    result = minimize(
        _neg_loglik,
        x0,
        method="L-BFGS-B",
        options={"maxiter": max_iter, "ftol": tol * 1e-3},
    )
    if not result.success:
        logger.warning("EM optimiser did not fully converge: %s", result.message)

    Q = _params_to_Q(result.x)
    return Q



# ---------------------------------------------------------------------------
# Fast panel likelihood and known-truth estimators
# ---------------------------------------------------------------------------

def panel_transition_probs(
    Q: np.ndarray, a: np.ndarray, b: np.ndarray, dt: np.ndarray
) -> np.ndarray:
    """Vectorised P(dt)[a, b] = expm(Q*dt)[a, b] for many observation pairs.

    Uses the eigendecomposition Q = V diag(w) V^-1, so that
    P(dt)[a, b] = sum_k V[a, k] exp(w_k dt) V^-1[k, b]; falls back to a
    batched matrix exponential when V is ill-conditioned (near-defective Q).
    """
    a = np.asarray(a, dtype=int); b = np.asarray(b, dtype=int)
    dt = np.asarray(dt, dtype=float)
    try:
        w, V = np.linalg.eig(Q)
        if np.linalg.cond(V) > 1e8:
            raise np.linalg.LinAlgError("ill-conditioned eigenbasis")
        Vinv = np.linalg.inv(V)
        E = np.exp(np.outer(dt, w))                      # (m, n)
        probs = np.einsum("mk,mk,km->m", V[a, :], E, Vinv[:, b]).real
    except np.linalg.LinAlgError:
        P = expm(Q[None, :, :] * dt[:, None, None])      # (m, n, n)
        probs = P[np.arange(len(dt)), a, b]
    return probs


def estimate_Q_crude(
    df, n_states: int, time_col: str = "disease_duration_yr",
    state_col: str = "state", patient_col: str = "patient_id",
) -> np.ndarray:
    """Crude estimator: lambda_ij = N_ij / T_i from consecutive visits.

    This is the MLE only for *continuously observed* paths. Applied to
    panel data it ignores unobserved intermediate jumps and is biased
    towards zero for fast, reversible transitions.
    """
    n = n_states
    cnt = np.zeros((n, n)); tot = np.zeros(n)
    for _, g in df.sort_values([patient_col, time_col]).groupby(patient_col, sort=False):
        s_ = g[state_col].to_numpy(dtype=int); t_ = g[time_col].to_numpy(dtype=float)
        for k in range(len(s_) - 1):
            cnt[s_[k], s_[k + 1]] += 1.0
            tot[s_[k]] += max(t_[k + 1] - t_[k], 1e-9)
    Q = np.zeros((n, n))
    for i in range(n):
        if tot[i] > 0:
            Q[i] = cnt[i] / tot[i]
    np.fill_diagonal(Q, 0.0); np.fill_diagonal(Q, -Q.sum(axis=1))
    return Q


def estimate_Q_panel(
    df, n_states: int, allowed: list | None = None,
    time_col: str = "disease_duration_yr", state_col: str = "state",
    patient_col: str = "patient_id", regularization: float = 0.0,
    max_iter: int = 500,
) -> np.ndarray:
    """Panel-data MLE (Kalbfleisch & Lawless 1985) of Q.

    Maximises sum log expm(Q dt)[s_k, s_{k+1}] over consecutive visits,
    optionally restricted to an ``allowed`` list of (i, j) pairs.
    """
    tmp = df.rename(columns={patient_col: DEFAULT_PATIENT_COL,
                             time_col: DEFAULT_TIME_COL,
                             state_col: DEFAULT_STATE_COL})
    return _em_estimate(tmp, n_states, regularization, allowed, max_iter=max_iter)


def estimate_Q_from_paths(latent: dict, n_states: int) -> np.ndarray:
    """Oracle estimator from exactly observed latent paths.

    ``latent`` maps patient id -> {"events": [(t, s), ...], "follow_up": T},
    as returned by ``generate_synthetic_ms_dataset(..., return_latent=True)``.
    With complete paths N_ij / T_i is the exact MLE, so its error reflects
    only finite-sample variability and between-patient heterogeneity.
    """
    n = n_states
    cnt = np.zeros((n, n)); tot = np.zeros(n)
    for rec in latent.values():
        ev = rec["events"]; T = float(rec["follow_up"])
        for k, (t_k, s_k) in enumerate(ev):
            if t_k >= T:
                break
            t_next = ev[k + 1][0] if k + 1 < len(ev) else np.inf
            tot[s_k] += min(t_next, T) - t_k
            if t_next <= T:
                cnt[s_k, ev[k + 1][1]] += 1.0
    Q = np.zeros((n, n))
    for i in range(n):
        if tot[i] > 0:
            Q[i] = cnt[i] / tot[i]
    np.fill_diagonal(Q, 0.0); np.fill_diagonal(Q, -Q.sum(axis=1))
    return Q

# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------


def _validate_Q(Q: np.ndarray) -> None:
    """Assert that Q is a valid infinitesimal generator.

    Parameters
    ----------
    Q : np.ndarray

    Raises
    ------
    ValueError
    """
    n = Q.shape[0]
    if Q.shape != (n, n):
        raise ValueError("Q must be square")
    if np.any(Q < -1e-6):
        neg_entries = np.argwhere(Q < -1e-6)
        # Off-diagonal negatives are invalid; diagonal is expected negative
        off_neg = [(i, j) for (i, j) in neg_entries if i != j]
        if off_neg:
            raise ValueError(
                f"Q has negative off-diagonal entries at positions: {off_neg[:5]}"
            )
    row_sums = Q.sum(axis=1)
    if np.any(np.abs(row_sums) > 1e-4):
        logger.warning("Q row sums are not all zero (max |sum| = %.2e)", np.max(np.abs(row_sums)))


def Q_to_transition_matrix(Q: np.ndarray, dt: float) -> np.ndarray:
    """Compute P(dt) = expm(Q · dt).

    Parameters
    ----------
    Q : np.ndarray, shape (n, n)
    dt : float

    Returns
    -------
    P : np.ndarray, shape (n, n)
    """
    P = expm(Q * dt)
    P = np.clip(P, 0.0, 1.0)
    # Re-normalise rows
    row_sums = P.sum(axis=1, keepdims=True)
    row_sums = np.where(row_sums > 0, row_sums, 1.0)
    return P / row_sums


def empirical_bayes_shrinkage(
    Q_pop: np.ndarray,
    count_matrix: np.ndarray,
    total_time: np.ndarray,
    shrinkage_weight: float = 0.3,
) -> np.ndarray:
    """Shrink a patient-specific Q towards the population Q.

    Parameters
    ----------
    Q_pop : np.ndarray
        Population-level generator matrix.
    count_matrix : np.ndarray, shape (n, n)
        Patient-level transition counts.
    total_time : np.ndarray, shape (n,)
        Patient-level time in each state.
    shrinkage_weight : float
        Weight in [0, 1] on the population estimate (0 = no shrinkage).

    Returns
    -------
    Q_shrunk : np.ndarray, shape (n, n)
    """
    n = Q_pop.shape[0]
    eps = 1e-6
    Q_pat = np.zeros((n, n))
    for i in range(n):
        Ti = max(float(total_time[i]), eps)
        for j in range(n):
            if i != j:
                Q_pat[i, j] = float(count_matrix[i, j]) / Ti
    np.fill_diagonal(Q_pat, 0.0)
    np.fill_diagonal(Q_pat, -Q_pat.sum(axis=1))

    Q_shrunk = (1.0 - shrinkage_weight) * Q_pat + shrinkage_weight * Q_pop
    np.fill_diagonal(Q_shrunk, 0.0)
    np.fill_diagonal(Q_shrunk, -Q_shrunk.sum(axis=1))
    return Q_shrunk
