"""
visualization.py
================
Publication-quality figures for the MS percolation framework.

All plotting functions accept a ``save_path`` parameter and an ``ax``
parameter (if None, a new figure is created).

Figures produced
----------------
plot_percolation_curve    : GSCC fraction vs θ with θ_c annotation.
plot_component_sizes      : GSCC / in / out component sizes vs θ.
plot_susceptibility       : Susceptibility χ(θ) with integral shading.
plot_network_snapshot     : DiGraph at a given θ (networkx layout).
plot_patient_delta        : Time-series of δ(t) for selected patients.
plot_model_comparison     : Bar chart of AIC/BIC for baseline vs augmented.
plot_phase_diagram        : 2-D grid of θ_c values across parameter sweeps
                            (analogous to Kannan et al. 2017 Fig. 7).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import networkx as nx
import numpy as np
import pandas as pd

from .percolation import PercolationCurve
from .utils import get_logger, ensure_dir

logger = get_logger(__name__)

# Use non-interactive backend if running headlessly
matplotlib.use("Agg")

_PALETTE = {
    "GSCC": "#2563EB",
    "in_comp": "#16A34A",
    "out_comp": "#DC2626",
    "susceptibility": "#9333EA",
    "theta_c": "#F59E0B",
    "RRMS": "#3B82F6",
    "SPMS": "#EF4444",
    "PPMS": "#F97316",
    "PRMS": "#8B5CF6",
}


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------


def _savefig(fig: plt.Figure, save_path: str | Path | None, dpi: int = 150) -> None:
    if save_path is not None:
        save_path = Path(save_path)
        ensure_dir(save_path.parent)
        fig.savefig(save_path, dpi=dpi, bbox_inches="tight")
        logger.info("Figure saved: %s", save_path)
    plt.close(fig)


# ---------------------------------------------------------------------------
# 1. Percolation curve (GSCC fraction vs θ)
# ---------------------------------------------------------------------------


def plot_percolation_curve(
    curve: PercolationCurve,
    ax: plt.Axes | None = None,
    save_path: str | Path | None = None,
    title: str = "Directed Percolation Curve (GSCC)",
    dpi: int = 150,
) -> plt.Axes:
    """Plot GSCC fraction S(θ) vs threshold θ with θ_c annotation.

    Parameters
    ----------
    curve : PercolationCurve
    ax : matplotlib.axes.Axes, optional
    save_path : str or Path, optional
    title : str
    dpi : int

    Returns
    -------
    matplotlib.axes.Axes
    """
    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=(8, 5))
    else:
        fig = ax.get_figure()

    ax.plot(
        curve.theta_grid,
        curve.gscc_fractions,
        color=_PALETTE["GSCC"],
        lw=2.0,
        label="GSCC fraction S(θ)",
    )
    # Shade area under curve (percolation integral)
    ax.fill_between(
        curve.theta_grid,
        curve.gscc_fractions,
        alpha=0.15,
        color=_PALETTE["GSCC"],
        label=f"I_perc = {curve.percolation_integral:.3f}",
    )
    # Mark θ_c
    ax.axvline(
        curve.theta_c,
        color=_PALETTE["theta_c"],
        lw=1.5,
        ls="--",
        label=f"θ_c = {curve.theta_c:.4f}",
    )
    ax.set_xlabel("Threshold θ  (transition intensity)", fontsize=12)
    ax.set_ylabel("Normalised GSCC size  S(θ) = |GSCC| / N", fontsize=12)
    ax.set_title(title, fontsize=13)
    ax.legend(fontsize=10)
    ax.set_xlim(curve.theta_grid[0], curve.theta_grid[-1])
    ax.set_ylim(-0.02, 1.05)
    ax.grid(True, alpha=0.3)

    if standalone:
        fig.tight_layout()
        _savefig(fig, save_path, dpi=dpi)
    return ax


# ---------------------------------------------------------------------------
# 2. Component sizes (GSCC / in / out)
# ---------------------------------------------------------------------------


def plot_component_sizes(
    curve: PercolationCurve,
    ax: plt.Axes | None = None,
    save_path: str | Path | None = None,
    dpi: int = 150,
) -> plt.Axes:
    """Plot absolute GSCC, in-component, and out-component sizes vs θ.

    Parameters
    ----------
    curve : PercolationCurve
    ax : matplotlib.axes.Axes, optional
    save_path : str or Path, optional
    dpi : int

    Returns
    -------
    matplotlib.axes.Axes
    """
    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=(8, 5))
    else:
        fig = ax.get_figure()

    ax.plot(curve.theta_grid, curve.gscc_sizes, color=_PALETTE["GSCC"],
            lw=2, label="GSCC")
    ax.plot(curve.theta_grid, curve.in_comp_sizes, color=_PALETTE["in_comp"],
            lw=1.5, ls="--", label="In-component")
    ax.plot(curve.theta_grid, curve.out_comp_sizes, color=_PALETTE["out_comp"],
            lw=1.5, ls=":", label="Out-component")
    ax.axvline(curve.theta_c, color=_PALETTE["theta_c"], lw=1.5, ls="--",
               label=f"θ_c = {curve.theta_c:.4f}")
    ax.set_xlabel("Threshold θ", fontsize=12)
    ax.set_ylabel("Component size (nodes)", fontsize=12)
    ax.set_title("Directed Component Sizes vs Threshold", fontsize=13)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    if standalone:
        fig.tight_layout()
        _savefig(fig, save_path, dpi=dpi)
    return ax


# ---------------------------------------------------------------------------
# 3. Susceptibility
# ---------------------------------------------------------------------------


def plot_susceptibility(
    curve: PercolationCurve,
    ax: plt.Axes | None = None,
    save_path: str | Path | None = None,
    dpi: int = 150,
) -> plt.Axes:
    """Plot susceptibility χ(θ) = Var(SCC sizes) vs θ.

    Parameters
    ----------
    curve : PercolationCurve
    ax, save_path, dpi : see other plot functions.

    Returns
    -------
    matplotlib.axes.Axes
    """
    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=(8, 4))
    else:
        fig = ax.get_figure()

    ax.plot(curve.theta_grid, curve.susceptibility, color=_PALETTE["susceptibility"],
            lw=2, label="χ(θ) = Var(SCC sizes)")
    ax.fill_between(curve.theta_grid, curve.susceptibility, alpha=0.2,
                    color=_PALETTE["susceptibility"],
                    label=f"∫χ dθ = {curve.susceptibility_integral:.4f}")
    ax.axvline(curve.theta_c, color=_PALETTE["theta_c"], lw=1.5, ls="--",
               label=f"θ_c = {curve.theta_c:.4f}")
    ax.set_xlabel("Threshold θ", fontsize=12)
    ax.set_ylabel("Susceptibility χ(θ)", fontsize=12)
    ax.set_title("Percolation Susceptibility", fontsize=13)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    if standalone:
        fig.tight_layout()
        _savefig(fig, save_path, dpi=dpi)
    return ax


# ---------------------------------------------------------------------------
# 4. Network snapshot at a given θ
# ---------------------------------------------------------------------------


def plot_network_snapshot(
    G,  # nx.DiGraph
    theta: float,
    theta_method: str = "absolute",
    state_labels: dict[int, str] | None = None,
    title: str | None = None,
    ax: plt.Axes | None = None,
    save_path: str | Path | None = None,
    dpi: int = 150,
    seed: int = 42,
) -> plt.Axes:
    """Draw the directed graph G after thresholding at θ.

    Nodes in the GSCC are highlighted.  Edge width encodes weight.

    Parameters
    ----------
    G : nx.DiGraph
    theta : float
    theta_method : str
    state_labels : dict[int, str], optional
    title : str, optional
    ax, save_path, dpi, seed : standard.

    Returns
    -------
    matplotlib.axes.Axes
    """
    from .graph_utils import threshold_graph, get_gscc

    G_t = threshold_graph(G, theta=theta, method=theta_method)
    gscc = get_gscc(G_t)

    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=(9, 7))
    else:
        fig = ax.get_figure()

    pos = nx.spring_layout(G_t, seed=seed)

    node_colors = [
        _PALETTE["GSCC"] if n in gscc else "#94A3B8"
        for n in G_t.nodes()
    ]
    labels = {n: (state_labels.get(n, str(n)) if state_labels else str(n)) for n in G_t.nodes()}

    weights = [d["weight"] for _, _, d in G_t.edges(data=True)]
    w_max = max(weights) if weights else 1.0
    edge_widths = [2.0 * w / max(w_max, 1e-9) for w in weights]

    nx.draw_networkx_nodes(G_t, pos, node_color=node_colors, node_size=300, ax=ax)
    nx.draw_networkx_labels(G_t, pos, labels=labels, font_size=7, ax=ax)
    nx.draw_networkx_edges(
        G_t, pos, width=edge_widths, edge_color="#64748B",
        arrows=True, arrowsize=10, ax=ax,
        connectionstyle="arc3,rad=0.1",
    )

    if title is None:
        title = f"State-Transition Graph at θ = {theta:.4f}  |GSCC|={len(gscc)}"
    ax.set_title(title, fontsize=12)
    ax.axis("off")

    # Legend patches
    import matplotlib.patches as mpatches
    gscc_patch = mpatches.Patch(color=_PALETTE["GSCC"], label="GSCC node")
    other_patch = mpatches.Patch(color="#94A3B8", label="Non-GSCC node")
    ax.legend(handles=[gscc_patch, other_patch], loc="lower right", fontsize=9)

    if standalone:
        fig.tight_layout()
        _savefig(fig, save_path, dpi=dpi)
    return ax


# ---------------------------------------------------------------------------
# 5. Patient δ(t) trajectories
# ---------------------------------------------------------------------------


def plot_patient_delta(
    patient_delta_df: pd.DataFrame,
    patient_ids: list | None = None,
    theta_c: float | None = None,
    ax: plt.Axes | None = None,
    save_path: str | Path | None = None,
    dpi: int = 150,
    max_patients: int = 20,
) -> plt.Axes:
    """Plot time-series of proximity-to-criticality δ(t) per patient.

    Parameters
    ----------
    patient_delta_df : pd.DataFrame
        Output of ``percolation.compute_patient_delta_series``.
    patient_ids : list, optional
        Subset of patient IDs to plot.  Defaults to first ``max_patients``.
    theta_c : float, optional
        Horizontal reference line at δ = 0 (criticality).
    ax, save_path, dpi, max_patients : standard.

    Returns
    -------
    matplotlib.axes.Axes
    """
    from .data import DEFAULT_PATIENT_COL, DEFAULT_TIME_COL

    if patient_ids is None:
        patient_ids = list(patient_delta_df[DEFAULT_PATIENT_COL].unique())[:max_patients]

    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=(10, 5))
    else:
        fig = ax.get_figure()

    cmap = cm.get_cmap("tab20", len(patient_ids))

    for k, pid in enumerate(patient_ids):
        sub = patient_delta_df[patient_delta_df[DEFAULT_PATIENT_COL] == pid]
        sub = sub.sort_values(DEFAULT_TIME_COL)
        color = _PALETTE.get("RRMS", cmap(k))
        if "trajectory_cluster" in sub.columns:
            cluster_colors = {0: _PALETTE["RRMS"], 1: _PALETTE["SPMS"], 2: _PALETTE["PPMS"]}
            cl = int(sub["trajectory_cluster"].iloc[0])
            color = cluster_colors.get(cl, cmap(k))
        ax.plot(sub[DEFAULT_TIME_COL], sub["delta"], lw=1.0, alpha=0.6, color=color)

    ax.axhline(0.0, color=_PALETTE["theta_c"], lw=1.5, ls="--",
               label="δ = 0  (critical threshold)")
    ax.set_xlabel("Time (years)", fontsize=12)
    ax.set_ylabel("Proximity to criticality  δ(t) = θ_c − θ_eff(t)", fontsize=12)
    ax.set_title("Patient-Level Proximity to Percolation Threshold", fontsize=13)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    if standalone:
        fig.tight_layout()
        _savefig(fig, save_path, dpi=dpi)
    return ax


# ---------------------------------------------------------------------------
# 6. Model comparison bar chart
# ---------------------------------------------------------------------------


def plot_model_comparison(
    comparison_result: dict[str, Any],
    ax: plt.Axes | None = None,
    save_path: str | Path | None = None,
    dpi: int = 150,
) -> plt.Axes:
    """Bar chart comparing AIC / BIC of baseline vs augmented MSM.

    Parameters
    ----------
    comparison_result : dict
        Output of ``augmentation.compare_models``.
    ax, save_path, dpi : standard.

    Returns
    -------
    matplotlib.axes.Axes
    """
    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=(7, 5))
    else:
        fig = ax.get_figure()

    metrics = ["aic", "bic"]
    x = np.arange(len(metrics))
    w = 0.35

    base_vals = [comparison_result["baseline"][m] for m in metrics]
    aug_vals = [comparison_result["augmented"][m] for m in metrics]

    ax.bar(x - w / 2, base_vals, width=w, label="Baseline MSM", color="#64748B")
    ax.bar(x + w / 2, aug_vals, width=w, label="PA-MSM (augmented)", color=_PALETTE["GSCC"])

    ax.set_xticks(x)
    ax.set_xticklabels(["AIC", "BIC"], fontsize=12)
    ax.set_ylabel("Information criterion (lower = better)", fontsize=11)
    ax.set_title(
        f"Model comparison  (ΔAIC={comparison_result['delta_aic']:.1f}, "
        f"p={comparison_result['lr_pvalue']:.3f})",
        fontsize=12,
    )
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, axis="y")

    if standalone:
        fig.tight_layout()
        _savefig(fig, save_path, dpi=dpi)
    return ax


# ---------------------------------------------------------------------------
# 7. Summary dashboard (all panels)
# ---------------------------------------------------------------------------


def plot_summary_dashboard(
    curve: PercolationCurve,
    G,
    patient_delta_df: pd.DataFrame | None = None,
    comparison_result: dict[str, Any] | None = None,
    save_path: str | Path | None = None,
    dpi: int = 150,
) -> plt.Figure:
    """Produce a 2×3 summary dashboard figure.

    Panels:
      [0,0] Percolation curve          [0,1] Component sizes
      [0,2] Susceptibility             [1,0] Network at sub-critical θ
      [1,1] Network at super-critical  [1,2] δ(t) or model comparison

    Parameters
    ----------
    curve : PercolationCurve
    G : nx.DiGraph
    patient_delta_df : pd.DataFrame, optional
    comparison_result : dict, optional
    save_path : str or Path, optional
    dpi : int

    Returns
    -------
    matplotlib.figure.Figure
    """
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    fig.suptitle(
        "MS Percolation Framework – Summary Dashboard",
        fontsize=15,
        fontweight="bold",
    )

    plot_percolation_curve(curve, ax=axes[0, 0])
    plot_component_sizes(curve, ax=axes[0, 1])
    plot_susceptibility(curve, ax=axes[0, 2])

    # Sub-critical snapshot (θ = θ_c * 0.5)
    theta_sub = 0.5 * curve.theta_c
    plot_network_snapshot(G, theta=theta_sub, ax=axes[1, 0],
                          title=f"Sub-critical: θ = {theta_sub:.3f}")

    # Super-critical snapshot (θ = θ_c * 1.5)
    theta_sup = 1.5 * curve.theta_c
    plot_network_snapshot(G, theta=theta_sup, ax=axes[1, 1],
                          title=f"Super-critical: θ = {theta_sup:.3f}")

    if comparison_result is not None:
        plot_model_comparison(comparison_result, ax=axes[1, 2])
    elif patient_delta_df is not None:
        plot_patient_delta(patient_delta_df, ax=axes[1, 2])
    else:
        axes[1, 2].text(0.5, 0.5, "No data", ha="center", va="center")
        axes[1, 2].axis("off")

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    _savefig(fig, save_path, dpi=dpi)
    return fig
