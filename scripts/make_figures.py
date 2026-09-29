#!/usr/bin/env python3
"""
Regenerate every manuscript figure as a vector PDF in results/figures/.

  Fig. 1  Model structure: state-transition graph of Q_base and the
          data-generating pipeline.
  Fig. 2  Latent paths versus observed data for one patient per latent class.
  Fig. 3  Cohort-level realism of the reference cohort (N = 500, seed 42).
  Fig. 4  Kaplan-Meier estimates: SPMS conversion and disability milestones.
  Fig. 5  Known-truth benchmark of generator-matrix estimators
          (requires results/benchmark_estimators.csv).

Definitions are shared with analyze_manuscript_stats.py via ms_synth.summaries.
Usage: python scripts/make_figures.py [--only 1 2 ...]
"""
from __future__ import annotations

import argparse
import logging
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from lifelines import KaplanMeierFitter  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

from ms_synth import summaries as S  # noqa: E402
from ms_synth.synthetic_data import (  # noqa: E402
    _parse_state_defs, build_Q_from_config, generate_synthetic_ms_dataset)
from ms_synth.utils import load_config  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results" / "figures"
CFG = load_config(str(ROOT / "configs" / "synthetic.yaml"))
SEED, N = 42, 500

# Okabe-Ito colour-blind-safe palette
C = dict(blue="#0072B2", orange="#E69F00", green="#009E73", red="#D55E00",
         purple="#CC79A7", sky="#56B4E9", yellow="#F0E442", grey="#7F7F7F", black="#000000")
CLASS_COL = {"stable": C["green"], "moderate": C["orange"], "aggressive": C["red"]}
EST_COL = {"oracle": C["grey"], "crude": C["red"], "panel": C["blue"]}
EST_LAB = {"oracle": "Oracle (latent paths)", "crude": "Crude (visit counts)",
           "panel": "Panel likelihood"}

plt.rcParams.update({
    "font.size": 8, "axes.titlesize": 8.5, "axes.labelsize": 8, "legend.fontsize": 7,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.spines.top": False,
    "axes.spines.right": False, "pdf.fonttype": 42, "font.family": "DejaVu Sans",
})


def panel_label(ax, s, x=-0.14, y=1.04):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom")


def save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / name, bbox_inches="tight")
    plt.close(fig)
    print("wrote", OUT / name)


def reference_cohort():
    logging.disable(logging.CRITICAL)
    df, lat = generate_synthetic_ms_dataset(cfg=CFG, n_patients=N, seed=SEED, return_latent=True)
    return df, lat, S.patient_table(df)


# --------------------------------------------------------------------------- Fig 1
def edge_type(i, j, sd):
    ei, ej = sd[i]["edss_centre"], sd[j]["edss_centre"]
    if sd[j]["phase"] == "ABS":
        return "absorption"
    if sd[i]["phase"] == "RRMS" and sd[j]["phase"] == "SPMS":
        return "conversion"
    if ej > ei + 0.5:
        return "progression"
    if ej < ei - 0.5:
        return "recovery"
    return "relapse"


def fig1():
    Q = build_Q_from_config(CFG)
    sd = _parse_state_defs(CFG)
    pos = {0: (0, 0), 1: (0, 1), 2: (1, 0), 3: (1, 1), 4: (2, 0), 5: (2, 1),
           7: (3.2, 0), 6: (3.2, 1), 9: (4.2, 0), 8: (4.2, 1), 10: (5.2, 0.5), 11: (6.2, 0.5)}
    short = {0: "0\nmild", 1: "1\nmild*", 2: "2\nmod", 3: "3\nmod*", 4: "4\nadv", 5: "5\nadv*",
             6: "6\nearly*", 7: "7\nearly", 8: "8\nlate*", 9: "9\nlate", 10: "10\nsevere",
             11: "11\nabs."}
    tcol = {"relapse": C["blue"], "progression": C["orange"], "conversion": C["red"],
            "recovery": C["green"], "absorption": C["grey"]}
    tlab = {"relapse": "Relapse onset / remission", "progression": "Progression (incl. PIRA)",
            "conversion": "Conversion RRMS \u2192 SPMS", "recovery": "Recovery",
            "absorption": "Absorption"}

    fig = plt.figure(figsize=(7.2, 4.9))
    ax = fig.add_axes([0.02, 0.40, 0.96, 0.56])
    ax.axvspan(-0.45, 2.45, color=C["sky"], alpha=0.10, lw=0)
    ax.axvspan(2.75, 5.65, color=C["red"], alpha=0.07, lw=0)
    ax.text(1.0, 1.55, "RRMS phase", ha="center", fontsize=8, color=C["blue"], fontweight="bold")
    ax.text(4.2, 1.55, "SPMS phase", ha="center", fontsize=8, color=C["red"], fontweight="bold")
    for lab, x in (("EDSS 0\u20131.5", 0), ("2\u20133.5", 1), ("4\u20135.5", 2),
                   ("6\u20137", 3.2), ("7.5\u20138.5", 4.2), ("9\u20139.5", 5.2), ("9.5", 6.2)):
        ax.text(x, -0.55, lab, ha="center", fontsize=7, color=C["grey"])
    ax.text(-0.75, 1.0, "relapse-\nactive", ha="center", va="center", fontsize=7, color=C["grey"])
    ax.text(-0.75, 0.0, "inactive", ha="center", va="center", fontsize=7, color=C["grey"])
    qmax = Q[Q > 0].max()
    used = set()
    for i in range(12):
        for j in range(12):
            if i == j or Q[i, j] <= 0:
                continue
            t = edge_type(i, j, sd); used.add(t)
            (x0, y0), (x1, y1) = pos[i], pos[j]
            rad = 0.25 if (j, i) in zip(*np.nonzero(Q > 0)) else 0.08
            lw = 0.6 + 3.2 * Q[i, j] / qmax
            ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), connectionstyle=f"arc3,rad={rad}",
                                         arrowstyle="-|>", mutation_scale=8, lw=lw,
                                         color=tcol[t], alpha=0.85, shrinkA=13, shrinkB=13))
    for s, (x, y) in pos.items():
        fc = "white" if sd[s]["phase"] != "ABS" else "#EEEEEE"
        ax.add_patch(plt.Circle((x, y), 0.2, fc=fc, ec=C["black"], lw=0.8, zorder=3))
        ax.text(x, y, short[s], ha="center", va="center", fontsize=6.2, zorder=4, linespacing=0.9)
    handles = [plt.Line2D([], [], color=tcol[t], lw=2, label=tlab[t]) for t in tlab if t in used]
    handles.append(plt.Line2D([], [], color="k", lw=0.6, label="width \u221d intensity"))
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=3,
              frameon=False, fontsize=7)
    ax.set_xlim(-1.1, 6.6); ax.set_ylim(-0.75, 1.75); ax.set_aspect("equal"); ax.axis("off")
    panel_label(ax, "A", x=0.0, y=0.97)

    # (B) pipeline
    bx = fig.add_axes([0.02, 0.0, 0.96, 0.25]); bx.axis("off"); bx.set_xlim(0, 10); bx.set_ylim(0, 2.2)
    boxes = [("Covariates\nclass, onset age,\nsex, DMT", 0.05),
             ("Patient generator\n$\\mathbf{Q}^{(p)}$ = scaled\n$\\mathbf{Q}^{\\mathrm{base}}$", 2.05),
             ("Exact CTMC path\n(Gillespie)", 4.05),
             ("Observation process\nirregular & post-relapse\nvisits, dropout", 6.05),
             ("Observed cohort\nEDSS (noisy, 0.5 grid),\nrelapses, covariates", 8.05)]
    for k, (txt, x) in enumerate(boxes):
        fc = "#F2F2F2" if k < 4 else "#DCEBF7"
        bx.add_patch(FancyBboxPatch((x, 0.95), 1.75, 1.1, boxstyle="round,pad=0.04",
                                    fc=fc, ec="k", lw=0.7))
        bx.text(x + 0.875, 1.5, txt, ha="center", va="center", fontsize=6.8)
        if k < 4:
            bx.annotate("", xy=(x + 2.0, 1.5), xytext=(x + 1.8, 1.5),
                        arrowprops=dict(arrowstyle="-|>", lw=0.9))
    bx.add_patch(FancyBboxPatch((2.05, 0.05), 5.75, 0.55, boxstyle="round,pad=0.04",
                                fc="#FFF4D6", ec=C["orange"], lw=0.7, ls="--"))
    bx.text(4.925, 0.33, "Known truth retained (return_latent=True): $\\mathbf{Q}^{(p)}$, "
            "latent path, latent class, follow-up", ha="center", va="center", fontsize=6.8)
    for x in (2.925, 4.925, 6.925):
        bx.annotate("", xy=(x, 0.62), xytext=(x, 0.93),
                    arrowprops=dict(arrowstyle="-|>", lw=0.7, color=C["orange"], ls="--"))
    panel_label(bx, "B", x=0.0, y=0.95)
    save(fig, "fig1_model_structure.pdf")


# --------------------------------------------------------------------------- Fig 2
def fig2(df, lat, pat):
    sd = _parse_state_defs(CFG)
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.4), sharey=True)
    for ax, cls, lab in zip(axes, ("stable", "moderate", "aggressive"), "ABC"):
        sub = pat[(pat.cls == cls) & (pat.fu >= 15)]
        target = sub.edss_last.median()
        pid = (sub.edss_last - target).abs().sort_values(kind="stable").index[0]
        ev, T = lat[pid]["events"], lat[pid]["follow_up"]
        tt = [t for t, _ in ev if t <= T] + [T]
        yy = [sd[s]["edss_centre"] for t, s in ev if t <= T]
        yy = yy + [yy[-1]]
        ax.step(tt, yy, where="post", color=CLASS_COL[cls], lw=1.4, label="latent state (band centre)")
        conv = pat.at[pid, "t_spms"]
        if np.isfinite(conv):
            ax.axvspan(conv, T, color=C["red"], alpha=0.07, lw=0)
            ax.text(min(conv + 0.2, 17.6), 9.1, "SPMS", fontsize=6.5, color=C["red"])
        g = df[df.patient_id == pid]
        ax.plot(g.disease_duration_yr, g.EDSS, "o", ms=2.6, mfc="white", mec="k", mew=0.6,
                label="observed EDSS at visit")
        r = g[g.relapse == 1]
        ax.plot(r.disease_duration_yr, np.full(len(r), -0.45), "^", ms=3.2, color=C["purple"],
                label="relapse recorded")
        ax.set_title(f"{cls.capitalize()} class", color=CLASS_COL[cls])
        ax.set_xlabel("Years since onset"); ax.set_xlim(0, 20.3); ax.set_ylim(-0.8, 9.8)
        panel_label(ax, lab, x=-0.08 if lab != "A" else -0.2)
    axes[0].set_ylabel("EDSS")
    axes[0].legend(loc="upper left", frameon=False, fontsize=6.3)
    fig.tight_layout(w_pad=0.6)
    save(fig, "fig2_latent_vs_observed.pdf")


# --------------------------------------------------------------------------- Fig 3
def fig3(df, pat):
    fig, ax = plt.subplots(2, 2, figsize=(7.2, 4.8))
    a = ax[0, 0]
    bins = np.arange(14, 62, 2.5)
    a.hist([pat.loc[pat.sex == "F", "age_at_onset"], pat.loc[pat.sex == "M", "age_at_onset"]],
           bins=bins, stacked=True, color=[C["purple"], C["sky"]], label=["Female", "Male"],
           edgecolor="white", lw=0.3)
    a.axvline(pat.age_at_onset.mean(), color="k", ls="--", lw=0.8)
    a.text(pat.age_at_onset.mean() + 0.8, a.get_ylim()[1] * 0.9,
           f"mean {pat.age_at_onset.mean():.1f} yr", fontsize=6.5)
    a.set_xlabel("Age at onset (years)"); a.set_ylabel("Patients"); a.legend(frameon=False)
    panel_label(a, "A")

    b = ax[0, 1]
    e = np.arange(-0.25, 10.0, 0.5)
    b.hist(pat.edss_base, bins=e, alpha=0.45, color=C["green"], label="Baseline visit")
    b.hist(pat.edss_last, bins=e, histtype="step", lw=1.4, color=C["red"], label="Last visit")
    b.set_xlabel("Observed EDSS"); b.set_ylabel("Patients"); b.legend(frameon=False)
    panel_label(b, "B")

    c = ax[1, 0]
    d = df.sort_values(["patient_id", "disease_duration_yr"]).copy()
    d["gap_m"] = d.groupby("patient_id")["disease_duration_yr"].diff() * 12
    d["phase_prev"] = d.groupby("patient_id")["phase"].shift()
    bins = np.arange(0, 30.5, 1)
    for ph, col in (("RRMS", C["blue"]), ("SPMS", C["red"])):
        v = d.loc[d.phase_prev == ph, "gap_m"].dropna()
        c.hist(v, bins=bins, density=True, histtype="step", lw=1.3, color=col,
               label=f"{ph} (median {v.median():.1f} mo)")
    c.set_xlabel("Inter-visit interval (months)"); c.set_ylabel("Density"); c.legend(frameon=False)
    panel_label(c, "C")

    dd = ax[1, 1]
    rr, sp = S.pooled_arr_by_phase(df, pat)
    arr = S.pooled_arr_by_group(pat, "dmt")
    labels = ["RRMS\nphase", "SPMS\nphase", "No\nDMT", "Moderate\nDMT", "High\nDMT"]
    vals = [rr, sp, arr["none"], arr["moderate_dmt"], arr["high_dmt"]]
    cols = [C["blue"], C["red"], C["grey"], C["orange"], C["purple"]]
    dd.bar(range(5), vals, color=cols, width=0.65)
    dd.errorbar([1.42], [0.32], yerr=[[0.09], [0.09]], fmt="none", ecolor="k", capsize=3, lw=0.9)
    dd.text(1.42, 0.425, "published\n0.23\u20130.41", fontsize=5.8, ha="center")
    for k, v in enumerate(vals):
        dd.text(k, v + 0.012, f"{v:.2f}", ha="center", fontsize=6.5)
    for k, key in ((3, "moderate_dmt"), (4, "high_dmt")):
        dd.text(k, 0.03, f"\u2212{100 * (1 - arr[key] / arr['none']):.0f}%", ha="center",
                fontsize=6.5, color="white", fontweight="bold")
    dd.axvline(1.5, color=C["grey"], lw=0.5, ls=":")
    dd.set_xticks(range(5)); dd.set_xticklabels(labels)
    dd.set_ylabel("Relapses per patient-year"); dd.set_ylim(0, 0.62)
    panel_label(dd, "D")
    fig.tight_layout(h_pad=1.2, w_pad=1.5)
    save(fig, "fig3_cohort_realism.pdf")


# --------------------------------------------------------------------------- Fig 4
def _km(ax, t, e, color, label, ls="-", lw=1.6, band=False):
    k = KaplanMeierFitter().fit(t, e)
    x = k.survival_function_.index.values
    y = 100 * (1 - k.survival_function_.values.ravel())
    ax.step(x, y, where="post", color=color, lw=lw, ls=ls, label=label)
    if band:
        ci = k.confidence_interval_survival_function_
        ax.fill_between(ci.index.values, 100 * (1 - ci.iloc[:, 1].values),
                        100 * (1 - ci.iloc[:, 0].values), step="post", color=color, alpha=0.15, lw=0)


def fig4(pat):
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.9))
    a = ax[0]
    _km(a, pat.spms_time, pat.spms_event, C["blue"], "All patients", lw=2.0, band=True)
    for cls, col in CLASS_COL.items():
        m = pat.cls == cls
        _km(a, pat.loc[m, "spms_time"], pat.loc[m, "spms_event"], col, cls.capitalize(), ls="--")
    a.set_xlabel("Years since onset"); a.set_ylabel("Cumulative SPMS conversion (%)")
    a.set_xlim(0, 20); a.set_ylim(0, 80); a.legend(frameon=False, loc="upper left")
    ticks = [0, 5, 10, 15, 20]
    a.set_xticks(ticks)
    a.set_xticklabels([f"{t}\n{int((pat.spms_time >= t).sum())}" for t in ticks])
    a.text(-0.5, -13.6, "at risk", fontsize=6, ha="right", color=C["grey"])
    panel_label(a, "A")
    b = ax[1]
    for k, col in ((3, C["sky"]), (6, C["red"]), (8, C["purple"])):
        t, e = S.milestone(pat, f"t_e{k}")
        _km(b, t, e, col, f"EDSS \u2265 {k}")
    for x, lab in ((10, "EDSS 3\nLondon Ont./Rennes"), (18, "EDSS 6\nLondon Ont.")):
        b.axvline(x, color=C["grey"], lw=0.6, ls=":")
        b.text(x + 0.2, 74, lab, fontsize=5.8, color=C["grey"], va="top")
    b.set_xlabel("Years since onset"); b.set_ylabel("Cumulative incidence (%)")
    b.set_xlim(0, 20); b.set_ylim(0, 80); b.legend(frameon=False, loc="upper left")
    panel_label(b, "B")
    fig.tight_layout(w_pad=2)
    save(fig, "fig4_kaplan_meier.pdf")


# --------------------------------------------------------------------------- Fig 5
def _logx(a):
    from matplotlib.ticker import FixedLocator, NullFormatter, NullLocator
    a.set_xscale("log")
    a.xaxis.set_major_locator(FixedLocator([250, 500, 1000, 2000, 4000]))
    a.xaxis.set_minor_locator(NullLocator()); a.xaxis.set_minor_formatter(NullFormatter())
    a.set_xticklabels(["250", "500", "1k", "2k", "4k"]); a.set_xlim(200, 5000)


def fig5():
    f = ROOT / "results" / "benchmark_summary.csv"
    fp = ROOT / "results" / "benchmark_estimators.csv"
    if not f.exists() or not fp.exists():
        print("skip Fig 5: run scripts/benchmark_estimators.py first"); return
    import json
    sm = pd.read_csv(f); raw = pd.read_csv(fp)
    fig = plt.figure(figsize=(7.2, 5.6))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1], hspace=0.75, wspace=0.42,
                          left=0.08, right=0.98, top=0.95, bottom=0.08)
    ax = [[fig.add_subplot(gs[r, c]) for c in range(3)] for r in range(2)]

    def curve(a, q, metric, est, ls, mk="o", lab=None, scale=1.0):
        d = q[q.estimator == est].sort_values("n_patients")
        a.errorbar(d.n_patients, scale * d[f"{metric}_mean"], yerr=scale * 1.96 * d[f"{metric}_mcse"],
                   color=EST_COL[est], ls=ls, marker=mk, ms=3, lw=1.2, capsize=1.5, label=lab)

    # A: relative bias vs true intensity at N = 4000 (homogeneous)
    a = ax[0][0]
    big = raw[(raw.n_patients == 4000) & (raw.population == "homogeneous")]
    big = big.assign(rb=100 * (big.q_hat - big.q_true) / big.q_true)
    pb = big.groupby(["visits", "estimator", "q_true", "from_state", "to_state"]).rb.mean().reset_index()
    for (vis, est, mk) in (("non-informative", "crude", "v"), ("non-informative", "panel", "o"),
                           ("informative", "panel", "o")):
        d = pb[(pb.visits == vis) & (pb.estimator == est)]
        a.semilogx(d.q_true, d.rb, mk, ms=3.2, color=EST_COL[est],
                   mfc=EST_COL[est] if vis == "non-informative" else "white", mew=0.9)
    a.axhline(0, color="k", lw=0.6)
    from matplotlib.ticker import FixedLocator, NullFormatter, NullLocator
    a.xaxis.set_major_locator(FixedLocator([0.05, 0.1, 0.2, 0.5]))
    a.xaxis.set_minor_locator(NullLocator()); a.xaxis.set_minor_formatter(NullFormatter())
    a.set_xticklabels(["0.05", "0.1", "0.2", "0.5"])
    a.set_xlabel("True intensity (yr$^{-1}$)"); a.set_ylabel("Relative bias (%)")
    a.set_title("Per-transition bias, N = 4000", fontsize=7)
    panel_label(a, "A", x=-0.3)

    hom = sm[(sm.population == "homogeneous") & (sm.interval_mult == 1.0)]
    for a, metric, ylab, lab, sc in ((ax[0][1], "abs_err", "Mean |relative error|", "B", 1.0),
                                     (ax[0][2], "bias", "Mean relative bias (%)", "C", 100.0)):
        for est in ("crude", "panel"):
            curve(a, hom[hom.visits == "informative"], metric, est, "-", "o", scale=sc)
            curve(a, hom[hom.visits == "non-informative"], metric, est, "--", "s", scale=sc)
        curve(a, hom[hom.visits == "non-informative"], metric, "oracle", ":", "D", scale=sc)
        _logx(a); a.set_xlabel("Cohort size N"); a.set_ylabel(ylab)
        if metric == "bias":
            a.axhline(0, color="k", lw=0.6)
        panel_label(a, lab, x=-0.3)
    ax[0][1].set_title("Homogeneous population", fontsize=7)
    ax[0][2].set_title("Homogeneous population", fontsize=7)
    H = [plt.Line2D([], [], color=EST_COL["crude"], marker="v", ls="-", ms=4, label="Crude (visit counts)"),
         plt.Line2D([], [], color=EST_COL["panel"], marker="o", ls="-", ms=4, label="Panel likelihood"),
         plt.Line2D([], [], color=EST_COL["oracle"], marker="D", ls=":", ms=3, label="Oracle (latent paths)"),
         plt.Line2D([], [], color="k", ls="-", label="outcome-dependent visits (solid / open)"),
         plt.Line2D([], [], color="k", ls="--", label="non-informative visits (dashed / filled)")]
    fig.legend(handles=H, loc="upper center", bbox_to_anchor=(0.53, 0.535), ncol=3, frameon=False, fontsize=6.3)

    d = ax[1][0]
    inf = sm[(sm.visits == "informative") & (sm.interval_mult == 1.0)]
    for pop, ls, mk in (("homogeneous", "-", "o"), ("heterogeneous", ":", "^")):
        for est in ("oracle", "panel", "crude"):
            curve(d, inf[inf.population == pop], "abs_err", est, ls, mk)
    _logx(d); d.set_xlabel("Cohort size N"); d.set_ylabel("Mean |relative error|")
    d.set_title("Homogeneous (solid) vs\nheterogeneous (dotted)", fontsize=7)
    panel_label(d, "D", x=-0.3)

    e = ax[1][1]
    vis = sm[(sm.population == "homogeneous") & (sm.visits == "non-informative") & (sm.n_patients == 1000)]
    for est, mk in (("crude", "v"), ("panel", "o"), ("oracle", "D")):
        q = vis[vis.estimator == est].sort_values("interval_mult")
        e.errorbar(7.0 * q.interval_mult, 100 * q.bias_mean, yerr=196 * q.bias_mcse, color=EST_COL[est],
                   marker=mk, ms=3, lw=1.2, capsize=1.5, ls="--" if est != "oracle" else ":")
    e.axhline(0, color="k", lw=0.6)
    e.set_xlabel("Mean RRMS visit interval (months)"); e.set_ylabel("Mean relative bias (%)")
    e.set_xticks([3.5, 7, 14]); e.set_title("Visit frequency (N = 1000)", fontsize=7)
    panel_label(e, "E", x=-0.3)

    ff = ax[1][2]
    R = json.loads((ROOT / "results" / "results.json").read_text())
    tab = pd.DataFrame(R["Q_table"])
    for est, mk in (("oracle", "D"), ("crude", "v"), ("panel", "o")):
        ff.loglog(tab.Q_true, tab[f"Q_{est}"].clip(lower=1e-3), mk, ms=3, color=EST_COL[est],
                  mfc="none" if est == "oracle" else EST_COL[est],
                  label=f"{est}: {R['Q_estimators'][est]['mean_abs_relerr']:.2f}")
    lim = [0.02, 1.5]; ff.plot(lim, lim, "k-", lw=0.6)
    ff.set_xlim(lim); ff.set_ylim(lim)
    ff.set_xlabel("True intensity (yr$^{-1}$)"); ff.set_ylabel("Estimate (yr$^{-1}$)")
    ff.set_title("Reference cohort (N = 500)", fontsize=7)
    ff.legend(frameon=False, fontsize=5.8, loc="upper left", title="mean |rel. err.|", title_fontsize=5.8)
    panel_label(ff, "F", x=-0.3)
    save(fig, "fig5_estimator_benchmark.pdf")


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--only", nargs="*", type=int)
    a = ap.parse_args(); want = set(a.only or [1, 2, 3, 4, 5])
    if want & {2, 3, 4}:
        df, lat, pat = reference_cohort()
    if 1 in want: fig1()
    if 2 in want: fig2(df, lat, pat)
    if 3 in want: fig3(df, pat)
    if 4 in want: fig4(pat)
    if 5 in want: fig5()


if __name__ == "__main__":
    main()
