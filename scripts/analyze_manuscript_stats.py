"""Generate the seed-42 N=500 cohort and compute every statistic the
manuscript reports, plus the new time-to-milestone validation. Writes
results.json and prints a readable summary."""
import sys, json, logging
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
logging.disable(logging.CRITICAL)
import numpy as np, pandas as pd
from src.utils import load_config
from src.synthetic_data import (
    generate_synthetic_ms_dataset, build_Q_from_config, compare_empirical_to_Q,
)
from lifelines import KaplanMeierFitter, CoxPHFitter

SEED, N = 42, 500
cfg = load_config(str(Path(__file__).resolve().parent.parent / 'configs' / 'synthetic.yaml'))
df = generate_synthetic_ms_dataset(cfg=cfg, n_patients=N, seed=SEED)

R = {}  # results

# ---- patient-level frame ----
g = df.sort_values(['patient_id','disease_duration_yr']).groupby('patient_id', sort=False)
pat = pd.DataFrame({
    'age_at_onset': g['age_at_onset'].first(),
    'sex': g['sex'].first(),
    'dmt': g['DMT_status'].first(),
    'cls': g['latent_class'].first(),
    'fu': g['disease_duration_yr'].max(),
    'n_visits': g.size(),
    'edss_base': g['EDSS'].first(),
    'edss_last': g['EDSS'].last(),
    'n_relapse': g['relapse'].sum(),
    'ever_spms': g['phase'].apply(lambda s: (s=='SPMS').any()),
})
# time to SPMS (first SPMS visit)
def first_time(gg, mask):
    sub = gg[mask(gg)]
    return sub['disease_duration_yr'].min() if len(sub) else np.nan
spms_t, e6_t, e3_t, e8_t = {}, {}, {}, {}
for pid, gg in g:
    spms = gg[gg['phase']=='SPMS']['disease_duration_yr']
    spms_t[pid] = spms.min() if len(spms) else np.nan
    for thr, store in [(3,e3_t),(6,e6_t),(8,e8_t)]:
        hit = gg[gg['EDSS_true']>=thr]['disease_duration_yr']
        store[pid] = hit.min() if len(hit) else np.nan
pat['t_spms']=pd.Series(spms_t); pat['t_e3']=pd.Series(e3_t)
pat['t_e6']=pd.Series(e6_t); pat['t_e8']=pd.Series(e8_t)

# EDSS slope per patient (OLS on observed EDSS vs duration)
def slope(gg):
    x=gg['disease_duration_yr'].to_numpy(); y=gg['EDSS'].to_numpy()
    if len(x)<2 or x.std()==0: return 0.0
    return float(np.polyfit(x,y,1)[0])
pat['slope']=g.apply(slope)

# ---- Demographics ----
R['n_patients']=int(N)
R['total_visits']=int(len(df))
R['visits_per_patient_mean']=round(float(pat['n_visits'].mean()),1)
R['visits_per_patient_median']=round(float(pat['n_visits'].median()),1)
R['female_pct']=round(100*float((pat['sex']=='F').mean()),1)
R['onset_age_mean']=round(float(pat['age_at_onset'].mean()),1)
R['onset_age_sd']=round(float(pat['age_at_onset'].std()),1)
# ---- Follow-up ----
R['fu_median']=round(float(pat['fu'].median()),1)
R['fu_mean']=round(float(pat['fu'].mean()),1)
R['fu_min']=round(float(pat['fu'].min()),1)
R['fu_max']=round(float(pat['fu'].max()),1)
# ---- Visit intervals (months) ----
iv=[]
for pid,gg in g:
    t=np.sort(gg['disease_duration_yr'].to_numpy())
    iv.extend(np.diff(t)*12.0)
iv=np.array(iv)
R['interval_mean_m']=round(float(iv.mean()),2)
R['interval_median_m']=round(float(np.median(iv)),2)
R['interval_sd_m']=round(float(iv.std(ddof=1)),2)
# ---- EDSS ----
R['edss_base_mean']=round(float(pat['edss_base'].mean()),2)
R['edss_base_sd']=round(float(pat['edss_base'].std()),2)
R['edss_last_mean']=round(float(pat['edss_last'].mean()),2)
R['slope_mean']=round(float(pat['slope'].mean()),3)
R['slope_median']=round(float(pat['slope'].median()),3)
R['pos_slope_pct']=round(100*float((pat['slope']>0).mean()),1)
R['fast_prog_pct']=round(100*float((pat['slope']>0.5).mean()),1)

# ---- SPMS conversion overall + by class (10/20 yr cumulative incidence via KM) ----
def cuminc(times, events, at):
    kmf=KaplanMeierFitter().fit(times, events)
    out={}
    for t in at:
        try: out[t]=round(100*float(1-kmf.survival_function_at_times(t).iloc[0]),1)
        except Exception: out[t]=None
    # median event time (where cuminc crosses 50%)
    med=kmf.median_survival_time_
    out['median']=(None if not np.isfinite(med) else round(float(med),1))
    return out

pat['spms_event']=pat['t_spms'].notna().astype(int)
pat['spms_time']=pat['t_spms'].where(pat['spms_event']==1, pat['fu'])
R['spms_overall']=cuminc(pat['spms_time'],pat['spms_event'],[5,10,15,20])
R['spms_overall_pct']=round(100*float(pat['ever_spms'].mean()),1)

R['by_class']={}
for c in ['stable','moderate','aggressive']:
    m=pat['cls']==c
    ci=cuminc(pat.loc[m,'spms_time'],pat.loc[m,'spms_event'],[10,20])
    R['by_class'][c]={
        'n':int(m.sum()),
        'pct':round(100*float(m.mean()),1),
        'spms10':ci[10],'spms20':ci[20],
        'edss_last':round(float(pat.loc[m,'edss_last'].mean()),2),
    }

# ---- Milestone KM medians (time from onset) ----
R['milestones']={}
for name,col in [('EDSS3','t_e3'),('EDSS6','t_e6'),('EDSS8','t_e8'),('SPMS','t_spms')]:
    ev=pat[col].notna().astype(int); tt=pat[col].where(ev==1,pat['fu'])
    ci=cuminc(tt,ev,[10,20])
    R['milestones'][name]={'median':ci['median'],'cum10':ci[10],'cum20':ci[20]}

# ---- Relapse rates by phase (pooled person-time) ----
def phase_arr():
    # RRMS person-time = time to SPMS (or last visit); SPMS = last - conversion
    rr_time=rr_rel=sp_time=sp_rel=0.0
    for pid,gg in g:
        conv=spms_t[pid]; last=gg['disease_duration_yr'].max()
        rrms_end = conv if not np.isnan(conv) else last
        rr_time += rrms_end
        rr_rel += gg[gg['phase']=='RRMS']['relapse'].sum()
        if not np.isnan(conv):
            sp_time += max(last-conv,0.0)
            sp_rel += gg[gg['phase']=='SPMS']['relapse'].sum()
    return rr_rel/max(rr_time,1e-9), sp_rel/max(sp_time,1e-9)
arr_rrms,arr_spms=phase_arr()
R['arr_rrms']=round(float(arr_rrms),3)
R['arr_spms']=round(float(arr_spms),3)

# ---- DMT ARR (pooled person-time, all phases) ----
R['dmt']={}
base_arr=None
for d in ['none','moderate_dmt','high_dmt']:
    m=pat['dmt']==d
    rel=pat.loc[m,'n_relapse'].sum(); ptime=pat.loc[m,'fu'].sum()
    arr=rel/max(ptime,1e-9)
    if d=='none': base_arr=arr
    R['dmt'][d]={'n':int(m.sum()),'pct':round(100*float(m.mean()),1),
                 'arr':round(float(arr),3),
                 'rel_reduction':(None if d=='none' else round(100*(1-arr/base_arr),1))}

# ---- Empirical Q recovery ----
Q_base=build_Q_from_config(cfg)
comp=compare_empirical_to_Q(df,Q_base)
R['Q_mean_relerr']=round(float(comp['relative_error'].mean()),3)
R['Q_median_relerr']=round(float(comp['relative_error'].median()),3)
R['Q_n_pairs']=int(len(comp))
R['Q_table']=comp.sort_values('relative_error').to_dict(orient='records')

# ---- Cox models: time to SPMS ----
pat['male']=(pat['sex']=='M').astype(int)
# early relapse in first 2 years
er_cnt={}
for pid,gg in g:
    er_cnt[pid]=int(gg[gg['disease_duration_yr']<=2.0]['relapse'].sum())
pat['early_relapse_cnt']=pd.Series(er_cnt)
pat['early_relapse_bin']=(pat['early_relapse_cnt']>=1).astype(int)
pat['age_z']=(pat['age_at_onset']-pat['age_at_onset'].mean())/pat['age_at_onset'].std()
pat['edssbase_z']=(pat['edss_base']-pat['edss_base'].mean())/pat['edss_base'].std()

def cox(cols):
    d=pat[cols+['spms_time','spms_event']].copy()
    cph=CoxPHFitter().fit(d,'spms_time','spms_event')
    s=cph.summary
    return {c:{'HR':round(float(np.exp(s.loc[c,'coef'])),2),
               'lo':round(float(np.exp(s.loc[c,'coef lower 95%'])),2),
               'hi':round(float(np.exp(s.loc[c,'coef upper 95%'])),2),
               'p':round(float(s.loc[c,'p']),4)} for c in cols}
R['cox']={}
try: R['cox']['early_bin']=cox(['early_relapse_bin'])
except Exception as e: R['cox']['early_bin']=str(e)
try: R['cox']['early_cnt']=cox(['early_relapse_cnt'])
except Exception as e: R['cox']['early_cnt']=str(e)
try: R['cox']['adjusted']=cox(['age_z','male','edssbase_z'])
except Exception as e: R['cox']['adjusted']=str(e)

out = Path(__file__).resolve().parent.parent / 'results' / 'results.json'
out.parent.mkdir(exist_ok=True)
json.dump(R, open(out,'w'), indent=2)

# ---- readable summary ----
p=lambda *a: print(*a)
p("="*64); p(f"SEED-{SEED} COHORT, N={N}"); p("="*64)
p(f"Total visits: {R['total_visits']}  | visits/pt mean {R['visits_per_patient_mean']} median {R['visits_per_patient_median']}")
p(f"Female: {R['female_pct']}%  | onset age {R['onset_age_mean']}±{R['onset_age_sd']}")
p(f"Follow-up median {R['fu_median']} mean {R['fu_mean']} range {R['fu_min']}-{R['fu_max']}")
p(f"Interval months: mean {R['interval_mean_m']} median {R['interval_median_m']} sd {R['interval_sd_m']}")
p(f"EDSS base {R['edss_base_mean']}±{R['edss_base_sd']} last {R['edss_last_mean']} | slope mean {R['slope_mean']} median {R['slope_median']}")
p(f"  pos slope {R['pos_slope_pct']}% | fast(>0.5) {R['fast_prog_pct']}%")
p(f"SPMS overall {R['spms_overall_pct']}% | KM cuminc {R['spms_overall']}")
for c,v in R['by_class'].items():
    p(f"  {c:10s} n={v['n']} ({v['pct']}%) SPMS10={v['spms10']}% SPMS20={v['spms20']}% edss_last={v['edss_last']}")
p("Milestone KM medians (yr):")
for k,v in R['milestones'].items():
    p(f"  {k:6s} median={v['median']} cum10={v['cum10']}% cum20={v['cum20']}%")
p(f"ARR RRMS {R['arr_rrms']} | ARR SPMS {R['arr_spms']}")
for d,v in R['dmt'].items():
    p(f"  DMT {d:12s} n={v['n']} ({v['pct']}%) ARR={v['arr']} redux={v['rel_reduction']}")
p(f"Q recovery: mean relerr {R['Q_mean_relerr']} median {R['Q_median_relerr']} over {R['Q_n_pairs']} pairs")
p("Cox (time-to-SPMS):")
for k,v in R['cox'].items(): p(f"  {k}: {v}")

# ============================================================
# KM figure (single source of truth: same definitions as the
# tables above -> figure and tables cannot drift). Panel A: time
# to SPMS overall + by latent class. Panel B: time to first-reaching
# EDSS 3/6/8 milestones. Saved to manuscript/figs/fig6_KM_curves.pdf
# ============================================================
def _make_km_figure():
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available; skipping KM figure")
        return

    def km_curve(times, events):
        k = KaplanMeierFitter().fit(times, events)
        t = k.survival_function_.index.values
        cuminc = 100.0 * (1.0 - k.survival_function_.values.flatten())
        ci = k.confidence_interval_survival_function_
        lo = 100.0 * (1.0 - ci.iloc[:, 1].values)  # 1 - upper(S)
        hi = 100.0 * (1.0 - ci.iloc[:, 0].values)  # 1 - lower(S)
        return t, cuminc, ci.index.values, lo, hi

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

    # -- Panel A: SPMS conversion, all + by class --
    ax = axes[0]
    t, c, tci, lo, hi = km_curve(pat["spms_time"], pat["spms_event"])
    ax.step(t, c, where="post", color="#2563EB", lw=2.5, label="All patients")
    ax.fill_between(tci, lo, hi, step="post", alpha=0.18, color="#2563EB")
    class_colors = {"stable": "#16A34A", "moderate": "#F59E0B", "aggressive": "#DC2626"}
    for cls, col in class_colors.items():
        m = pat["cls"] == cls
        tc, cc, *_ = km_curve(pat.loc[m, "spms_time"], pat.loc[m, "spms_event"])
        ax.step(tc, cc, where="post", color=col, lw=1.8, ls="--", alpha=0.85,
                label=cls.capitalize())
    for y in (30, 50):
        ax.axhline(y, color="grey", ls=":", lw=1.0, alpha=0.6)
    ax.set_xlim(0, 20); ax.set_ylim(-1, 80)
    ax.set_xlabel("Disease duration (years)", fontsize=11)
    ax.set_ylabel("Cumulative SPMS conversion (%)", fontsize=11)
    ax.set_title("(A)  Time to SPMS conversion", fontweight="bold", fontsize=12)
    ax.legend(fontsize=9, loc="upper left"); ax.grid(True, alpha=0.25)
    for t_nr in (0, 5, 10, 15, 20):
        ax.text(t_nr, -6, str(int((pat["spms_time"] >= t_nr).sum())),
                ha="center", fontsize=7.5, color="#2563EB")
    ax.text(-1.4, -6, "N at risk:", ha="right", fontsize=7.5)

    # -- Panel B: EDSS milestones (first-reaching), all patients --
    ax2 = axes[1]
    for label, col, color in [("EDSS \u2265 3", "t_e3", "#0891B2"),
                              ("EDSS \u2265 6", "t_e6", "#DC2626"),
                              ("EDSS \u2265 8", "t_e8", "#7C3AED")]:
        ev = pat[col].notna().astype(int)
        tt = pat[col].where(ev == 1, pat["fu"])
        t, c, *_ = km_curve(tt, ev)
        ax2.step(t, c, where="post", color=color, lw=2.2, label=label)
    ax2.set_xlim(0, 20); ax2.set_ylim(-1, 80)
    ax2.set_xlabel("Disease duration (years)", fontsize=11)
    ax2.set_ylabel("Cumulative incidence (%)", fontsize=11)
    ax2.set_title("(B)  Time to disability milestones", fontweight="bold", fontsize=12)
    ax2.legend(fontsize=9, loc="upper left"); ax2.grid(True, alpha=0.25)
    ax2.annotate("EDSS\u22656 \u2261 SPMS\n(coincident by construction)",
                 xy=(15, 100 * (1 - KaplanMeierFitter().fit(
                     pat["t_e6"].where(pat["t_e6"].notna(), pat["fu"]),
                     pat["t_e6"].notna().astype(int)
                 ).survival_function_at_times(15).iloc[0])),
                 xytext=(9.0, 55), fontsize=8, color="#DC2626",
                 arrowprops=dict(arrowstyle="->", color="#DC2626", lw=1))

    plt.suptitle("Kaplan\u2013Meier estimates \u2014 synthetic MS cohort (N = 500, seed 42)",
                 fontsize=12, fontweight="bold", y=1.01)
    plt.tight_layout()
    figpath = Path(__file__).resolve().parent.parent / "manuscript" / "figs" / "fig6_KM_curves.pdf"
    figpath.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(figpath, bbox_inches="tight")
    plt.close()
    print("Figure saved:", figpath)

_make_km_figure()
