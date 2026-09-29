"""Tests for the known-truth interface, CT-MSM estimators and summaries."""
from __future__ import annotations

import copy
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.linalg import expm

from ms_synth import summaries as S
from ms_synth.msm_fit import (
    estimate_Q_crude, estimate_Q_from_paths, estimate_Q_panel, panel_transition_probs)
from ms_synth.synthetic_data import (
    build_Q_from_config, compare_empirical_to_Q, generate_synthetic_ms_dataset)
from ms_synth.utils import load_config

logging.disable(logging.CRITICAL)
CONFIG = Path(__file__).resolve().parent.parent / "configs" / "synthetic.yaml"


@pytest.fixture(scope="module")
def cfg():
    return load_config(str(CONFIG))


@pytest.fixture(scope="module")
def homogeneous(cfg):
    c = copy.deepcopy(cfg)
    c["cohort"]["class_fractions"] = {"stable": 0.0, "moderate": 1.0, "aggressive": 0.0}
    c["age_sex_effects"]["age_progression_multiplier_per_decade"] = 1.0
    c["age_sex_effects"]["male_progression_multiplier"] = 1.0
    c["demographics"]["dmt_status_probs"] = {"none": 1.0, "moderate_dmt": 0.0, "high_dmt": 0.0}
    df, lat = generate_synthetic_ms_dataset(cfg=c, n_patients=400, seed=7, return_latent=True)
    return c, df, lat


class TestKnownTruthInterface:
    def test_return_latent_does_not_change_observed_data(self, cfg):
        a = generate_synthetic_ms_dataset(cfg=cfg, n_patients=40, seed=3)
        b, lat = generate_synthetic_ms_dataset(cfg=cfg, n_patients=40, seed=3, return_latent=True)
        pd.testing.assert_frame_equal(a, b)
        assert set(lat) == set(b["patient_id"].unique())

    def test_latent_paths_consistent_with_visits(self, cfg):
        df, lat = generate_synthetic_ms_dataset(cfg=cfg, n_patients=30, seed=5, return_latent=True)
        for pid, g in df.groupby("patient_id"):
            ev = lat[pid]["events"]
            assert g["disease_duration_yr"].max() <= lat[pid]["follow_up"] + 1e-6
            for t, s in zip(g["disease_duration_yr"], g["state"]):
                cur = [st for (tt, st) in ev if tt <= t + 1e-9][-1]
                assert cur == s

    def test_patient_Q_rows_sum_to_zero(self, cfg):
        _, lat = generate_synthetic_ms_dataset(cfg=cfg, n_patients=20, seed=1, return_latent=True)
        for rec in lat.values():
            np.testing.assert_allclose(rec["Q_patient"].sum(axis=1), 0.0, atol=1e-12)


class TestEstimators:
    def test_vectorised_panel_probs_match_expm(self, cfg):
        Q = build_Q_from_config(cfg)
        rng = np.random.default_rng(0)
        a = rng.integers(0, 11, 200); b = rng.integers(0, 12, 200); dt = rng.uniform(0.01, 3, 200)
        ref = np.array([expm(Q * d)[i, j] for i, j, d in zip(a, b, dt)])
        np.testing.assert_allclose(panel_transition_probs(Q, a, b, dt), ref, atol=1e-12)

    def test_crude_matches_legacy_comparison(self, cfg):
        df = generate_synthetic_ms_dataset(cfg=cfg, n_patients=60, seed=11)
        Q = build_Q_from_config(cfg)
        comp = compare_empirical_to_Q(df, Q)
        Qc = estimate_Q_crude(df, 12)
        for _, r in comp.iterrows():
            assert Qc[int(r.from_state), int(r.to_state)] == pytest.approx(r.Q_empirical, abs=1e-5)

    def test_oracle_recovers_truth_in_homogeneous_cohort(self, homogeneous):
        c, _, lat = homogeneous
        Q = build_Q_from_config(c)
        Qo = estimate_Q_from_paths(lat, 12)
        rel = [abs(Qo[i, j] - Q[i, j]) / Q[i, j]
               for i in range(12) for j in range(12) if i != j and Q[i, j] > 0]
        assert np.median(rel) < 0.2

    def test_panel_mle_less_biased_than_crude(self, homogeneous):
        c, df, _ = homogeneous
        Q = build_Q_from_config(c)
        nz = [(i, j) for i in range(12) for j in range(12) if i != j and Q[i, j] > 0]
        bias = lambda Qh: np.mean([(Qh[i, j] - Q[i, j]) / Q[i, j] for i, j in nz])
        b_crude = bias(estimate_Q_crude(df, 12))
        b_panel = bias(estimate_Q_panel(df, 12, allowed=nz))
        assert b_crude < -0.05                 # visit-level counting under-estimates
        assert abs(b_panel) < abs(b_crude)


class TestSummaries:
    def test_flat_trajectory_not_counted_as_worsening(self):
        df = pd.DataFrame({
            "patient_id": ["p"] * 4, "disease_duration_yr": [0.0, 0.7, 1.3, 2.9],
            "EDSS": [2.0] * 4, "EDSS_true": [2.1] * 4, "phase": ["RRMS"] * 4,
            "relapse": [0] * 4, "age_at_onset": [30.0] * 4, "sex": ["F"] * 4,
            "DMT_status": ["none"] * 4, "latent_class": ["stable"] * 4})
        pat = S.patient_table(df)
        assert abs(pat["slope"].iloc[0]) < S.SLOPE_TOL
        assert S.positive_slope_fraction(pat) == 0.0

    def test_reference_cohort_headline_numbers(self, cfg):
        df = generate_synthetic_ms_dataset(cfg=cfg, n_patients=500, seed=42)
        pat = S.patient_table(df)
        assert len(df) == 14637
        assert round(100 * S.positive_slope_fraction(pat), 1) == 73.6
        rr, sp = S.pooled_arr_by_phase(df, pat)
        assert (round(rr, 3), round(sp, 3)) == (0.396, 0.281)
