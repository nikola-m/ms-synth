"""
tests/test_synthetic_data.py
============================
Unit and integration tests for src/synthetic_data.py.

Test categories
---------------
1. Q matrix construction  – shape, validity, sparse entries.
2. Patient scaling         – Q scaling by class / age / sex / DMT.
3. Gillespie simulator    – state validity, EDSS range, visit count.
4. Dataset generator       – schema, reproducibility, epidemiology.
5. Validation helpers      – empirical Q comparison, transition counts.
6. Edge cases              – single patient, absorbing states, tiny cohort.
"""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest


from ms_synth.synthetic_data import (
    build_Q_from_config,
    scale_Q_for_patient,
    simulate_patient_trajectory,
    generate_synthetic_ms_dataset,
    compare_empirical_to_Q,
    _parse_state_defs,
    _snap_to_grid,
    _state_at_time,
    _relapse_in_window,
    _DEFAULT_STATE_DEFS,
    _DEFAULT_Q_ENTRIES,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def minimal_cfg() -> dict:
    """Minimal valid config for fast tests."""
    return {
        "reproducibility": {"seed": 0},
        "cohort": {
            "n_patients": 20,
            "class_fractions": {"stable": 0.33, "moderate": 0.34, "aggressive": 0.33},
        },
        "states": {
            "n_states": 12,
            "onset_state_probs": [0.30, 0.35, 0.15, 0.10, 0.05, 0.05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            "definitions": None,  # use defaults
        },
        "Q_base": None,  # use built-in defaults
        "class_Q_factors": {
            "stable":     {"progression_scale": 0.35, "relapse_scale": 0.6, "recovery_scale": 1.4},
            "moderate":   {"progression_scale": 1.00, "relapse_scale": 1.0, "recovery_scale": 1.0},
            "aggressive": {"progression_scale": 2.20, "relapse_scale": 1.5, "recovery_scale": 0.55},
        },
        "age_sex_effects": {
            "age_progression_multiplier_per_decade": 1.25,
            "male_progression_multiplier": 1.15,
            "dmt_effects": {
                "none":         {"relapse_mult": 1.0, "progression_mult": 1.0},
                "moderate_dmt": {"relapse_mult": 0.65, "progression_mult": 0.85},
                "high_dmt":     {"relapse_mult": 0.40, "progression_mult": 0.70},
            },
        },
        "demographics": {
            "onset_age": {"distribution": "normal", "mean": 32.0, "std": 9.0, "min": 16.0, "max": 60.0},
            "sex_female_prob": 0.68,
            "dmt_status_probs": {"none": 0.2, "moderate_dmt": 0.35, "high_dmt": 0.45},
        },
        "visits": {
            "max_followup_years": 10.0,
            "min_followup_years": 2.0,
            "dropout_rate_per_year": 0.03,
            "mean_interval_rrms_months": 7.0,
            "mean_interval_spms_months": 9.0,
            "post_relapse_visit_months": 1.5,
            "jitter_months": 2.0,
        },
        "edss_noise": {"sigma": 0.30, "round_to": 0.5},
    }


@pytest.fixture
def Q_base(minimal_cfg) -> np.ndarray:
    return build_Q_from_config(minimal_cfg)


@pytest.fixture
def state_defs(minimal_cfg) -> dict:
    return _parse_state_defs(minimal_cfg)


@pytest.fixture
def small_df(minimal_cfg) -> pd.DataFrame:
    return generate_synthetic_ms_dataset(cfg=minimal_cfg)


# ---------------------------------------------------------------------------
# 1. Q matrix construction
# ---------------------------------------------------------------------------


class TestBuildQFromConfig:
    def test_shape(self, Q_base, minimal_cfg):
        n = minimal_cfg["states"]["n_states"]
        assert Q_base.shape == (n, n)

    def test_row_sums_near_zero(self, Q_base):
        np.testing.assert_allclose(Q_base.sum(axis=1), np.zeros(12), atol=1e-9)

    def test_diagonal_nonpositive(self, Q_base):
        assert np.all(np.diag(Q_base) <= 1e-9)

    def test_offdiag_nonneg(self, Q_base):
        n = Q_base.shape[0]
        for i in range(n):
            for j in range(n):
                if i != j:
                    assert Q_base[i, j] >= -1e-9

    def test_default_entries_present(self, Q_base):
        # At least one known entry: (0,1) = 0.60 from defaults
        assert abs(Q_base[0, 1] - 0.60) < 1e-9

    def test_custom_sparse_entries(self, minimal_cfg):
        cfg = dict(minimal_cfg)
        cfg["Q_base"] = [[0, 1, 0.9], [1, 0, 0.5]]
        cfg["states"] = dict(minimal_cfg["states"])
        cfg["states"]["n_states"] = 3
        Q = build_Q_from_config(cfg)
        assert abs(Q[0, 1] - 0.9) < 1e-9
        assert abs(Q[1, 0] - 0.5) < 1e-9

    def test_absorbing_state_zero_row(self, Q_base):
        # State 11 (absorbed) should have no out-transitions in defaults
        n = Q_base.shape[0]
        if n > 11:
            # diagonal only non-zero; off-diag should be ~0 for fully absorbing
            row_offdiag = Q_base[11].copy()
            row_offdiag[11] = 0.0
            # May not be fully 0 but should be very small (only entry is [10,11])
            assert row_offdiag.sum() < 1e-9


# ---------------------------------------------------------------------------
# 2. Patient Q scaling
# ---------------------------------------------------------------------------


class TestScaleQForPatient:
    def test_stable_slower_than_aggressive(self, Q_base, minimal_cfg, state_defs):
        Q_stable = scale_Q_for_patient(
            Q_base, "stable", 32.0, "F", "none", minimal_cfg
        )
        Q_aggr = scale_Q_for_patient(
            Q_base, "aggressive", 32.0, "F", "none", minimal_cfg
        )
        # Aggressive should have higher progression intensities overall
        # Compare a known progression edge: (2,4) index
        assert Q_aggr[2, 4] >= Q_stable[2, 4] - 1e-9

    def test_q_still_valid_after_scaling(self, Q_base, minimal_cfg):
        for cls in ["stable", "moderate", "aggressive"]:
            Q_pat = scale_Q_for_patient(Q_base, cls, 30.0, "F", "none", minimal_cfg)
            np.testing.assert_allclose(Q_pat.sum(axis=1), np.zeros(12), atol=1e-6)

    def test_older_onset_faster_progression(self, Q_base, minimal_cfg):
        Q_young = scale_Q_for_patient(Q_base, "moderate", 20.0, "F", "none", minimal_cfg)
        Q_old = scale_Q_for_patient(Q_base, "moderate", 50.0, "F", "none", minimal_cfg)
        # Sum of progression entries should be higher for older onset
        prog_young = sum(Q_young[i, j] for i in range(12) for j in range(i + 2, 12))
        prog_old = sum(Q_old[i, j] for i in range(12) for j in range(i + 2, 12))
        assert prog_old > prog_young

    def test_dmt_reduces_relapse(self, Q_base, minimal_cfg):
        Q_no_dmt = scale_Q_for_patient(Q_base, "moderate", 32.0, "F", "none", minimal_cfg)
        Q_high_dmt = scale_Q_for_patient(Q_base, "moderate", 32.0, "F", "high_dmt", minimal_cfg)
        # Relapse edge (0→1): high DMT should reduce it
        assert Q_high_dmt[0, 1] <= Q_no_dmt[0, 1] + 1e-9

    def test_male_faster_progression(self, Q_base, minimal_cfg):
        Q_female = scale_Q_for_patient(Q_base, "moderate", 32.0, "F", "none", minimal_cfg)
        Q_male = scale_Q_for_patient(Q_base, "moderate", 32.0, "M", "none", minimal_cfg)
        prog_f = sum(Q_female[i, j] for i in range(12) for j in range(i + 2, 12))
        prog_m = sum(Q_male[i, j] for i in range(12) for j in range(i + 2, 12))
        assert prog_m > prog_f


# ---------------------------------------------------------------------------
# 3. Gillespie simulator
# ---------------------------------------------------------------------------


class TestSimulateTrajectory:
    def _run(self, Q_base, state_defs, seed=0):
        rng = np.random.default_rng(seed)
        return simulate_patient_trajectory(
            Q_pat=Q_base,
            state_defs=state_defs,
            rng=rng,
            initial_state=0,
            max_time_years=10.0,
            min_time_years=2.0,
            dropout_rate=0.02,
            mean_interval_rrms=7.0,
            mean_interval_spms=9.0,
            post_relapse_interval=1.5,
            jitter_months=2.0,
            edss_noise_sigma=0.30,
        )

    def test_returns_nonempty(self, Q_base, state_defs):
        records = self._run(Q_base, state_defs)
        assert len(records) >= 1

    def test_states_in_valid_range(self, Q_base, state_defs):
        n = len(state_defs)
        records = self._run(Q_base, state_defs)
        for r in records:
            assert 0 <= r["state"] < n

    def test_edss_in_valid_range(self, Q_base, state_defs):
        records = self._run(Q_base, state_defs)
        for r in records:
            assert 0.0 <= r["EDSS"] <= 9.5
            assert 0.0 <= r["EDSS_true"] <= 9.5

    def test_edss_on_half_grid(self, Q_base, state_defs):
        records = self._run(Q_base, state_defs)
        for r in records:
            # Should be multiple of 0.5
            assert abs(round(r["EDSS"] * 2) - r["EDSS"] * 2) < 0.01

    def test_visit_times_monotone(self, Q_base, state_defs):
        records = self._run(Q_base, state_defs)
        times = [r["time_years"] for r in records]
        assert times == sorted(times)

    def test_relapse_binary(self, Q_base, state_defs):
        records = self._run(Q_base, state_defs)
        for r in records:
            assert r["relapse"] in {0, 1}

    def test_reproducibility_same_seed(self, Q_base, state_defs):
        r1 = self._run(Q_base, state_defs, seed=7)
        r2 = self._run(Q_base, state_defs, seed=7)
        assert len(r1) == len(r2)
        for a, b in zip(r1, r2):
            assert a["state"] == b["state"]
            assert abs(a["EDSS_true"] - b["EDSS_true"]) < 1e-9

    def test_different_seeds_differ(self, Q_base, state_defs):
        r1 = self._run(Q_base, state_defs, seed=1)
        r2 = self._run(Q_base, state_defs, seed=2)
        # Very unlikely to be identical
        states_1 = [r["state"] for r in r1]
        states_2 = [r["state"] for r in r2]
        assert states_1 != states_2 or len(r1) != len(r2)


# ---------------------------------------------------------------------------
# 4. Full dataset generator
# ---------------------------------------------------------------------------


class TestGenerateDataset:
    def test_schema_columns(self, small_df):
        required = {
            "patient_id", "visit_num", "visit_date", "days_from_onset",
            "EDSS", "EDSS_true", "state", "state_label", "relapse", "phase",
            "age_at_visit", "sex", "age_at_onset", "disease_duration_yr",
            "DMT_status", "latent_class",
        }
        assert required.issubset(set(small_df.columns))

    def test_patient_count(self, small_df, minimal_cfg):
        assert small_df["patient_id"].nunique() == minimal_cfg["cohort"]["n_patients"]

    def test_edss_range(self, small_df):
        assert small_df["EDSS"].min() >= 0.0
        assert small_df["EDSS"].max() <= 9.5

    def test_edss_half_grid(self, small_df):
        vals = small_df["EDSS"].to_numpy()
        np.testing.assert_allclose(vals * 2, np.round(vals * 2), atol=0.01)

    def test_state_range(self, small_df, minimal_cfg):
        n = minimal_cfg["states"]["n_states"]
        assert small_df["state"].min() >= 0
        assert small_df["state"].max() < n

    def test_relapse_binary(self, small_df):
        assert set(small_df["relapse"].unique()).issubset({0, 1})

    def test_sex_values(self, small_df):
        assert set(small_df["sex"].unique()).issubset({"F", "M"})

    def test_sex_ratio_approx(self, small_df):
        # ~68% F – allow ±10% for small cohort
        female_frac = (small_df.drop_duplicates("patient_id")["sex"] == "F").mean()
        assert 0.50 < female_frac < 0.88

    def test_phase_values(self, small_df):
        assert set(small_df["phase"].unique()).issubset({"RRMS", "SPMS", "ABS"})

    def test_latent_class_values(self, small_df):
        assert set(small_df["latent_class"].unique()).issubset(
            {"stable", "moderate", "aggressive"}
        )

    def test_dmt_values(self, small_df):
        assert set(small_df["DMT_status"].unique()).issubset(
            {"none", "moderate_dmt", "high_dmt"}
        )

    def test_age_at_visit_increases(self, small_df):
        for _, grp in small_df.groupby("patient_id"):
            ages = grp.sort_values("disease_duration_yr")["age_at_visit"].to_numpy()
            assert np.all(np.diff(ages) >= -1e-6)

    def test_visit_date_is_datetime(self, small_df):
        assert pd.api.types.is_datetime64_any_dtype(small_df["visit_date"])

    def test_days_from_onset_nonneg(self, small_df):
        assert (small_df["days_from_onset"] >= 0).all()

    def test_reproducibility(self, minimal_cfg):
        df1 = generate_synthetic_ms_dataset(cfg=minimal_cfg, seed=42)
        df2 = generate_synthetic_ms_dataset(cfg=minimal_cfg, seed=42)
        pd.testing.assert_frame_equal(
            df1.drop(columns=["visit_date"]).reset_index(drop=True),
            df2.drop(columns=["visit_date"]).reset_index(drop=True),
        )

    def test_different_seeds_differ(self, minimal_cfg):
        df1 = generate_synthetic_ms_dataset(cfg=minimal_cfg, seed=1)
        df2 = generate_synthetic_ms_dataset(cfg=minimal_cfg, seed=2)
        assert not df1["EDSS"].equals(df2["EDSS"])

    def test_aggressive_higher_edss_than_stable(self, minimal_cfg):
        df = generate_synthetic_ms_dataset(cfg=minimal_cfg, seed=7)
        stable_edss = df[df["latent_class"] == "stable"]["EDSS"].mean()
        aggr_edss = df[df["latent_class"] == "aggressive"]["EDSS"].mean()
        # Aggressive should have higher mean EDSS
        assert aggr_edss > stable_edss

    def test_visits_per_patient_reasonable(self, small_df):
        visits_per_pat = small_df.groupby("patient_id").size()
        assert float(visits_per_pat.median()) >= 5   # at least ~5 visits
        assert float(visits_per_pat.max()) <= 300    # upper bound sanity

    def test_onset_age_distribution(self, small_df):
        ages = small_df.drop_duplicates("patient_id")["age_at_onset"]
        assert 16 <= float(ages.min()) <= 30
        assert 40 <= float(ages.max()) <= 60
        assert 25 <= float(ages.mean()) <= 40


# ---------------------------------------------------------------------------
# 5. Epidemiological plausibility (moderate cohort)
# ---------------------------------------------------------------------------


class TestEpidemiologicalPlausibility:
    """Run on slightly larger cohort to get meaningful statistics."""

    @pytest.fixture
    def medium_df(self, minimal_cfg):
        cfg = dict(minimal_cfg)
        cfg["cohort"] = dict(minimal_cfg["cohort"])
        cfg["cohort"]["n_patients"] = 100
        cfg["visits"] = dict(minimal_cfg["visits"])
        cfg["visits"]["max_followup_years"] = 15.0
        return generate_synthetic_ms_dataset(cfg=cfg, seed=0)

    def test_spms_conversion_nonzero(self, medium_df):
        spms_pats = medium_df[medium_df["phase"] == "SPMS"]["patient_id"].nunique()
        total = medium_df["patient_id"].nunique()
        # At least 5% convert in 15 years (realistically 30–50%)
        assert spms_pats / total > 0.05

    def test_relapse_rate_decreases_with_progression(self, medium_df):
        # RRMS visits should have higher relapse fraction than SPMS visits
        rrms_relapse = medium_df[medium_df["phase"] == "RRMS"]["relapse"].mean()
        spms_df = medium_df[medium_df["phase"] == "SPMS"]
        if len(spms_df) > 10:
            spms_relapse = spms_df["relapse"].mean()
            # RRMS relapse rate should be >= SPMS (on average)
            assert rrms_relapse >= spms_relapse - 0.1

    def test_follow_up_within_bounds(self, medium_df, minimal_cfg):
        fu = medium_df.groupby("patient_id")["disease_duration_yr"].max()
        # The hard invariant is the upper follow-up horizon. The last *visit*
        # time can fall below the min-follow-up floor by up to one (sparse)
        # inter-visit interval, so we check the horizon strictly and the
        # typical patient loosely rather than asserting a tight lower bound.
        assert float(fu.max()) <= 15.0 + 0.1
        assert float(fu.median()) >= 1.0


# ---------------------------------------------------------------------------
# 6. Validation helpers
# ---------------------------------------------------------------------------


class TestValidationHelpers:
    def test_compare_empirical_to_Q_shape(self, small_df, minimal_cfg, Q_base):
        result = compare_empirical_to_Q(small_df, Q_base)
        assert "from_state" in result.columns
        assert "Q_true" in result.columns
        assert "Q_empirical" in result.columns
        assert len(result) > 0

    def test_empirical_Q_nonneg(self, small_df, Q_base):
        result = compare_empirical_to_Q(small_df, Q_base)
        assert (result["Q_empirical"] >= 0).all()


# ---------------------------------------------------------------------------
# 7. Helper function unit tests
# ---------------------------------------------------------------------------


class TestHelperFunctions:
    def test_snap_to_grid_0_5(self):
        assert _snap_to_grid(2.3, 0.5) == 2.5
        assert _snap_to_grid(2.2, 0.5) == 2.0
        assert _snap_to_grid(0.0, 0.5) == 0.0
        assert _snap_to_grid(9.5, 0.5) == 9.5

    def test_snap_to_grid_unity(self):
        assert _snap_to_grid(3.7, 1.0) == 4.0

    def test_snap_to_grid_zero_step(self):
        # zero step → no rounding
        assert _snap_to_grid(2.37, 0.0) == 2.37

    def test_state_at_time_before_first_event(self):
        events = [(0.0, 0), (2.0, 1), (5.0, 3)]
        assert _state_at_time(events, 0.0) == 0

    def test_state_at_time_between_events(self):
        events = [(0.0, 0), (2.0, 1), (5.0, 3)]
        assert _state_at_time(events, 3.0) == 1

    def test_state_at_time_after_last(self):
        events = [(0.0, 0), (2.0, 1), (5.0, 3)]
        assert _state_at_time(events, 8.0) == 3

    def test_relapse_in_window_true(self):
        state_defs = _DEFAULT_STATE_DEFS
        # State 1 is relapse-active
        events = [(0.0, 0), (1.0, 1), (2.0, 0)]
        assert _relapse_in_window(events, 1.1, state_defs, window=0.25)

    def test_relapse_in_window_false(self):
        state_defs = _DEFAULT_STATE_DEFS
        # Only inactive states
        events = [(0.0, 0), (2.0, 2)]
        assert not _relapse_in_window(events, 1.0, state_defs, window=0.25)

    def test_relapse_outside_window(self):
        state_defs = _DEFAULT_STATE_DEFS
        events = [(0.0, 0), (5.0, 1)]
        # Visit at t=1.0, relapse at t=5.0 → outside window
        assert not _relapse_in_window(events, 1.0, state_defs, window=0.25)


# ---------------------------------------------------------------------------
# 8. Edge cases
# ---------------------------------------------------------------------------


class TestEdgeCases:
    def test_single_patient(self, minimal_cfg):
        cfg = dict(minimal_cfg)
        cfg["cohort"] = dict(minimal_cfg["cohort"])
        cfg["cohort"]["n_patients"] = 1
        df = generate_synthetic_ms_dataset(cfg=cfg, seed=0)
        assert df["patient_id"].nunique() == 1
        assert len(df) >= 1

    def test_very_short_followup(self, minimal_cfg):
        cfg = dict(minimal_cfg)
        cfg["visits"] = dict(minimal_cfg["visits"])
        cfg["visits"]["max_followup_years"] = 2.0
        cfg["visits"]["min_followup_years"] = 1.0
        df = generate_synthetic_ms_dataset(cfg=cfg, seed=0)
        fu = df.groupby("patient_id")["disease_duration_yr"].max()
        assert float(fu.max()) <= 2.1

    def test_override_n_patients(self, minimal_cfg):
        df = generate_synthetic_ms_dataset(cfg=minimal_cfg, n_patients=5, seed=0)
        assert df["patient_id"].nunique() == 5

    def test_override_seed(self, minimal_cfg):
        df1 = generate_synthetic_ms_dataset(cfg=minimal_cfg, n_patients=5, seed=10)
        df2 = generate_synthetic_ms_dataset(cfg=minimal_cfg, n_patients=5, seed=10)
        pd.testing.assert_frame_equal(
            df1.drop(columns=["visit_date"]).reset_index(drop=True),
            df2.drop(columns=["visit_date"]).reset_index(drop=True),
        )

    def test_all_stable_class(self, minimal_cfg):
        cfg = dict(minimal_cfg)
        cfg["cohort"] = dict(minimal_cfg["cohort"])
        cfg["cohort"]["class_fractions"] = {"stable": 1.0, "moderate": 0.0, "aggressive": 0.0}
        cfg["cohort"]["n_patients"] = 10
        df = generate_synthetic_ms_dataset(cfg=cfg, seed=0)
        assert set(df["latent_class"].unique()) == {"stable"}
