# Changelog

## [1.0.0] - 2026-09-28
### Added
- Installable package `ms_synth` (`pyproject.toml`, src layout) with extras
  `analysis`, `percolation`, `parquet`, `dev`; console command `ms-synth-generate`.
- Known-truth interface: `generate_synthetic_ms_dataset(..., return_latent=True)`
  returns exact latent paths and patient generators without altering the cohort.
- Estimators in `ms_synth.msm_fit`: `estimate_Q_crude`, `estimate_Q_panel`
  (panel-likelihood MLE, vectorised; ~35x faster), `estimate_Q_from_paths` (oracle).
- `ms_synth.summaries`: single source of truth for all reported definitions.
- `scripts/benchmark_estimators.py` (resumable ADEMP simulation study),
  `scripts/make_figures.py` (all manuscript figures), `scripts/verify_reproduction.py`.
- `reference/` expected outputs + SHA-256 checksum; `requirements-lock.txt`;
  `Makefile`; `reproduce.sh`; `CITATION.cff`; GitHub Actions CI; Dockerfile.
- Option to disable outcome-dependent (post-relapse) visits
  (`visits.post_relapse_visit_months: 0`).
- 10 new tests (68 total).
### Changed
- Default CLI output is CSV (Parquet previously required an undeclared dependency).
- Positive EDSS slope now requires slope > 1e-9 yr^-1: a numerically flat
  trajectory previously flipped sign across BLAS builds (73.6% vs 73.8%).
- `msm_fit` method `"em"` renamed `"panel"` (alias kept): it is a direct
  numerical MLE, not an EM algorithm.
- Removed a pandas-deprecated `groupby.apply` call.
### Fixed
- Removed incorrect/non-existent literature citations and the claim that the
  generator matrix was fitted to registry data (it is calibrated to aggregates).
