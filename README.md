# ms-synth

**A known-truth generator of synthetic longitudinal multiple sclerosis cohorts**

`ms-synth` simulates registry-like MS cohorts — irregular clinic visits, noisy
EDSS scores, relapses, disease-modifying therapy, dropout — from a
continuous-time multi-state Markov model (CT-MSM) whose parameters are
calibrated to published natural-history cohorts. Because every cohort is
generated from a fully specified process, the **truth is known**: the
generator matrix of every patient, the exact latent disease path, and each
patient's latent progression class. Methods fitted to the observed data can
therefore be scored against the truth — something no real registry allows.

Typical uses: benchmarking CT-MSM / panel-data estimators, survival and
milestone analyses, trajectory clustering and latent-class recovery,
testing pipelines before access to restricted registry data is granted,
teaching, and the percolation analysis of state-transition graphs that
originally motivated the project.

> No real patient data are used or contained. All cohorts are synthetic.

## Install

Python ≥ 3.10.

```bash
pip install -e ".[dev]"                       # latest compatible dependencies
# or, the exact environment used for the paper:
pip install -r requirements-lock.txt && pip install --no-deps -e .
```

The core generator needs only `numpy`, `scipy`, `pandas`, `pyyaml`. Extras:
`analysis` (lifelines, matplotlib), `percolation` (networkx), `parquet` (pyarrow).

## Quick start

```bash
ms-synth-generate --n-patients 500 --seed 42 --output data/cohort.csv
```

```python
from ms_synth.utils import load_config
from ms_synth.synthetic_data import generate_synthetic_ms_dataset, build_Q_from_config
from ms_synth.msm_fit import estimate_Q_panel, estimate_Q_crude, estimate_Q_from_paths

cfg = load_config("configs/synthetic.yaml")
df, truth = generate_synthetic_ms_dataset(cfg=cfg, n_patients=500, seed=42, return_latent=True)
# df:    one row per visit (observed data)
# truth: patient_id -> {"events": exact latent path, "follow_up": T, "Q_patient": generator}

Q_true = build_Q_from_config(cfg)
Q_panel = estimate_Q_panel(df, n_states=12)        # panel-likelihood MLE
Q_oracle = estimate_Q_from_paths(truth, n_states=12)  # upper bound from exact paths
```

`return_latent=True` consumes no random numbers: the observed cohort is
byte-identical with or without it.

## Reproducing the paper

```bash
bash reproduce.sh            # everything (~20 min on one core)
QUICK=1 bash reproduce.sh    # skip the 15-min estimator benchmark (~1 min)
```

or step by step: `make test`, `make stats`, `make benchmark`, `make figures`,
`make verify`. `make verify` regenerates the reference cohort (N = 500,
seed 42), checks its SHA-256 against `reference/SHA256SUMS`, and compares
every reported statistic with `reference/expected_results.json`.
The reference cohort is byte-identical across the dependency stacks we tested
(NumPy 2.4–2.5, SciPy 1.17–1.18, pandas 2.3–3.0); only panel-likelihood
estimates vary, in the fourth decimal (optimiser termination).

| Output | Script |
|---|---|
| `results/results.json` — every number in the manuscript | `scripts/analyze_manuscript_stats.py` |
| `results/benchmark_*.csv` — estimator benchmark | `scripts/benchmark_estimators.py` |
| `results/figures/fig1–5*.pdf` | `scripts/make_figures.py` |
| verification | `scripts/verify_reproduction.py` |

## Repository layout

```
src/ms_synth/        synthetic_data (generator), msm_fit (estimators),
                     summaries (shared definitions), cli, data, utils,
                     graph_utils / percolation / augmentation / visualization (optional)
configs/synthetic.yaml   12-state generator, covariate effects, visit process
scripts/             analysis, benchmark, figures, verification, percolation pipeline
reference/           expected outputs and checksums
tests/               68 unit tests
Validation.md        traceable natural-history parameter compilation
```

## Parameter provenance

Transition intensities are **hand-specified and calibrated** so that
aggregate outputs match published natural-history targets (London Ontario,
Lyon, British Columbia, Rennes, Big MS Data, MSBase, pivotal trials; see
`Validation.md` and the manuscript). They are not fitted to any registry's
patient-level data.

## Citation

See `CITATION.cff`. Please cite the software (with version) and the article.

## License

MIT
