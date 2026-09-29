# ms-synth: one command per reproducibility step. `make all` reruns everything.
PY ?= python

.PHONY: install install-locked test cohort stats benchmark figures verify verify-all all clean

install:            ## editable install with analysis + test extras
	$(PY) -m pip install -e ".[dev]"

install-locked:     ## exact dependency versions used for the paper
	$(PY) -m pip install -r requirements-lock.txt && $(PY) -m pip install --no-deps -e .

test:               ## unit tests (~20 s)
	$(PY) -m pytest

cohort:             ## reference cohort, N = 500, seed 42 -> data/
	ms-synth-generate --config configs/synthetic.yaml --n-patients 500 --seed 42 \
	    --output data/reference_cohort_seed42_n500.csv --validate

stats:              ## every statistic in the manuscript -> results/results.json
	$(PY) scripts/analyze_manuscript_stats.py

benchmark:          ## estimator benchmark, 85 tasks (~15 min on one core; resumable)
	$(PY) scripts/benchmark_estimators.py --reps 5

figures:            ## all manuscript figures -> results/figures/*.pdf
	$(PY) scripts/make_figures.py

verify:             ## compare with reference/ (cohort checksum + statistics)
	$(PY) scripts/verify_reproduction.py

verify-all: benchmark
	$(PY) scripts/verify_reproduction.py --benchmark

all: test stats benchmark figures verify-all

clean:
	rm -rf results/bench_parts data/*.csv .pytest_cache
