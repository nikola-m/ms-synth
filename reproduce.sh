#!/usr/bin/env bash
# One-shot reproduction of every result, table and figure in the manuscript.
# Usage: bash reproduce.sh            (full, ~20 min on one CPU core)
#        QUICK=1 bash reproduce.sh    (skip the 15-min estimator benchmark)
set -euo pipefail
cd "$(dirname "$0")"
python -m venv .venv && . .venv/bin/activate
python -m pip install --upgrade pip >/dev/null
make install-locked
make test
make stats
if [ "${QUICK:-0}" = "1" ]; then
  make figures && make verify
else
  make benchmark && make figures && make verify-all
fi
