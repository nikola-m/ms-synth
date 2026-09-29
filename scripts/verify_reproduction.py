#!/usr/bin/env python3
"""
Verify that this installation reproduces the published outputs.

1. Regenerates the reference cohort (N = 500, seed 42) and checks its
   SHA-256 against reference/SHA256SUMS (byte-identical output).
2. Recomputes every manuscript statistic and compares with
   reference/expected_results.json: exact equality, except panel-likelihood
   estimates (numerical optimisation; abs. tolerance 1e-3).
3. Optionally (--benchmark) compares results/benchmark_summary.csv with
   reference/expected_benchmark_summary.csv (tolerance 1e-3).

Exit status 0 = reproduced; 1 = mismatch.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import logging
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REF = ROOT / "reference"
TOL_PANEL = 1e-3


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def compare(exp, got, path=""):
    if isinstance(exp, dict):
        for k in exp:
            if k == "ms_synth_version":
                continue
            yield from compare(exp[k], got.get(k) if isinstance(got, dict) else None, f"{path}/{k}")
    elif isinstance(exp, list):
        if not isinstance(got, list) or len(got) != len(exp):
            yield (path, exp, got); return
        for i, (e, g) in enumerate(zip(exp, got)):
            yield from compare(e, g, f"{path}[{i}]")
    elif isinstance(exp, float) and isinstance(got, (int, float)) and "panel" in path.lower():
        if abs(exp - got) > TOL_PANEL:
            yield (path, exp, got)
    elif exp != got:
        yield (path, exp, got)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--benchmark", action="store_true", help="also check the estimator benchmark summary")
    a = ap.parse_args()
    ok = True
    logging.disable(logging.CRITICAL)
    from ms_synth.synthetic_data import generate_synthetic_ms_dataset
    from ms_synth.utils import load_config

    # 1. byte-identical cohort
    expected = dict(line.split()[::-1] for line in (REF / "SHA256SUMS").read_text().splitlines() if line.strip())
    out = ROOT / "data" / "reference_cohort_seed42_n500.csv"
    out.parent.mkdir(exist_ok=True)
    df = generate_synthetic_ms_dataset(cfg=load_config(str(ROOT / "configs" / "synthetic.yaml")),
                                       n_patients=500, seed=42)
    df.to_csv(out, index=False)
    got = sha256(out)
    same = got == expected[out.name]
    ok &= same
    print(f"[{'PASS' if same else 'FAIL'}] cohort SHA-256 {got[:16]}... ({len(df)} visit records)")

    # 2. statistics
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td) / "results.json"
        subprocess.run([sys.executable, str(ROOT / "scripts" / "analyze_manuscript_stats.py"),
                        "--out", str(tmp)], check=True, capture_output=True)
        diffs = list(compare(json.loads((REF / "expected_results.json").read_text()),
                             json.loads(tmp.read_text())))
    ok &= not diffs
    print(f"[{'PASS' if not diffs else 'FAIL'}] manuscript statistics ({len(diffs)} mismatches)")
    for d in diffs[:20]:
        print("       ", d)

    # 3. benchmark (optional; requires scripts/benchmark_estimators.py to have run)
    if a.benchmark:
        import pandas as pd
        e = pd.read_csv(REF / "expected_benchmark_summary.csv")
        g = pd.read_csv(ROOT / "results" / "benchmark_summary.csv")
        keys = ["scenario", "population", "visits", "n_patients", "interval_mult", "estimator"]
        m = e.merge(g, on=keys, suffixes=("_e", "_g"))
        bad = m[(m.bias_mean_e - m.bias_mean_g).abs().gt(1e-3) | (m.abs_err_mean_e - m.abs_err_mean_g).abs().gt(1e-3)]
        good = len(m) == len(e) and bad.empty
        ok &= good
        print(f"[{'PASS' if good else 'FAIL'}] estimator benchmark summary ({len(m)}/{len(e)} cells matched, {len(bad)} outside tolerance)")

    print("REPRODUCED" if ok else "NOT REPRODUCED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
