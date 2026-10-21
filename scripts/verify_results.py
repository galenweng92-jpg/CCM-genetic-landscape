#!/usr/bin/env python3
"""Compare every released reference table with a fresh aggregate-count rerun."""
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REF, OUT = ROOT / "results/reference", ROOT / "results/recomputed"
checks = []
for source in sorted(REF.glob("*.csv")):
    target = OUT / source.name
    assert target.exists(), target
    with source.open(encoding="utf-8-sig", newline="") as a, target.open(encoding="utf-8-sig", newline="") as b:
        first, second = list(csv.reader(a)), list(csv.reader(b))
    assert first == second, source.name
    checks.append({"file": source.name, "rows": len(first)-1, "all_CSV_values_identical": True})
assert len(checks) >= 15
audit = json.loads((OUT / "analysis_audit.json").read_text())
for name, expected in audit["input_hashes"].items():
    assert hashlib.sha256((ROOT / "data/aggregate_inputs" / name).read_bytes()).hexdigest() == expected
report = {"all_passed": True, "scope": "Aggregate counts to estimates; no case extraction rerun.",
          "comparisons": checks, "input_checksums_unchanged": True}
(ROOT / "validation/RESULTS_VALIDATION.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps({"all_passed": True, "reference_tables": len(checks),
                  "rows_compared": sum(c["rows"] for c in checks)}, indent=2))
