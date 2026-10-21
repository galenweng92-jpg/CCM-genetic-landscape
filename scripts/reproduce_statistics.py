#!/usr/bin/env python3
"""Recompute CCM descriptive and association statistics from aggregate counts.

This entry point does not perform source-paper extraction or record selection.
No participant-level data, credentials, or network connection are required.
"""
import argparse
import csv
import hashlib
import json
import math
import platform
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--input-dir", type=Path, default=ROOT / "data/aggregate_inputs")
parser.add_argument("--output-dir", type=Path, default=ROOT / "results/recomputed")
parser.add_argument("--rscript", default="Rscript")
args = parser.parse_args()
SRC, OUT = args.input_dir.resolve(), args.output_dir.resolve()
OUT.mkdir(parents=True, exist_ok=True)
Z = 1.959963984540054


def read_csv(name):
    with (SRC / name).open(encoding="utf-8-sig", newline="") as f:
        rows = list(csv.DictReader(f))
    for r in rows:
        for k, v in r.items():
            if k in {"positive", "eligible_n", "total_records", "not_detected", "unavailable",
                     "a", "b", "c", "d", "n", "a_driver_pos_ccm_pos", "b_driver_pos_ccm_neg",
                     "c_driver_neg_ccm_pos", "d_driver_neg_ccm_neg"}:
                r[k] = int(v)
            elif k in {"at_least_23_source_records", "informative_for_CMH"}:
                assert v in {"True", "False"}
                r[k] = v == "True"
    return rows


def save_csv(name, rows):
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with (OUT / name).open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def wilson(x, n):
    assert isinstance(x, int) and isinstance(n, int) and 0 <= x <= n
    if n == 0:
        return dict(positive=x, eligible_n=n, fraction=None, ci95_lower=None, ci95_upper=None)
    p = x / n
    div = 1 + Z * Z / n
    center = (p + Z * Z / (2 * n)) / div
    half = Z / div * math.sqrt(p * (1 - p) / n + Z * Z / (4 * n * n))
    return dict(positive=x, eligible_n=n, fraction=p,
                ci95_lower=center-half, ci95_upper=center+half)


rate_names = [
    "familial_composition_wilson.csv", "paired_endpoints_wilson.csv",
    "sporadic_detection_sensitivity_wilson.csv",
    "map3k3_pik3ca_cofinding_sensitivity_wilson.csv",
    "study_specific_detection_wilson.csv",
]
rates = {}
for name in rate_names:
    rows = read_csv(name.replace("_wilson", "_counts"))
    for r in rows:
        r.update(wilson(r["positive"], r["eligible_n"]))
    rates[name] = rows
    save_csv(name, rows)
family_rows, paired_rows, rate_rows, co_rows, study_rows = [rates[n] for n in rate_names]
range_rows = []
for gene in ["MAP3K3", "PIK3CA"]:
    rr = [r for r in study_rows if r["gene"] == gene and r["at_least_23_source_records"] and r["eligible_n"]]
    lo, hi = min(rr, key=lambda r: r["fraction"]), max(rr, key=lambda r: r["fraction"])
    range_rows.append(dict(gene=gene, threshold="study source records >=23; gene evaluable N>0",
                           contributing_studies=len(rr), min_study=lo["study"], min_positive=lo["positive"],
                           min_eligible_n=lo["eligible_n"], min_fraction=lo["fraction"],
                           max_study=hi["study"], max_positive=hi["positive"], max_eligible_n=hi["eligible_n"],
                           max_fraction=hi["fraction"]))
save_csv("between_study_ranges.csv", range_rows)

analysis_rows = read_csv("association_2x2_counts.csv")
strata_rows = read_csv("association_study_strata.csv")
for r in analysis_rows:
    ss = [s for s in strata_rows if s["analysis_id"] == r["analysis_id"]]
    for aggregate_col, stratum_col in [
        ("a_driver_pos_ccm_pos", "a"), ("b_driver_pos_ccm_neg", "b"),
        ("c_driver_neg_ccm_pos", "c"), ("d_driver_neg_ccm_neg", "d"), ("n", "n")
    ]:
        assert r[aggregate_col] == sum(s[stratum_col] for s in ss), r["analysis_id"]
for s in strata_rows:
    assert s["n"] == sum(s[k] for k in "abcd")
    assert s["informative_for_CMH"] == (min(s["a"]+s["b"], s["c"]+s["d"], s["a"]+s["c"], s["b"]+s["d"]) > 0)
save_csv("association_2x2_counts.csv", analysis_rows)
save_csv("association_study_strata.csv", strata_rows)
subprocess.run([args.rscript, str(Path(__file__).with_name("exact_stratified_tests.R")), str(OUT)], check=True)

results = list(csv.DictReader((OUT / "association_exact_stratified_results.csv").open()))
validation, association_ci, publication = [], [], []
for r in results:
    a, b, c, d = [int(r[k]) for k in ["a_driver_pos_ccm_pos", "b_driver_pos_ccm_neg",
                                     "c_driver_neg_ccm_pos", "d_driver_neg_ccm_neg"]]
    n, n1, m1 = a+b+c+d, a+b, a+c
    den = math.comb(n, n1)
    probs = [math.comb(m1, k)*math.comb(n-m1, n1-k)/den
             for k in range(max(0, n1-(n-m1)), min(n1, m1)+1)]
    observed = math.comb(m1, a)*math.comb(n-m1, n1-a)/den
    pf = sum(p for p in probs if p <= observed*(1+1e-7))
    assert abs(pf-float(r["fisher_two_sided_p"])) < 1e-10
    ss = [s for s in strata_rows if s["analysis_id"] == r["analysis_id"] and s["informative_for_CMH"]]
    mh = sum(s["a"]*s["d"]/s["n"] for s in ss)/sum(s["b"]*s["c"]/s["n"] for s in ss)
    assert abs(mh-float(r["mh_odds_ratio"])) < 1e-10
    validation.append(dict(check=r["analysis_id"], fisher_agrees=True, mh_agrees=True))
    for group, x, nn in [("Driver-positive", a, a+b), ("Driver-negative", c, c+d)]:
        association_ci.append(dict(analysis_id=r["analysis_id"], group=group, **wilson(x, nn)))
    if r["population"] not in {"all_source_records", "known_absent_family_history"}:
        continue
    publication.append(dict(
        driver=r["driver"], ccm_assay_scope=r["scope"], population=r["population"],
        ccm_positive_among_driver_positive=a, driver_positive_evaluable_n=a+b,
        fraction_among_driver_positive=a/(a+b),
        ccm_positive_among_driver_negative=c, driver_negative_evaluable_n=c+d,
        fraction_among_driver_negative=c/(c+d),
        fisher_two_sided_p=r["fisher_two_sided_p"],
        fisher_conditional_odds_ratio=r["fisher_conditional_odds_ratio"],
        fisher_exact_ci95_lower=r["fisher_exact_ci95_lower"],
        fisher_exact_ci95_upper=r["fisher_exact_ci95_upper"],
        mh_odds_ratio=r["mh_odds_ratio"],
        mh_asymptotic_ci95_lower=r["mh_ci95_lower"] or "Not estimable",
        mh_asymptotic_ci95_upper=r["mh_ci95_upper"] or "Not estimable",
        cmh_p_no_continuity_correction=r["cmh_p_uncorrected"],
        stratified_exact_conditional_or=r["stratified_exact_conditional_or"],
        stratified_exact_ci95_lower=r["stratified_exact_ci95_lower"],
        stratified_exact_ci95_upper=r["stratified_exact_ci95_upper"],
        stratified_exact_two_sided_p=r["stratified_exact_two_sided_p"],
        informative_studies=r["strata_informative"], informative_records=r["informative_n"],
        analysis_role="Complete-panel analysis" if r["scope"] == "all_three_CCM_genes_evaluable"
                      else "Partial-panel sensitivity analysis",
        unit="Source-linked lesion-testing record",
        interpretation="Exploratory source-record association; no causal or independent-patient inference",
    ))
save_csv("association_group_wilson.csv", association_ci)
save_csv("TableS8_association_results.csv", publication)
save_csv("TableS8_detection_sensitivity.csv", rate_rows)
save_csv("TableS8_family_and_paired_CI.csv", family_rows+paired_rows)
save_csv("TableS8_study_detection.csv", study_rows)
save_csv("TableS8_study_ranges.csv", range_rows)
save_csv("TableS8_driver_cofindings.csv", co_rows)
(OUT / "independent_validation.json").write_text(json.dumps(validation, indent=2) + "\n")
(OUT / "analysis_audit.json").write_text(json.dumps({
    "source_version": "2026-10-16 frozen-source tables; aggregates prepared 2026-10-20",
    "scope": "Aggregate counts to statistical estimates only; no source extraction or record filtering rerun.",
    "input_hashes": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(SRC.glob("*.csv"))},
    "python_version": platform.python_version(),
    "independent_Fisher_and_MH_checks": len(validation),
    "limitations": [
        "Source-linked records are not certified independent patients.",
        "Predefined subsets were aggregated before release; this package cannot audit record-level eligibility.",
        "Stratification does not remove unresolved cross-report overlap or within-study selection.",
        "Complete-panel eligibility describes recorded assay scope, not uniform sequencing sensitivity.",
        "Ren/Hong exclusions are sensitivity analyses, not established deduplication.",
        "Descriptive binomial confidence intervals are not population prevalence estimates.",
    ],
}, indent=2) + "\n")
print(json.dumps({"aggregate_reanalysis_completed": True, "association_checks": len(validation)}, indent=2))
