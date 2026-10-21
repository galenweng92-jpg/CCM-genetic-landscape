# CCM genetic landscape: aggregate statistical reproducibility

This repository candidate accompanies *The genetic landscape of familial and sporadic cerebral cavernous malformations*. It recomputes descriptive and association statistics from frozen aggregate counts corresponding to Supplemental Table S8.

**Status: private-repository initial-deposit package, prepared on 2026-10-21; no public release or archival DOI.** The authorized target is [galenweng92-jpg/CCM-genetic-landscape](https://github.com/galenweng92-jpg/CCM-genetic-landscape), with private visibility. This local package does not by itself confirm that the remote upload has completed. Code and data licensing require author approval before public deposition.

## Scope

The inputs are study-level or predefined-subset counts. The package reproduces Wilson confidence intervals, gene-specific detection summaries and sensitivity analyses, paired-lesion summaries, Fisher exact tests, Mantel–Haenszel odds ratios, CMH tests, and exact conditional stratified tests.

It does **not** redo literature screening, case extraction, record-level eligibility decisions, variant annotation, overlap adjudication, genomic plotting, transcriptomic modeling, or mechanism-figure construction. It does not constitute release of all analysis code for the study. The full author-controlled audit package retains record-level workbooks and the associated screening scripts separately.

No participant/specimen/family identifiers, individual-level genetic records, clinical free text, manuscript comments, credentials, private filesystem paths, third-party gene-set memberships, article PDFs, or BioRender assets are included.

## Run

Requirements:

- Python 3.12 or later; standard library only.
- R 4.5 with `jsonlite`; statistical tests use base R `stats`.
- `Rscript` on `PATH`, or provide its location with `--rscript`.

```sh
Rscript -e 'install.packages("jsonlite", repos="https://cloud.r-project.org")'
python3 scripts/reproduce_statistics.py
python3 scripts/verify_results.py
python3 scripts/audit_release.py
```

The installation command is for a new environment. The analysis itself operates offline. Input files are read-only; results are written to `results/recomputed`. A custom output directory is available:

```sh
python3 scripts/reproduce_statistics.py --output-dir results/recomputed --rscript Rscript
```

`verify_results.py` compares all 16 reference CSV tables cell by cell, including all six Table S8 tables. `audit_release.py` scans the release tree for disallowed file types and obvious sensitive-content patterns, and writes its file manifest and checksums. These scans support author review; they are not a legal or absolute privacy certification.

## Files

| Path | Contents |
|---|---|
| `data/aggregate_inputs` | Seven count-only inputs. No individual case rows. |
| `data/source_publications.csv` | Public study names and bibliographic URLs. |
| `data/AGGREGATE_PROVENANCE.json` | Frozen source version, aggregate source hashes and field-selection provenance. |
| `scripts/reproduce_statistics.py` | Aggregate-count analysis and independent Python cross-checks. |
| `scripts/exact_stratified_tests.R` | Fisher, MH, CMH, exact stratified analyses, and Wilson cross-checks. |
| `scripts/verify_results.py` | Exact comparison with frozen reference tables. |
| `scripts/audit_release.py` | File-content checks and SHA-256 manifest generation. |
| `results/reference` | Sixteen frozen aggregate result tables from the author-controlled analysis. |
| `results/recomputed` | Freshly computed results. |
| `validation` | Executed environment and validation reports. |
| `metadata` | Local Zenodo draft metadata and release status. No publication has been initiated. |

## Reading the data

`positive` is the number of positive records and `eligible_n` is the gene- or endpoint-specific denominator. `total_records` can exceed `eligible_n` because unavailable/untested records are not treated as negatives. Different sensitivity subsets overlap and must not be added together.

Association cells are `a = driver-positive/classical-CCM-positive`, `b = driver-positive/classical-CCM-negative`, `c = driver-negative/classical-CCM-positive`, and `d = driver-negative/classical-CCM-negative`. The same conventions apply to each study stratum. Complete-panel analyses require all three classical CCM genes to be evaluable; partial-panel sensitivity analyses require at least one. A negative composite in a partial panel is not proven negativity in all three genes.

Source-linked lesion-testing records are not certified independent patients. Study stratification does not resolve within-study selection or cross-report overlap. Excluding Ren 2025 and/or Hong 2021 is a sensitivity analysis, not established deduplication. Descriptive binomial confidence intervals are not population prevalence estimates. Exploratory associations do not establish molecular exclusivity or cooperation.

## Rights and citation

Read `RIGHTS_AND_PROVENANCE.md` before redistribution. No blanket license has been assigned. The project-authored code and the literature-derived aggregate data have separate provenance and require separate licensing review. Third-party source articles and database memberships are not relicensed by this package.

`CITATION.cff` retains the manuscript author names/order from the local manuscript package. It supplies neither ORCIDs nor affiliations that were not verified. The repository URL identifies the authorized private-deposit target, not a public release. After a real public release, insert the assigned archival DOI and update the manuscript's Data and code availability statement to describe this package's limited scope.
