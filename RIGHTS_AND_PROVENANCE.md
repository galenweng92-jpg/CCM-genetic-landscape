# Rights and provenance

This repository was made publicly accessible on 2026-10-21 following the project owner's authorization. No software or data reuse license is assigned. A public repository alone does not establish permission for unrestricted reuse.

## Project-authored code

- `scripts/reproduce_statistics.py` adapts the project's 2026-10-20 statistical workflow to aggregate-only inputs.
- `scripts/exact_stratified_tests.R` is retained from that workflow.
- `scripts/verify_results.py` and `scripts/audit_release.py` were prepared for this release candidate.

These scripts were prepared within this project, including AI-assisted drafting and subsequent numerical checks. No claim is made that manuscript authorship establishes each person's software contribution or resolves institutional ownership. The authors must approve the code's attribution and license. The code uses Python standard-library functions and R's `stats`/`jsonlite` APIs; those dependencies are not redistributed here.

## Literature-derived data and results

`data/aggregate_inputs` contains author-computed aggregate counts derived from the frozen 2026-10-16 supplemental record tables. It contains no case-level rows, source case/family identifiers, individual variants, or clinical narrative. `data/source_publications.csv` records public bibliographic links. `data/AGGREGATE_PROVENANCE.json` identifies the aggregate source snapshots and transformations.

The source papers, their tables, and their licenses remain with their respective rightsholders. This repository candidate does not include or grant rights over their PDFs, figures, full-text content, or original tables. The rights to redistribute the derived aggregate dataset should be reviewed independently from the code license. Some strata have small counts; the aggregates are not a formal anonymization or privacy guarantee.

`results/reference` and `results/recomputed` contain aggregate estimates only. Their provenance is the inputs and statistical methods described in the README.

## Deliberately excluded

Individual-level source workbooks and analysis ledgers; specimen/patient/family IDs; clinical free text; source case locators; original article PDFs; third-party gene-set memberships including MSigDB/Reactome-derived lists; BioRender elements, artwork, or licenses; unpublished manuscript comments; local account paths; credentials; raw sequencing and cell-level expression data.

The full local audit package is retained separately for author-controlled verification. This public candidate does not reproduce source extraction or individual-record eligibility and is not a complete study-wide code release.

## Tagged and archival releases

The GitHub repository is publicly accessible. For a tagged or Zenodo archival release, confirm the creator list and ownership, select appropriate separate code/data reuse licenses, and verify that the README and manuscript describe the aggregate-only scope accurately. Insert an archival DOI only after it has actually been assigned. No DOI or license is reserved or invented here.
