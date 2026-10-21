#!/usr/bin/env python3
"""Audit this release tree and write its environment, inventory and checksums.

Pattern checks are a review aid, not an absolute privacy or rights guarantee.
No source files, remote accounts or directories outside this tree are modified.
"""
import argparse
import ast
import csv
import hashlib
import json
import platform
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--rscript", default="Rscript")
args = parser.parse_args()
VALIDATION = ROOT / "validation"
VALIDATION.mkdir(exist_ok=True)
GENERATED_MANIFESTS = {"FILE_MANIFEST.csv", "CHECKSUMS.sha256"}
ALLOW_SUFFIXES = {".py", ".r", ".csv", ".json", ".md", ".cff", ".txt", ".sha256"}
PATTERNS = {
    "private_absolute_path": r"/(?:Users|home)/[^/\s]+/|[A-Za-z]:\\Users\\",
    "credential_value": r"(?:ghp_[A-Za-z0-9]{20,}|github_pat_[A-Za-z0-9_]{30,}|AKIA[A-Z0-9]{16}|-----BEGIN (?:RSA |OPENSSH |EC )?PRIVATE KEY-----|(?i:access_token|api_key|client_secret)\s*[:=]\s*[\x22\x27][A-Za-z0-9_/-]{12,})",
}
CASE_HEADERS = {"record id", "record_id", "patient id", "patient_id", "original patient id",
                "original family id", "patient / specimen id", "source_excel_row", "source locator",
                "evidence / notes", "variant reported", "protein reported"}


def payload_files():
    return sorted(p for p in ROOT.rglob("*") if (p.is_file() or p.is_symlink())
                  and ".git" not in p.relative_to(ROOT).parts
                  and str(p.relative_to(ROOT)) not in GENERATED_MANIFESTS)


env = {
    "checked_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    "python_version": platform.python_version(),
    "python_implementation": platform.python_implementation(),
    "operating_system": platform.system(), "operating_system_version": platform.mac_ver()[0] or platform.release(),
    "architecture": platform.machine(),
    "python_dependencies": "Standard library only",
}
try:
    value = subprocess.check_output([args.rscript, "-e",
        'cat(R.version.string,"\\n",as.character(packageVersion("jsonlite")),"\\n",sep="")'],
        text=True, stderr=subprocess.PIPE).splitlines()
    env.update(R_version=value[0], jsonlite_version=value[1], R_available=True)
except (OSError, subprocess.CalledProcessError, IndexError) as exc:
    env.update(R_available=False, R_check_error_type=type(exc).__name__)
(VALIDATION / "ENVIRONMENT.json").write_text(json.dumps(env, indent=2) + "\n")

findings, csv_checks, syntax_checks = [], [], []
for p in payload_files():
    rel = str(p.relative_to(ROOT))
    if p.is_symlink():
        findings.append({"file": rel, "issue": "symlink_not_permitted"})
        continue
    if p.name != ".gitignore" and p.suffix.lower() not in ALLOW_SUFFIXES:
        findings.append({"file": rel, "issue": "unsupported_or_prohibited_asset_type"})
    try:
        content = p.read_text(encoding="utf-8-sig")
    except (UnicodeError, OSError):
        findings.append({"file": rel, "issue": "not_inspectable_UTF8_text"})
        continue
    for name, pattern in PATTERNS.items():
        if re.search(pattern, content):
            findings.append({"file": rel, "issue": name})
    if p.suffix.lower() == ".csv":
        rows = list(csv.reader(content.splitlines()))
        forbidden = sorted({h.strip().lower() for h in rows[0]} & CASE_HEADERS) if rows else []
        csv_checks.append({"file": rel, "rows": max(0, len(rows)-1), "case_headers_found": forbidden})
        if forbidden:
            findings.append({"file": rel, "issue": "case_level_header"})
    if p.suffix == ".py":
        try:
            ast.parse(content, filename=rel)
            syntax_checks.append({"file": rel, "parsed": True})
        except SyntaxError:
            findings.append({"file": rel, "issue": "python_syntax_error"})

report = {
    "checked_at_utc": env["checked_at_utc"],
    "scope": "This release directory only; .git internals and the two generated checksum/manifest files are excluded.",
    "all_pattern_checks_passed": not findings,
    "findings": findings, "csv_checks": csv_checks, "python_syntax_checks": syntax_checks,
    "limitations": [
        "Pattern checks do not certify anonymity, absence of every possible secret, or redistribution rights.",
        "Only aggregate inputs are included; source extraction and individual-record filtering are not audited here.",
        "Some published study strata have small counts; author rights/privacy review remains required before public release.",
        "No license, repository publication or archival DOI is established by this audit.",
    ],
}
(VALIDATION / "RELEASE_AUDIT.json").write_text(json.dumps(report, indent=2) + "\n")
manifest = []
for p in payload_files():
    if p.is_symlink():
        continue
    rel = str(p.relative_to(ROOT))
    blob = p.read_bytes()
    manifest.append({"path": rel, "bytes": len(blob), "sha256": hashlib.sha256(blob).hexdigest(),
                     "role": rel.split("/")[0] if "/" in rel else "repository_metadata"})
with (ROOT / "FILE_MANIFEST.csv").open("w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["path", "bytes", "sha256", "role"])
    w.writeheader()
    w.writerows(manifest)
(ROOT / "CHECKSUMS.sha256").write_text("".join(f'{r["sha256"]}  {r["path"]}\n' for r in manifest))
print(json.dumps({"all_pattern_checks_passed": not findings, "payload_files": len(manifest),
                  "payload_bytes": sum(r["bytes"] for r in manifest), "manifest_files": 2,
                  "R_available": env["R_available"], "findings": findings}, indent=2))
sys.exit(1 if findings or not env["R_available"] else 0)
