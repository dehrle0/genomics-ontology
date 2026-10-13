#!/usr/bin/env python3
"""
lib/zero_curation_guard.py
Permanent Automated Gate Enforcing Zero Static Curation & Zero PII in Code and Skills.

Scans all Python source files (*.py, lib/*.py, scripts/*.py) and custom skills
(~/.gemini/config/skills/**/*) against:
  1. Regex patterns for rsIDs (rs\\d+), HGVS (p./c.), genomic coordinates (chrN:pos),
     and OMIM/MIM identifiers.
  2. Dynamic vocabulary extracted at runtime from local cohort JSON files and OS environment
     (patient names, sample IDs, home directory usernames, and all HGNC gene symbols in the cohort),
     ensuring the guard itself contains zero hardcoded PII or gene symbols.

Exits with code 1 and blocks execution / git commits if any violation is detected.
"""

import glob
import json
import os
import re
import sys

SCRIPT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HOME_DIR = os.path.expanduser("~")


def _build_dynamic_forbidden_patterns(reports_dir):
    """Build regex patterns dynamically from local environment and cohort JSONs (zero static strings)."""
    static_patterns = {
        "rsid_literal": re.compile(r"\brs\d+\b"),
        "hgvs_literal": re.compile(r"\b[pc]\.[A-Z][a-z]{2}\d+|\b[pc]\.\d+[ACGT>\-+_]+"),
        "coord_literal": re.compile(r"\bchr(?:\d+|X|Y|M|MT):\d+"),
        "omim_literal": re.compile(r"\b(?:OMIM|MIM):\d+"),
    }

    pii_tokens = set()
    user_name = os.path.basename(HOME_DIR)
    if user_name and len(user_name) >= 3:
        pii_tokens.add(re.escape(user_name))
        for part in re.split(r"[-_.\s]+", user_name):
            if len(part) >= 4:
                pii_tokens.add(re.escape(part))

    gene_tokens = set()
    # Common uppercase words in Python/HTML/CSS that are not cohort-specific curations
    allowed_upper = {
        "HTML", "JSON", "HTTP", "HTTPS", "UTF", "ASCII", "TRUE", "FALSE", "NONE", "NULL",
        "AVOID", "CAUTION", "MONITOR", "FAVOUR", "DOSE", "DOWN", "HIGH", "LOW", "INFO",
        "WARN", "ERROR", "CRITICAL", "PASS", "FAIL", "TODO", "NOTE", "TIP", "WARNING",
        "IMPORTANT", "HEAD", "BODY", "DIV", "SPAN", "TABLE", "THEAD", "TBODY", "TR", "TD",
        "TH", "UL", "LI", "OL", "BR", "HR", "IMG", "SVG", "PATH", "RECT", "CIRCLE", "TEXT",
        "LINE", "POLYLINE", "POLYGON", "STOP", "SPLICE", "START", "LOST", "TIER", "VEP",
        "WGS", "BAM", "VCF", "CADD", "REVEL", "AVI", "CPIC", "DPWG", "FDA", "NCCN", "ACMG",
        "AMP", "GWAS", "EHR", "PCP", "ENT", "MRI", "ECG", "CBC", "CMP", "AST", "ALT", "INR",
        "VTE", "APC", "MMA", "MCV", "SPF", "DNA", "RNA", "ATP", "ADP", "AMP", "NAD", "FAD", "MAF",
        "PROBAND", "SAMPLE", "PATIENT", "REPORT", "BRIEF", "DOSSIER", "PORTRAIT", "LANDSCAPE",
    }

    if os.path.exists(reports_dir):
        for json_f in glob.glob(os.path.join(reports_dir, "**/*_ontology_pharma_alphagenome.json"), recursive=True):
            bname = os.path.basename(json_f)
            pid = bname.split("_ontology_pharma_alphagenome.json")[0]
            for part in re.split(r"[-_.\s]+", pid):
                if len(part) >= 4:
                    pii_tokens.add(re.escape(part))
            try:
                with open(json_f, encoding="utf-8") as f:
                    data = json.load(f)
                for g in data.get("genes", []):
                    sym = str(g.get("symbol") or "").strip()
                    if sym and len(sym) >= 3 and sym.isupper() and sym not in allowed_upper:
                        gene_tokens.add(re.escape(sym))
                raw_pgx = (data.get("raw_pharmacogenomics") or {}).get("pharmacogenes") or {}
                for pg_sym in raw_pgx:
                    if pg_sym and len(pg_sym) >= 3 and pg_sym not in allowed_upper:
                        gene_tokens.add(re.escape(pg_sym))
            except Exception:
                pass

    if pii_tokens:
        static_patterns["pii_literal"] = re.compile(r"(?i)\b(?:" + "|".join(sorted(pii_tokens)) + r")\b")
    if gene_tokens:
        static_patterns["gene_symbol_literal"] = re.compile(r"\b(?:" + "|".join(sorted(gene_tokens)) + r")\b")

    return static_patterns


def audit_codebase_and_skills(root_dir=None, quiet=False):
    """
    Scan all Python files in root_dir and all skill files in ~/.gemini/config/skills.
    Returns (passed: bool, violations: list of (file_basename, line_no, category)).
    Never includes the matched PII or gene string in output.
    """
    root_dir = root_dir or SCRIPT_DIR
    reports_dir = os.path.join(root_dir, "reports")
    patterns = _build_dynamic_forbidden_patterns(reports_dir)

    target_files = sorted(
        glob.glob(os.path.join(root_dir, "*.py"))
        + glob.glob(os.path.join(root_dir, "lib/*.py"))
        + glob.glob(os.path.join(root_dir, "scripts/*.py"))
        + glob.glob(os.path.join(HOME_DIR, ".gemini/config/skills/**/*"), recursive=True)
    )
    target_files = [f for f in target_files if os.path.isfile(f) and not f.endswith(".pyc")]

    violations = []
    for fpath in target_files:
        try:
            with open(fpath, encoding="utf-8") as fp:
                lines = fp.readlines()
        except Exception:
            continue
        for line_no, line in enumerate(lines, 1):
            for cat, pat in patterns.items():
                if pat.search(line):
                    violations.append((os.path.basename(fpath), line_no, cat))

    if violations and not quiet:
        sys.stderr.write(f"[ZERO-CURATION GUARD BLOCKED] Found {len(violations)} static identifier violation(s):\n")
        for bname, lno, cat in violations[:25]:
            sys.stderr.write(f"  - {bname}:L{lno} [{cat}]\n")
    return len(violations) == 0, violations


def enforce_or_exit(root_dir=None):
    passed, _ = audit_codebase_and_skills(root_dir=root_dir, quiet=False)
    if not passed:
        sys.exit(1)


if __name__ == "__main__":
    passed, violations = audit_codebase_and_skills(quiet=False)
    if passed:
        print("ZERO_CURATION_GUARD_PASSED")
        sys.exit(0)
    sys.exit(1)
