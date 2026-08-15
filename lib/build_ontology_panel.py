#!/usr/bin/env python3
"""
build_ontology_panel.py
Derive a gene panel from HPO + GO ontology annotations (OpenCRAVAT modules),
driven entirely by the domain config (config/<domain>.yaml). No hard-coded gene list.

Output: panel.json
{
  "genes": {"MYH7": {"support": 2, "hpo": [...], "go": [...], "forced": false}, ...},
  "hpo_terms_used": {...}, "go_terms_used": {...},
  "counts": {...}, "config_domain": "cardiology"
}
"""
import argparse
import json
import os
import sqlite3
import sys

try:
    import yaml
except ImportError:
    sys.exit("[panel] PyYAML required (micromamba activate cravat_env)")

DEFAULT_MODULES_DIR = "/data/opencravat/modules"


def modules_dir():
    return os.environ.get("OC_MODULES_DIR", DEFAULT_MODULES_DIR)


def hpo_db():
    return os.path.join(modules_dir(), "annotators", "hpo", "data", "hpo.sqlite")


def go_db():
    return os.path.join(modules_dir(), "annotators", "go", "data", "go.sqlite")


def _genes_for_hpo(cur, term_ids, keywords):
    """Return {gene: set(matched_term_names)} for HPO seeds."""
    hits = {}
    for tid in term_ids:
        for hugo, term in cur.execute(
            "SELECT DISTINCT hugo, hpo_term FROM genes WHERE hpo_id=?", (tid,)
        ):
            if hugo and hugo != "-":
                hits.setdefault(hugo, set()).add(f"{tid}:{term}")
    for kw in keywords:
        like = f"%{kw.lower()}%"
        for hugo, term in cur.execute(
            "SELECT DISTINCT hugo, hpo_term FROM genes WHERE lower(hpo_term) LIKE ?",
            (like,),
        ):
            if hugo and hugo != "-":
                hits.setdefault(hugo, set()).add(f"kw[{kw}]:{term}")
    return hits


def _genes_for_go(cur, term_ids, keywords, drop_iea_only):
    """Return {gene: set(matched_term_names)} for GO seeds."""
    hits = {}
    # Map GO id -> name for readable output.
    id2name = {r[0]: r[1] for r in cur.execute("SELECT go_id, name FROM go_name")}

    def _keep(evidence_code):
        if not drop_iea_only:
            return True
        return (evidence_code or "").upper() != "IEA"

    for gid in term_ids:
        for hugo, ev in cur.execute(
            "SELECT DISTINCT hugo, evidence_code FROM go_annotation WHERE go_id=?", (gid,)
        ):
            if hugo and hugo != "-" and _keep(ev):
                hits.setdefault(hugo, set()).add(f"{gid}:{id2name.get(gid, gid)}")
    if keywords:
        matched_ids = {
            gid
            for gid, name in id2name.items()
            if any(kw.lower() in name.lower() for kw in keywords)
        }
        if matched_ids:
            qmarks = ",".join("?" * len(matched_ids))
            for hugo, gid, ev in cur.execute(
                f"SELECT DISTINCT hugo, go_id, evidence_code FROM go_annotation "
                f"WHERE go_id IN ({qmarks})",
                tuple(matched_ids),
            ):
                if hugo and hugo != "-" and _keep(ev):
                    hits.setdefault(hugo, set()).add(f"{gid}:{id2name.get(gid, gid)}")
    return hits


def build_panel(config_path, out_path):
    with open(config_path) as f:
        cfg = yaml.safe_load(f)

    hpo_cfg = cfg.get("hpo", {}) or {}
    go_cfg = cfg.get("go", {}) or {}
    panel_cfg = cfg.get("panel", {}) or {}

    if not os.path.exists(hpo_db()):
        sys.exit(f"[panel] HPO data not found: {hpo_db()}")
    if not os.path.exists(go_db()):
        sys.exit(f"[panel] GO data not found: {go_db()}")

    hcur = sqlite3.connect(f"file:{hpo_db()}?mode=ro", uri=True).cursor()
    gcur = sqlite3.connect(f"file:{go_db()}?mode=ro", uri=True).cursor()

    hpo_hits = _genes_for_hpo(
        hcur, hpo_cfg.get("term_ids", []) or [], hpo_cfg.get("term_keywords", []) or []
    )
    go_hits = _genes_for_go(
        gcur,
        go_cfg.get("term_ids", []) or [],
        go_cfg.get("term_keywords", []) or [],
        bool(go_cfg.get("drop_iea_only", False)),
    )

    min_support = int(panel_cfg.get("min_ontology_support", 1))
    force_include = set(panel_cfg.get("force_include", []) or [])
    force_exclude = set(panel_cfg.get("force_exclude", []) or [])

    genes = {}
    all_genes = set(hpo_hits) | set(go_hits) | force_include
    for g in sorted(all_genes):
        if g in force_exclude:
            continue
        hpo_terms = sorted(hpo_hits.get(g, set()))
        go_terms = sorted(go_hits.get(g, set()))
        support = (1 if hpo_terms else 0) + (1 if go_terms else 0)
        forced = g in force_include
        if not forced and support < min_support:
            continue
        genes[g] = {
            "support": support,
            "forced": forced,
            "hpo": hpo_terms,
            "go": go_terms,
        }

    panel = {
        "config_domain": cfg.get("domain", "unknown"),
        "genes": genes,
        "counts": {
            "total_genes": len(genes),
            "hpo_only": sum(1 for v in genes.values() if v["hpo"] and not v["go"]),
            "go_only": sum(1 for v in genes.values() if v["go"] and not v["hpo"]),
            "both": sum(1 for v in genes.values() if v["hpo"] and v["go"]),
            "forced": sum(1 for v in genes.values() if v["forced"]),
        },
        "min_ontology_support": min_support,
    }

    with open(out_path, "w") as f:
        json.dump(panel, f, indent=2)
    return panel


def main():
    ap = argparse.ArgumentParser(description="Build ontology-driven gene panel")
    ap.add_argument("--config", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    panel = build_panel(args.config, args.out)
    c = panel["counts"]
    print(f"[panel] domain={panel['config_domain']} min_support={panel['min_ontology_support']}")
    print(
        f"[panel] genes={c['total_genes']} "
        f"(hpo_only={c['hpo_only']}, go_only={c['go_only']}, both={c['both']}, forced={c['forced']})"
    )
    print(f"[panel] wrote {args.out}")


if __name__ == "__main__":
    main()
