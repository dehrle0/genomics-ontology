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

    if not os.path.exists(hpo_db()):
        sys.exit(f"[panel] HPO data not found: {hpo_db()}")
    if not os.path.exists(go_db()):
        sys.exit(f"[panel] GO data not found: {go_db()}")

    hcur = sqlite3.connect(f"file:{hpo_db()}?mode=ro", uri=True).cursor()
    gcur = sqlite3.connect(f"file:{go_db()}?mode=ro", uri=True).cursor()

    # Determine if single-domain or multi-domain registry
    domain_sections = []
    if "hpo" in cfg or "go" in cfg or "panel" in cfg:
        domain_sections.append(cfg)
    else:
        systems = cfg.get("level1_systems", cfg)
        for k, v in systems.items():
            if isinstance(v, dict):
                domain_sections.append(v)
                for l2k, l2v in (v.get("level2_subcategories") or {}).items():
                    if isinstance(l2v, dict):
                        domain_sections.append(l2v)

    hpo_hits = {}
    go_hits = {}
    force_include = set()
    force_exclude = set()
    min_support = 1
    config_domain = cfg.get("domain", "all_domains") if "domain" in cfg else "master_hub"

    for sec in domain_sections:
        hpo_cfg = sec.get("hpo", {}) or {}
        go_cfg = sec.get("go", {}) or {}
        panel_cfg = sec.get("panel", {}) or {}

        h_terms = (hpo_cfg.get("term_ids", []) or []) + (sec.get("hpo_terms", []) or [])
        h_kws = (hpo_cfg.get("term_keywords", []) or []) + (sec.get("keywords", []) or [])

        g_terms = (go_cfg.get("term_ids", []) or []) + (sec.get("go_terms", []) or [])
        g_kws = (go_cfg.get("term_keywords", []) or [])

        h_hits = _genes_for_hpo(hcur, h_terms, h_kws)
        g_hits = _genes_for_go(
            gcur,
            g_terms,
            g_kws,
            bool(go_cfg.get("drop_iea_only", False)),
        )

        for g, terms in h_hits.items():
            hpo_hits.setdefault(g, set()).update(terms)
        for g, terms in g_hits.items():
            go_hits.setdefault(g, set()).update(terms)

        force_include.update(panel_cfg.get("force_include", []) or [])
        force_exclude.update(panel_cfg.get("force_exclude", []) or [])
        if "min_ontology_support" in panel_cfg:
            min_support = int(panel_cfg["min_ontology_support"])

    genes = {}
    all_genes = (set(hpo_hits) | set(go_hits) | force_include) - force_exclude
    for g in sorted(all_genes):
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
        "config_domain": config_domain,
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
