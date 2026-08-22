#!/usr/bin/env python3
"""
demo_autoimmune.py
Build a realistic mock OpenCRAVAT-style SQLite for the AUTOIMMUNITY domain and
run the full downstream pipeline (schema probe -> filter -> enrich -> render)
WITHOUT needing OpenCRAVAT, micromamba, or the HPO/GO module databases.

Why this exists: the panel builder needs the installed hpo/go annotator SQLite
files, which aren't present in a bare dev/CI box. Here we synthesise a panel
directly from the config's force_include genes (plus the demo genes) so the rest
of the pipeline — the parts we actually changed — can be exercised and a genuine
sample report produced.

Usage:
  python3 tests/demo_autoimmune.py OUTDIR [--offline]

Outputs into OUTDIR: DEMO_autoimmunity_report.{html,tsv,txt}, the actionable
JSON, and the enrichment cache. With network access it fills in live NCBI gene
descriptions and current GWAS Catalog study evidence; --offline uses cache only.
"""
import argparse
import json
import os
import sqlite3
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "lib"))

import ontology_filter as of        # noqa: E402
import schema_probe as sp           # noqa: E402
import enrich_report as er          # noqa: E402
import render_autoimmune as ra      # noqa: E402

try:
    import yaml
except ImportError:
    sys.exit("PyYAML required")

CONFIG = os.path.join(ROOT, "config", "autoimmunity.yaml")

# Variant table columns for the mock DB (a realistic subset).
VCOLS = [
    "base__uid", "base__hugo", "base__so", "base__coding", "base__achange",
    "base__cchange", "base__transcript", "base__chrom", "base__pos",
    "base__ref_base", "base__alt_base",
    "dbsnp__rsid",
    "vcfinfo__zygosity", "vcfinfo__alt_reads", "vcfinfo__tot_reads", "vcfinfo__af",
    "gnomad4__af", "allofus250k__gvs_all_af",
    "clinvar__sig", "clinvar__id", "clinvar__disease_names", "clinvar__rev_stat",
    "alphamissense__am_pathogenicity", "alphamissense__am_class", "revel__score",
    "metarnn__score", "esm1b__score", "varity_r__varity_r",
    "bayesdel__bayesdel_addAF_score",
    "spliceai__ds_ag", "spliceai__ds_al", "spliceai__ds_dg", "spliceai__ds_dl",
    "cadd__phred", "linsight__value", "ncer__score", "regulomedb__ra",
    # GWAS catalog annotator columns (autoimmune-risk evidence)
    "gwas_catalog__disease", "gwas_catalog__trait", "gwas_catalog__pmid",
    "gwas_catalog__risk_allele",
]
GCOLS = ["base__hugo", "hpo__id", "hpo__term", "go__bpo_name", "go__mfo_name",
         "go__cco_name", "clingen__classification", "clingen__disease"]


def V(uid, hugo, so, **kw):
    row = {c: None for c in VCOLS}
    row.update({
        "base__uid": uid, "base__hugo": hugo, "base__so": so,
        "base__chrom": kw.get("chrom", "chr1"), "base__pos": 1000 + uid,
        "base__ref_base": kw.get("ref", "A"), "base__alt_base": kw.get("alt", "G"),
    })
    alias = {
        "rsid": "dbsnp__rsid", "zyg": "vcfinfo__zygosity",
        "alt_reads": "vcfinfo__alt_reads", "tot_reads": "vcfinfo__tot_reads",
        "vaf": "vcfinfo__af", "af": "gnomad4__af", "aou": "allofus250k__gvs_all_af",
        "clinvar": "clinvar__sig", "clinvar_id": "clinvar__id",
        "am": "alphamissense__am_pathogenicity", "revel": "revel__score",
        "metarnn": "metarnn__score", "esm1b": "esm1b__score",
        "varity": "varity_r__varity_r", "achange": "base__achange",
        "gwas_disease": "gwas_catalog__disease", "gwas_trait": "gwas_catalog__trait",
        "gwas_pmid": "gwas_catalog__pmid", "gwas_risk": "gwas_catalog__risk_allele",
        "ds_ag": "spliceai__ds_ag", "cadd": "cadd__phred",
    }
    for k, v in kw.items():
        if k in ("chrom", "ref", "alt"):
            continue
        row[alias.get(k, k)] = v
    return [row[c] for c in VCOLS]


def build_db(path):
    if os.path.exists(path):
        os.remove(path)
    conn = sqlite3.connect(path)
    cur = conn.cursor()
    cur.execute("CREATE TABLE variant (%s)" % ", ".join('"%s"' % c for c in VCOLS))
    cur.execute("CREATE TABLE gene (%s)" % ", ".join('"%s"' % c for c in GCOLS))

    variants = [
        # --- catalogued common GWAS risk alleles (kept via gwas_catalog bypass) ---
        V(1, "PTPN22", "MIS", rsid="rs2476601", af=0.09, zyg="het",
          alt_reads=54, tot_reads=110, achange="p.R620W", revel=0.6,
          gwas_disease="Type 1 diabetes", gwas_trait="type 1 diabetes mellitus",
          gwas_pmid="17554260", gwas_risk="rs2476601-T", chrom="chr1"),
        V(2, "CTLA4", "INT", rsid="rs3087243", af=0.45, zyg="hom",
          gwas_disease="Autoimmune thyroid disease", gwas_trait="Graves disease",
          gwas_pmid="12724780", gwas_risk="rs3087243-G", chrom="chr2"),
        V(3, "IL23R", "MIS", rsid="rs11209026", af=0.07, zyg="het",
          achange="p.R381Q", gwas_disease="Inflammatory bowel disease",
          gwas_trait="inflammatory bowel disease", gwas_pmid="17068223",
          gwas_risk="rs11209026-A", chrom="chr1"),
        V(4, "STAT4", "INT", rsid="rs7574865", af=0.22, zyg="het",
          gwas_disease="Systemic lupus erythematosus", gwas_trait="systemic lupus erythematosus",
          gwas_pmid="17804842", gwas_risk="rs7574865-T", chrom="chr2"),
        V(5, "TNFAIP3", "INT", rsid="rs2230926", af=0.03, zyg="het",
          gwas_disease="Rheumatoid arthritis", gwas_trait="rheumatoid arthritis",
          gwas_pmid="18345021", gwas_risk="rs2230926-G", chrom="chr6"),
        V(6, "IL7R", "MIS", rsid="rs6897932", af=0.25, zyg="hom",
          achange="p.T244I", gwas_disease="Multiple sclerosis", gwas_trait="multiple sclerosis",
          gwas_pmid="17660530", gwas_risk="rs6897932-C", chrom="chr5"),
        V(7, "HLA-DRB1", "INT", rsid="rs660895", af=0.18, zyg="het",
          gwas_disease="Rheumatoid arthritis", gwas_trait="rheumatoid arthritis",
          gwas_pmid="22446963", gwas_risk="rs660895-G", chrom="chr6"),
        V(8, "TYK2", "MIS", rsid="rs34536443", af=0.04, zyg="het",
          achange="p.P1104A", revel=0.7, gwas_disease="Multiple autoimmune",
          gwas_trait="rheumatoid arthritis", gwas_pmid="26301688",
          gwas_risk="rs34536443-C", chrom="chr19"),

        # --- monogenic / high-impact (kept via ClinVar / LoF / predictors) ---
        V(9, "AIRE", "MIS", rsid="rs1800522", af=1e-4, zyg="hom",
          clinvar="Pathogenic", clinvar_id="3308", achange="p.R257X",
          revel=0.9, am=0.95, chrom="chr21"),
        V(10, "FOXP3", "STG", rsid=None, af=None, zyg="hemi",
          achange="p.Q210*", chrom="chrX"),
        V(11, "FAS", "SPL", rsid=None, af=1e-5, zyg="het", ds_ag=0.85,
          achange="c.1-2A>G", chrom="chr10"),
        V(12, "LRBA", "MIS", rsid=None, af=1e-5, zyg="het",
          revel=0.9, am=0.9, metarnn=0.9, esm1b=-9.0, varity=0.9,
          achange="p.G1234R", chrom="chr4"),

        # --- should be DROPPED ---
        # 13: common, no GWAS catalog entry, no coding/clinical signal
        V(13, "STAT4", "INT", rsid="rs99999999", af=0.30, zyg="het", chrom="chr2"),
        # 14: not in panel at all
        V(14, "BRCA1", "STG", af=None, zyg="het", chrom="chr17"),
    ]
    cur.executemany("INSERT INTO variant VALUES (%s)" % ",".join("?" * len(VCOLS)),
                    variants)

    genes = [
        ("PTPN22", "HP:0002960", "Autoimmunity", "", "T cell activation", "", None, None),
        ("CTLA4", "HP:0002960", "Autoimmunity", "regulation of immune response", "", "", "Definitive", "Immune dysregulation"),
        ("IL23R", "HP:0002715", "Abnormality of the immune system", "inflammatory response", "interleukin-23 receptor activity", "", None, None),
        ("STAT4", "HP:0002960", "Autoimmunity", "cytokine-mediated signaling pathway", "", "", None, None),
        ("TNFAIP3", "HP:0002960", "Autoimmunity", "regulation of immune response", "", "", None, None),
        ("IL7R", "HP:0002715", "Abnormality of the immune system", "T cell differentiation", "", "", None, None),
        ("HLA-DRB1", "HP:0002960", "Autoimmunity", "antigen processing and presentation", "MHC class II protein complex", "", None, None),
        ("TYK2", "HP:0002960", "Autoimmunity", "cytokine-mediated signaling pathway", "", "", None, None),
        ("AIRE", "HP:0002715", "Abnormality of the immune system", "tolerance induction", "", "", "Definitive", "APECED"),
        ("FOXP3", "HP:0002960", "Autoimmunity", "regulatory T cell differentiation", "", "", "Definitive", "IPEX"),
        ("FAS", "HP:0002960", "Autoimmunity", "regulation of apoptotic process", "", "", "Definitive", "ALPS"),
        ("LRBA", "HP:0002715", "Abnormality of the immune system", "regulation of immune response", "", "", "Definitive", "LRBA deficiency"),
    ]
    cur.executemany("INSERT INTO gene VALUES (%s)" % ",".join("?" * len(GCOLS)), genes)
    conn.commit()
    conn.close()


def build_panel(path):
    """Synthesise a panel.json (config force_include + demo genes) so we don't
    need the installed HPO/GO databases in this environment."""
    cfg = yaml.safe_load(open(CONFIG))
    genes = set(cfg.get("panel", {}).get("force_include", []) or [])
    genes |= {"PTPN22", "CTLA4", "IL23R", "STAT4", "TNFAIP3", "IL7R",
              "HLA-DRB1", "TYK2", "AIRE", "FOXP3", "FAS", "LRBA"}
    panel = {"config_domain": "autoimmunity",
             "genes": {g: {"support": 2, "forced": True, "hpo": ["x"], "go": ["y"]}
                       for g in sorted(genes)},
             "counts": {"total_genes": len(genes)},
             "min_ontology_support": 1}
    json.dump(panel, open(path, "w"), indent=2)
    return panel


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("outdir")
    ap.add_argument("--offline", action="store_true")
    ap.add_argument("--patient", default="DEMO")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    db = os.path.join(args.outdir, "DEMO.sqlite")
    panel_p = os.path.join(args.outdir, "DEMO_autoimmunity_panel.json")
    schema_p = os.path.join(args.outdir, "DEMO_schema.json")
    act_sqlite = os.path.join(args.outdir, "DEMO_autoimmunity_actionable.sqlite")
    act_json = os.path.join(args.outdir, "DEMO_autoimmunity_actionable.json")
    cache = os.path.join(args.outdir, "DEMO_autoimmunity_enrich_cache.json")
    html_p = os.path.join(args.outdir, "DEMO_autoimmunity_report.html")
    tsv_p = os.path.join(args.outdir, "DEMO_autoimmunity_report.tsv")
    text_p = os.path.join(args.outdir, "DEMO_autoimmunity_report.txt")

    build_db(db)
    build_panel(panel_p)
    schema = sp.probe(db)
    json.dump(schema, open(schema_p, "w"), indent=2)

    of.run(db, panel_p, schema_p, CONFIG, act_sqlite, act_json, args.patient)

    data = json.load(open(act_json))
    data, meta = er.enrich(data, cache, do_genes=True, do_studies=True,
                           offline=args.offline)
    json.dump(data, open(act_json, "w"), indent=2)

    ra.write_html(data, html_p)
    ra.rr.write_tsv(data["records"], tsv_p)
    ra.rr.write_text(data, text_p)

    print(f"\n[demo] enrichment: {meta}")
    print(f"[demo] report -> {html_p}")


if __name__ == "__main__":
    main()
