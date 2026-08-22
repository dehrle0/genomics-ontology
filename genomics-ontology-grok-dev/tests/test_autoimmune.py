#!/usr/bin/env python3
"""
Offline validation for the autoimmunity domain + the new engine features
(zygosity, gene-description enrichment, live-study enrichment, autoimmune
renderer). Fully offline and deterministic: the enrichment cache is pre-seeded
so no network is required and results are reproducible. Runs in <1s.

Usage:  python3 tests/test_autoimmune.py
"""
import json
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "lib"))
sys.path.insert(0, HERE)

import ontology_filter as of        # noqa: E402
import schema_probe as sp           # noqa: E402
import enrich_report as er          # noqa: E402
import render_autoimmune as ra      # noqa: E402
import demo_autoimmune as demo       # noqa: E402

CONFIG = os.path.join(ROOT, "config", "autoimmunity.yaml")


def expect(cond, msg, failures):
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        failures.append(msg)


def seed_cache(path):
    """Pre-seed the enrichment cache so the test never touches the network."""
    cache = {
        "_version": er.CACHE_VERSION,
        # gene descriptions (final memoised form used by gene_info)
        "gene:PTPN22": {"ok": {
            "ncbi_gene_id": "26191",
            "description": "protein tyrosine phosphatase non-receptor type 22",
            "summary": "Lymphoid-specific intracellular phosphatase; risk alleles "
                       "associate with multiple autoimmune diseases.",
            "map_location": "1p13.2", "aliases": "LYP PEP PTPN8"}},
        # per-variant GWAS studies (final memoised form used by snp_studies)
        "gwas_snp:rs2476601": {"ok": {
            "rsid": "rs2476601", "n_associations": 163, "top": [
                {"traits": ["type 1 diabetes mellitus"], "pvalue": "2e-80",
                 "or_beta": 1.9, "risk_allele": "rs2476601-T", "risk_freq": "0.09",
                 "pubmed": {"pubmed_id": "17554260", "title": "T1D GWAS",
                            "journal": "Nat Genet", "date": "2007-06-06",
                            "author": "Todd JA"}},
                {"traits": ["rheumatoid arthritis"], "pvalue": "5e-40",
                 "or_beta": 1.5, "risk_allele": "rs2476601-A", "risk_freq": None,
                 "pubmed": {"pubmed_id": "20453842", "title": "RA GWAS",
                            "journal": "Nat Genet", "date": "2010", "author": "Stahl EA"}},
            ]}},
    }
    json.dump(cache, open(path, "w"))


def main():
    failures = []
    tmp = tempfile.mkdtemp(prefix="autotest_")
    db = os.path.join(tmp, "mock.sqlite")
    panel_p = os.path.join(tmp, "panel.json")
    schema_p = os.path.join(tmp, "schema.json")
    act_sqlite = os.path.join(tmp, "act.sqlite")
    act_json = os.path.join(tmp, "act.json")
    cache = os.path.join(tmp, "cache.json")
    html_p = os.path.join(tmp, "report.html")

    # Reuse the demo's realistic mock builder (variants + panel).
    demo.build_db(db)
    demo.build_panel(panel_p)
    schema = sp.probe(db)
    json.dump(schema, open(schema_p, "w"))

    # --- schema probe picked up genotype + rsID columns ---
    expect(schema.get("zygosity") == "vcfinfo__zygosity",
           "schema_probe detects vcfinfo zygosity column", failures)
    expect(schema.get("rsid") == "dbsnp__rsid",
           "schema_probe detects dbSNP rsID column", failures)

    kept = of.run(db, panel_p, schema_p, CONFIG, act_sqlite, act_json, "TESTAI")
    by_gene = {}
    for k in kept:
        by_gene.setdefault(k["hugo"], []).append(k)

    # --- GWAS common risk allele bypasses the frequency ceiling ---
    ptpn22 = by_gene.get("PTPN22", [{}])[0]
    expect("GWAS_RISK_ALLELE" in ptpn22.get("reason_codes", ""),
           "common GWAS risk allele (PTPN22 rs2476601) kept via gwas_catalog", failures)
    expect("RISK_ALLELE_COMMON" in ptpn22.get("reason_codes", ""),
           "common risk allele flagged RISK_ALLELE_COMMON (not COMMON_AF_FLAG)", failures)
    expect(ptpn22.get("evidence", {}).get("zygosity") == "Heterozygous",
           "zygosity label derived for PTPN22 (Heterozygous)", failures)
    expect(ptpn22.get("evidence", {}).get("vaf") == round(54 / 110, 4),
           "variant allele fraction computed from read depths", failures)

    # --- zygosity normalisation across het/hom/hemi ---
    expect(by_gene.get("CTLA4", [{}])[0].get("evidence", {}).get("zygosity") == "Homozygous",
           "hom -> Homozygous", failures)
    expect(by_gene.get("FOXP3", [{}])[0].get("evidence", {}).get("zygosity") == "Hemizygous",
           "hemi (chrX) -> Hemizygous", failures)

    # --- monogenic high-impact still tiered correctly ---
    expect(by_gene.get("AIRE", [{}])[0].get("tier") == "Tier1",
           "AIRE ClinVar Pathogenic rare -> Tier1", failures)
    expect("PVS1_HAPLOINSUFFICIENT" in by_gene.get("FOXP3", [{}])[0].get("reason_codes", ""),
           "FOXP3 LoF in haploinsufficient gene -> PVS1", failures)

    # --- negatives ---
    expect("BRCA1" not in by_gene, "non-panel BRCA1 dropped", failures)
    stat4 = by_gene.get("STAT4", [])
    expect(all("rs99999999" != s.get("rsid") for s in stat4),
           "common STAT4 SNP with no GWAS/clinical evidence dropped", failures)

    # --- enrichment (offline, pre-seeded cache) injects fields ---
    seed_cache(cache)
    data = json.load(open(act_json))
    data, meta = er.enrich(data, cache, do_genes=True, do_studies=True, offline=True)
    ptp = next((r for r in data["records"] if r["hugo"] == "PTPN22"), {})
    expect(ptp.get("gene_info", {}).get("description", "").startswith("protein tyrosine"),
           "NCBI gene description injected from cache (PTPN22)", failures)
    expect(bool(ptp.get("study_evidence", {}).get("snp", {}).get("top")),
           "live GWAS study evidence injected from cache (PTPN22)", failures)
    expect(meta["remote_calls"] == 0,
           "offline enrichment made zero network calls", failures)

    # --- autoimmune renderer produces the visualization + study tables ---
    ra.write_html(data, html_p)
    h = open(html_p).read()
    expect('class="trait-chart"' in h and "<rect" in h,
           "autoimmune renderer emits SVG trait-burden chart", failures)
    expect('class="study-tbl"' in h,
           "autoimmune renderer emits GWAS study-evidence table", failures)
    expect("type 1 diabetes mellitus" in h and "pubmed.ncbi.nlm.nih.gov" in h,
           "study table shows trait + PubMed citation", failures)
    expect('class="gene-summary"' in h,
           "gene NCBI summary shown on card", failures)

    # --- trait aggregation for the chart ---
    traits = ra.collect_traits(data["records"])
    tnames = {t["trait"] for t in traits}
    expect("type 1 diabetes mellitus" in tnames and "rheumatoid arthritis" in tnames,
           "trait aggregation rolls up GWAS traits across variants", failures)

    print()
    if failures:
        print(f"RESULT: {len(failures)} FAILURE(S)")
        for f in failures:
            print(f"   - {f}")
        sys.exit(1)
    print("RESULT: ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
