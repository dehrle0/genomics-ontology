#!/usr/bin/env python3
"""
schema_probe.py
Inspect an OpenCRAVAT result SQLite and report which relevant columns exist in
the `variant` and `gene` tables. Column names vary by annotator version, so the
downstream filter references only columns that are actually present.
"""
import argparse
import json
import sqlite3
import sys


def table_columns(cur, table):
    try:
        return [r[1] for r in cur.execute(f"PRAGMA table_info({table})")]
    except sqlite3.OperationalError:
        return []


def has(cols, name):
    return name if name in cols else None


def first_present(cols, candidates):
    for c in candidates:
        if c in cols:
            return c
    return None


def probe(db_path):
    conn = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    cur = conn.cursor()
    vcols = table_columns(cur, "variant")
    gcols = table_columns(cur, "gene")
    scols = table_columns(cur, "sample")  # per-sample genotype (zygosity/reads)
    vg = vcols + gcols  # some annotators land in either table by version
    tables = [r[0] for r in cur.execute("SELECT name FROM sqlite_master WHERE type='table'")]

    schema = {
        "db": db_path,
        "tables": tables,
        "variant_col_count": len(vcols),
        "gene_col_count": len(gcols),
        "sample_col_count": len(scols),
        # --- identity / consequence ---
        "uid": has(vcols, "base__uid"),
        "hugo": has(vcols, "base__hugo"),
        "so": has(vcols, "base__so"),
        "coding": has(vcols, "base__coding"),
        "achange": has(vcols, "base__achange"),
        "cchange": has(vcols, "base__cchange"),
        "transcript": has(vcols, "base__transcript"),
        "chrom": has(vcols, "base__chrom"),
        "pos": has(vcols, "base__pos"),
        "ref": has(vcols, "base__ref_base"),
        "alt": has(vcols, "base__alt_base"),
        # --- dbSNP rsID (used to fetch live study evidence) ---
        "rsid": first_present(vcols, ["dbsnp__rsid", "dbsnp__snp", "base__dbsnp"]),
        # --- genotype / zygosity ---------------------------------------------
        # Single-sample OC databases usually carry genotype on the variant table
        # via the vcfinfo annotator; multi-sample databases carry it per-sample
        # in the `sample` table. We record whichever is present. The filter
        # prefers the variant-level columns and falls back to the sample table.
        "zygosity": first_present(vcols, ["vcfinfo__zygosity"]),
        "alt_reads": first_present(vcols, ["vcfinfo__alt_reads"]),
        "tot_reads": first_present(vcols, ["vcfinfo__tot_reads"]),
        "vaf": first_present(vcols, ["vcfinfo__af"]),
        "sample_uid": has(scols, "base__uid"),
        "sample_zygosity": has(scols, "base__zygosity"),
        "sample_alt_reads": has(scols, "base__alt_reads"),
        "sample_tot_reads": has(scols, "base__tot_reads"),
        "sample_vaf": has(scols, "base__af"),
        # --- Phasing / Haplotype blocks ---
        "hap_block": first_present(vcols, ["vcfinfo__hap_block", "vcfinfo__ps", "vcfinfo__phase_set"]),
        "hap_strand": first_present(vcols, ["vcfinfo__hap_strand", "vcfinfo__strand"]),
        "sample_hap_block": has(scols, "base__hap_block"),
        "sample_hap_strand": has(scols, "base__hap_strand"),
        # --- population ---
        "gnomad4_af": has(vcols, "gnomad4__af"),
        "allofus_af": first_present(vcols, ["allofus250k__gvs_all_af", "allofus250k__af"]),
        # --- clinical / phenotype ---
        "clinvar_sig": has(vcols, "clinvar__sig"),
        "clinvar_id": has(vcols, "clinvar__id"),
        "clinvar_disease": first_present(vcols, ["clinvar__disease_names", "clinvar__disease_name"]),
        "clinvar_rev": has(vcols, "clinvar__rev_stat"),
        # clingen is a gene-level annotator -> columns live in the gene table
        "clingen_class": first_present(gcols, ["clingen__classification"]),
        "clingen_disease": first_present(gcols, ["clingen__disease"]),
        "omim_id": first_present(vcols, ["omim__omim_id", "omim__id"]),
        "clinvar_acmg_ps1": first_present(vcols, ["clinvar_acmg__ps1", "clinvar_acmg__ps1_id"]),
        "clinvar_acmg_pm5": first_present(vcols, ["clinvar_acmg__pm5", "clinvar_acmg__pm5_id"]),
        # --- coding predictors (generic, pan-disease) ---
        "revel": first_present(vcols, ["revel__score"]),
        "am_path": has(vcols, "alphamissense__am_pathogenicity"),
        "am_class": has(vcols, "alphamissense__am_class"),
        "bayesdel": first_present(vcols, ["bayesdel__bayesdel_addAF_score", "bayesdel__bayesdel_noAF_score"]),
        "metarnn": first_present(vcols, ["metarnn__score"]),
        "esm1b": first_present(vcols, ["esm1b__score"]),
        "varity": first_present(vcols, ["varity_r__varity_r"]),
        # --- splice ---
        "spliceai_ds_ag": has(vcols, "spliceai__ds_ag"),
        "spliceai_ds_al": has(vcols, "spliceai__ds_al"),
        "spliceai_ds_dg": has(vcols, "spliceai__ds_dg"),
        "spliceai_ds_dl": has(vcols, "spliceai__ds_dl"),
        # --- non-coding regulatory ---
        "cadd_phred": has(vcols, "cadd__phred"),
        "linsight": has(vcols, "linsight__value"),
        "ncer": has(vcols, "ncer__score"),
        "regulomedb_ra": has(vcols, "regulomedb__ra"),
        "ccre_group": has(vcols, "ccre_screen___group"),
        # --- gene table ontology (may be joined from gene level) ---
        "gene_hpo_id": has(gcols, "hpo__id"),
        "gene_hpo_term": has(gcols, "hpo__term"),
        "gene_go_bpo": has(gcols, "go__bpo_name"),
        "gene_go_mfo": has(gcols, "go__mfo_name"),
        "gene_go_cco": has(gcols, "go__cco_name"),
        # --- variant-level ontology (present only if --nogenelevelonvariantlevel off) ---
        "var_hpo_term": has(vcols, "hpo__term"),
        "var_go_bpo": has(vcols, "go__bpo_name"),
        # --- Structural Variants & CNVs (SURVIVOR, CNV, Manta) ---
        "svtype": first_present(vcols, ["vcfinfo__svtype", "survivor__svtype", "cnv__svtype", "manta__svtype", "base__type"]),
        "svlen": first_present(vcols, ["vcfinfo__svlen", "survivor__svlen", "cnv__svlen", "manta__svlen"]),
        "cnv_copy_number": first_present(vcols, ["cnv__copy_number", "cnv__cn", "vcfinfo__cn"]),
        # --- Deep in silico predictors (EVE, PrimateAI, GERP, phyloP) ---
        "eve_score": first_present(vcols, ["eve__score", "eve__eve_score", "eve__pathogenicity"]),
        "primateai_score": first_present(vcols, ["primateai__score", "primateai__primateai_score"]),
        "gerp_score": first_present(vcols, ["gerp__score", "gerp__gerp_rs"]),
        "phylop_score": first_present(vcols, ["phylop__score", "phyloP__score"]),
        # literature / tissue
        "pubmed_n": first_present(vcols, ["pubmed__n"]),
        "gtex_tissue": has(vcols, "gtex__gtex_tissue"),
    }
    conn.close()
    return schema


def main():
    ap = argparse.ArgumentParser(description="Probe OpenCRAVAT sqlite schema")
    ap.add_argument("db")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    schema = probe(args.db)
    with open(args.out, "w") as f:
        json.dump(schema, f, indent=2)
    present = sum(1 for k, v in schema.items() if isinstance(v, str) and v and k not in ("db",))
    print(f"[schema] {args.db}")
    print(f"[schema] variant cols={schema['variant_col_count']} "
          f"gene cols={schema['gene_col_count']} sample cols={schema['sample_col_count']}")
    missing = [k for k, v in schema.items() if v is None]
    if missing:
        print(f"[schema] absent fields: {', '.join(missing)}")
    print(f"[schema] wrote {args.out}")


if __name__ == "__main__":
    main()
