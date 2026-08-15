#!/usr/bin/env python3
"""
ontology_filter.py
Select ACTIONABLE variants from an OpenCRAVAT result SQLite, restricted to the
ontology-derived gene panel, then assign tier + explainable reason codes.

The engine is DOMAIN-AGNOSTIC. Everything disease-specific (which predictors,
which extra evidence columns, which LoF-sensitive genes) comes from the config
YAML. Absent annotator columns are silently ignored, so one config can run
against databases annotated with different module sets.

Input : raw OpenCRAVAT .sqlite, panel.json, schema.json, config yaml
Output: <prefix>_actionable.sqlite  (tables: variant, report_meta, panel_gene)
        + <prefix>_actionable.json   (structured records for the renderer)
"""
import argparse
import json
import os
import sqlite3
import sys

try:
    import yaml
except ImportError:
    sys.exit("[filter] PyYAML required (micromamba activate cravat_env)")

# OpenCRAVAT hg38 sequence-ontology short codes (see mappers/hg38/hg38.py)
LOF_CODES = {"STG", "FSD", "FSI", "SPL", "MLO", "STL", "EXL", "TAB",
             "stop_gained", "frameshift_variant", "frameshift_elongation",
             "frameshift_truncation", "splice_acceptor_variant",
             "splice_donor_variant", "start_lost", "stop_lost",
             "exon_loss_variant", "transcript_ablation"}
MISSENSE_CODES = {"MIS", "missense_variant"}
INFRAME_CODES = {"IND", "INI", "inframe_deletion", "inframe_insertion", "CSS", "complex_substitution"}
CODING_ALTERING = LOF_CODES | MISSENSE_CODES | INFRAME_CODES

# Logical fields whose schema column lives in the gene-level table.
GENE_TABLE_KEYS = {"gene_hpo_id", "gene_hpo_term", "gene_go_bpo", "gene_go_mfo",
                   "gene_go_cco", "clingen_class", "clingen_disease"}

# Fixed logical fields pulled for every variant (mapped via schema.json).
PULL_KEYS = [
    "uid", "hugo", "so", "coding", "achange", "cchange", "transcript",
    "chrom", "pos", "ref", "alt", "rsid",
    "zygosity", "alt_reads", "tot_reads", "vaf",
    "gnomad4_af", "allofus_af",
    "clinvar_sig", "clinvar_id", "clinvar_disease", "clinvar_rev",
    "clingen_class", "omim_id", "clinvar_acmg_ps1", "clinvar_acmg_pm5",
    "revel", "am_path", "am_class",
    "bayesdel", "metarnn", "esm1b", "varity",
    "spliceai_ds_ag", "spliceai_ds_al", "spliceai_ds_dg", "spliceai_ds_dl",
    "cadd_phred", "linsight", "ncer", "regulomedb_ra", "ccre_group",
    "gene_hpo_id", "gene_hpo_term", "gene_go_bpo", "gene_go_mfo", "gene_go_cco",
    "gtex_tissue",
]


def _num(x):
    try:
        if x is None or x == "":
            return None
        return float(x)
    except (TypeError, ValueError):
        return None


def zygosity_label(raw):
    """Normalize the many zygosity spellings OpenCRAVAT/VCF tools emit into a
    single human label. Returns None when unknown."""
    if raw is None:
        return None
    s = str(raw).strip().lower()
    if not s or s in ("-", "na", "none", "unknown", "."):
        return None
    if s in ("het", "heterozygous", "0/1", "1/0", "0|1", "1|0"):
        return "Heterozygous"
    if s in ("hom", "homozygous", "1/1", "1|1"):
        return "Homozygous"
    if s in ("hemi", "hemizygous", "1", "1/.", "./1"):
        return "Hemizygous"
    if s in ("ref", "0/0", "0|0", "homref"):
        return "Reference"
    return str(raw)


def compute_vaf(vaf, alt_reads, tot_reads):
    """Prefer an explicit VAF; otherwise derive it from allele read depths."""
    v = _num(vaf)
    if v is not None:
        return v
    a = _num(alt_reads)
    t = _num(tot_reads)
    if a is not None and t not in (None, 0):
        return round(a / t, 4)
    return None


def _clinvar_class(sig):
    if not sig:
        return None
    s = sig.lower()
    if "conflicting" in s:
        return "CONFLICT"
    if "pathogenic" in s and "likely" in s and "/" in s:
        return "PLP"
    if s.startswith("pathogenic") or s == "pathogenic":
        return "PLP"
    if "likely pathogenic" in s:
        return "PLP"
    if "uncertain" in s:
        return "VUS"
    return None


def _ontology_reasons(text, keywords):
    """Return matched keyword tokens found in an ontology term-name string."""
    if not text:
        return []
    low = text.lower()
    hits = [kw.replace(" ", "_").upper() for kw in (keywords or []) if kw.lower() in low]
    return sorted(set(hits))


def evaluate_variant(row, cfg, panel, haploinsufficient, runtime):
    """Given a dict row of pulled fields, return (keep, tier, reasons, evidence)."""
    p = cfg["predictors"]
    freq = cfg["frequency"]

    reasons = []
    pheno = []
    geno = []

    hugo = row.get("hugo")
    so = (row.get("so") or "").strip()

    gnomad = _num(row.get("gnomad4_af"))
    aou = _num(row.get("allofus_af"))

    def le(v, thr):
        return (v is None) or (v <= thr)
    is_rare_t1 = le(gnomad, freq["tier1_af"]) and le(aou, freq["tier1_af"])
    is_rare_t2 = le(gnomad, freq["tier2_af"]) and le(aou, freq["tier2_af"])
    over_ceiling = (gnomad is not None and gnomad > freq["max_af"]) or \
                   (aou is not None and aou > freq["max_af"])

    # ---------------- Phenotype evidence ----------------
    cvc = _clinvar_class(row.get("clinvar_sig"))
    if cvc == "PLP":
        pheno.append("CLINVAR_PLP")
    elif cvc == "VUS":
        pheno.append("CLINVAR_VUS")
    elif cvc == "CONFLICT":
        pheno.append("CLINVAR_CONFLICT")

    clingen = row.get("clingen_class")
    if clingen and str(clingen).lower() not in ("no known disease relationship", "", "none"):
        pheno.append("CLINGEN_VALIDITY")
    if row.get("omim_id"):
        pheno.append("OMIM_DISEASE")

    hpo_hits = _ontology_reasons(row.get("gene_hpo_term"), cfg["hpo"].get("term_keywords"))
    for h in hpo_hits:
        pheno.append(f"HPO_{h}")
    go_text = " ".join(t for t in (row.get("gene_go_bpo"), row.get("gene_go_mfo"),
                                    row.get("gene_go_cco")) if t)
    go_hits = _ontology_reasons(go_text, cfg["go"].get("term_keywords"))
    for g in go_hits:
        pheno.append(f"GO_{g}")

    # Config-driven domain evidence columns (presence -> reason)
    domain_pheno = False
    domain_bypass = False   # a firing "bypass_frequency" evidence keeps common variants
    for ev in runtime["domain_evidence"]:
        if any(row.get(a) not in (None, "") for a in ev["aliases"]):
            reason = ev["code"]
            (pheno if ev.get("kind", "phenotype") == "phenotype" else geno).append(reason)
            if ev.get("kind", "phenotype") == "phenotype":
                domain_pheno = True
            if ev.get("bypass_frequency"):
                domain_bypass = True

    # ---------------- Genotype evidence ----------------
    is_lof = so in LOF_CODES
    is_missense = so in MISSENSE_CODES
    is_coding_altering = so in CODING_ALTERING
    if is_lof:
        geno.append(f"LOF_{so}")
    elif is_missense:
        geno.append("MISSENSE")
    elif so in INFRAME_CODES:
        geno.append("INFRAME_INDEL")

    # Core (pan-disease) predictors
    pred_hits = []
    if _num(row.get("revel")) is not None and _num(row.get("revel")) >= p["revel_min"]:
        pred_hits.append("REVEL")
    am = _num(row.get("am_path"))
    am_class = (row.get("am_class") or "").lower()
    if (am is not None and am >= p["alphamissense_min"]) or "pathogenic" in am_class:
        pred_hits.append("ALPHAMISSENSE")
    bd = _num(row.get("bayesdel"))
    if bd is not None and bd >= p["bayesdel_min"]:
        pred_hits.append("BAYESDEL")
    mr = _num(row.get("metarnn"))
    if mr is not None and mr >= p["metarnn_min"]:
        pred_hits.append("METARNN")
    es = _num(row.get("esm1b"))
    if es is not None and es <= p["esm1b_max"]:
        pred_hits.append("ESM1B")
    vr = _num(row.get("varity"))
    if vr is not None and vr >= p["varity_min"]:
        pred_hits.append("VARITY")

    # Config-driven domain predictors (disease-tuned models)
    for dp in runtime["domain_predictors"]:
        hit = False
        for a in dp["score_aliases"]:
            v = _num(row.get(a))
            if v is not None and v >= dp["min"]:
                hit = True
                break
        if not hit:
            for a in dp["text_aliases"]:
                t = row.get(a)
                if t and "pathogenic" in str(t).lower():
                    hit = True
                    break
        if hit:
            pred_hits.append(dp["code"])

    for ph in pred_hits:
        geno.append(f"PP3_{ph}")
    consensus = len(pred_hits) >= p["pp3_consensus_n"]
    if consensus:
        geno.append("PP3_CONSENSUS")

    # Splice
    ds = [_num(row.get(k)) for k in ("spliceai_ds_ag", "spliceai_ds_al",
                                     "spliceai_ds_dg", "spliceai_ds_dl")]
    ds = [d for d in ds if d is not None]
    spliceai_max = max(ds) if ds else None
    if spliceai_max is not None and spliceai_max >= p["spliceai_tier1"]:
        geno.append("SPLICEAI_HIGH")
    elif spliceai_max is not None and spliceai_max >= p["spliceai_min"]:
        geno.append("SPLICEAI_MOD")

    # Rarity
    if is_rare_t1:
        geno.append("PM2_RARE")
    elif is_rare_t2:
        geno.append("RARE_TIER2")

    # Non-coding regulatory (only for non-coding-altering variants)
    nc = cfg["noncoding"]
    noncoding_hits = []
    if not is_coding_altering:
        cadd = _num(row.get("cadd_phred"))
        if cadd is not None and cadd >= nc["cadd_phred_min"]:
            noncoding_hits.append("NONCODING_CADD")
        lin = _num(row.get("linsight"))
        if lin is not None and lin >= nc["linsight_min"]:
            noncoding_hits.append("NONCODING_LINSIGHT")
        ncer = _num(row.get("ncer"))
        if ncer is not None and ncer >= nc["ncer_percentile_min"]:
            noncoding_hits.append("NONCODING_NCER")
        rdb = row.get("regulomedb_ra")
        if rdb and str(rdb).strip() and str(rdb).strip()[0] in ("1", "2", "3"):
            noncoding_hits.append("NONCODING_REGULOMEDB")
    for h in noncoding_hits:
        geno.append(h)

    # ---------------- Actionability gate ----------------
    clinvar_plp = (cvc == "PLP")
    clinvar_vus = cvc in ("VUS", "CONFLICT")
    splice_signal = ("SPLICEAI_HIGH" in geno) or ("SPLICEAI_MOD" in geno)

    keep = False
    if clinvar_plp:
        keep = True
    elif domain_bypass:
        # e.g. an established GWAS risk allele: keep even though it is common.
        keep = True
    elif over_ceiling:
        keep = False
    elif is_coding_altering:
        keep = True
    elif splice_signal:
        keep = True
    elif clinvar_vus:
        keep = True
    elif domain_pheno:
        keep = True
    elif noncoding_hits and is_rare_t2 and len(noncoding_hits) >= 2:
        keep = True

    if not keep:
        return False, "Filtered", [], {}

    # ---------------- Tiering ----------------
    lof_in_hi = is_lof and (hugo in haploinsufficient)
    strong_missense = consensus and (is_rare_t2 or is_rare_t1)
    if clinvar_plp:
        tier = "Tier3" if over_ceiling else "Tier1"
    elif lof_in_hi:
        tier = "Tier1"
        geno.append("PVS1_HAPLOINSUFFICIENT")
    elif "SPLICEAI_HIGH" in geno:
        tier = "Tier1"
    elif strong_missense:
        tier = "Tier1"
    elif clinvar_vus:
        tier = "Tier2"
    elif (pred_hits and (is_rare_t2 or is_rare_t1)) or "SPLICEAI_MOD" in geno:
        tier = "Tier2"
    elif is_lof or (is_missense and (is_rare_t2 or is_rare_t1)):
        tier = "Tier2"
    else:
        tier = "Tier3"

    if over_ceiling:
        reasons.append("RISK_ALLELE_COMMON" if (domain_bypass and not clinvar_plp)
                       else "COMMON_AF_FLAG")

    reasons = pheno + geno + reasons
    zyg = zygosity_label(row.get("zygosity"))
    evidence = {
        "gnomad4_af": gnomad,
        "allofus_af": aou,
        "spliceai_max": spliceai_max,
        "predictors": pred_hits,
        "clinvar_class": cvc,
        "hpo_context": hpo_hits,
        "go_context": go_hits,
        "panel_support": panel.get(hugo, {}).get("support"),
        "zygosity": zyg,
        "vaf": compute_vaf(row.get("vaf"), row.get("alt_reads"), row.get("tot_reads")),
    }
    return True, tier, reasons, evidence


def _build_runtime(cfg, vset, gset):
    """Resolve config-declared domain columns to those actually present."""
    def present(col):
        return col in vset or col in gset

    domain_predictors = []
    for dp in cfg.get("domain_predictors", []) or []:
        score = [c for c in dp.get("score_cols", []) if present(c)]
        text = [c for c in dp.get("text_cols", []) if present(c)]
        if score or text:
            domain_predictors.append({
                "code": dp["code"], "min": float(dp.get("min", 0.9)),
                "score_cols": score, "text_cols": text,
            })
    domain_evidence = []
    for ev in cfg.get("domain_evidence", []) or []:
        cols = [c for c in ev.get("columns", []) if present(c)]
        if cols:
            domain_evidence.append({
                "code": ev["code"], "kind": ev.get("kind", "phenotype"),
                "columns": cols,
                "bypass_frequency": bool(ev.get("bypass_frequency", False)),
            })
    return domain_predictors, domain_evidence


def build_select(schema, domain_pull):
    cols = []
    for key in PULL_KEYS:
        col = schema.get(key)
        if not col:
            cols.append(f"NULL AS {key}")
        elif key in GENE_TABLE_KEYS:
            cols.append(f'g."{col}" AS {key}')
        else:
            cols.append(f'v."{col}" AS {key}')
    for alias, col, table in domain_pull:
        cols.append(f'{table}."{col}" AS {alias}')
    return ", ".join(cols)


def run(raw_db, panel_path, schema_path, config_path, out_sqlite, out_json, patient):
    cfg = yaml.safe_load(open(config_path))
    panel = json.load(open(panel_path))["genes"]
    schema = json.load(open(schema_path))
    haploinsufficient = set(cfg.get("haploinsufficient_genes", []) or [])

    hugo_col = schema["hugo"]
    gene_hugo = "base__hugo"

    conn = sqlite3.connect(f"file:{raw_db}?mode=ro", uri=True)
    conn.row_factory = sqlite3.Row
    cur = conn.cursor()
    vset = {r[1] for r in cur.execute("PRAGMA table_info(variant)")}
    gset = {r[1] for r in cur.execute("PRAGMA table_info(gene)")}

    # Multi-sample databases keep genotype in the `sample` table rather than on
    # the variant row. When the variant table has no vcfinfo zygosity but a
    # sample table does, build a uid -> genotype map (first sample per uid).
    sample_geno = {}
    if not schema.get("zygosity") and schema.get("sample_uid") and schema.get("sample_zygosity"):
        s_map = {
            "uid": schema["sample_uid"],
            "zygosity": schema.get("sample_zygosity"),
            "alt_reads": schema.get("sample_alt_reads"),
            "tot_reads": schema.get("sample_tot_reads"),
            "vaf": schema.get("sample_vaf"),
        }
        sel = ", ".join(f'"{col}" AS {alias}' for alias, col in s_map.items() if col)
        try:
            for sr in cur.execute(f"SELECT {sel} FROM sample"):
                d = {k: sr[k] for k in sr.keys()}
                uid = d.pop("uid", None)
                if uid is not None and uid not in sample_geno:
                    sample_geno[uid] = d
        except sqlite3.OperationalError:
            sample_geno = {}

    domain_predictors, domain_evidence = _build_runtime(cfg, vset, gset)

    # Assign stable aliases to each present domain column and remember them.
    domain_pull = []
    alias_i = 0

    def _alias(col):
        nonlocal alias_i
        a = f"dom{alias_i}"
        alias_i += 1
        table = "g" if (col in gset and col not in vset) else "v"
        domain_pull.append((a, col, table))
        return a

    for dp in domain_predictors:
        dp["score_aliases"] = [_alias(c) for c in dp["score_cols"]]
        dp["text_aliases"] = [_alias(c) for c in dp["text_cols"]]
    for ev in domain_evidence:
        ev["aliases"] = [_alias(c) for c in ev["columns"]]
    runtime = {"domain_predictors": domain_predictors, "domain_evidence": domain_evidence}

    select_cols = build_select(schema, domain_pull)
    gene_join = ""
    need_gene = any(schema.get(k) for k in GENE_TABLE_KEYS) or \
        any(t == "g" for _, _, t in domain_pull)
    if need_gene and "gene" in schema.get("tables", []):
        gene_join = f'LEFT JOIN gene g ON g."{gene_hugo}" = v."{hugo_col}"'

    genes = list(panel.keys())
    placeholders = ",".join("?" * len(genes))
    sql = (f"SELECT {select_cols} FROM variant v {gene_join} "
           f'WHERE v."{hugo_col}" IN ({placeholders})')
    cur.execute(sql, genes)

    kept = []
    scanned = 0
    for r in cur:
        scanned += 1
        row = {k: r[k] for k in r.keys()}
        if sample_geno:
            sg = sample_geno.get(row.get("uid"))
            if sg:
                for k, v in sg.items():
                    if row.get(k) in (None, ""):
                        row[k] = v
        keep, tier, reasons, evidence = evaluate_variant(
            row, cfg, panel, haploinsufficient, runtime)
        if keep:
            row["tier"] = tier
            row["reason_codes"] = ";".join(reasons)
            row["evidence"] = evidence
            kept.append(row)
    conn.close()

    tier_order = {"Tier1": 0, "Tier2": 1, "Tier3": 2}
    kept.sort(key=lambda x: (tier_order.get(x["tier"], 9),
                             -(x["evidence"].get("panel_support") or 0),
                             x.get("hugo") or ""))

    _write_sqlite(out_sqlite, kept, panel, patient)
    _write_json(out_json, kept, panel, patient, scanned, cfg)

    counts = {t: sum(1 for k in kept if k["tier"] == t) for t in ("Tier1", "Tier2", "Tier3")}
    active = [dp["code"] for dp in domain_predictors] + [ev["code"] for ev in domain_evidence]
    print(f"[filter] domain={cfg.get('domain')} active_domain_signals={active or 'none'}")
    print(f"[filter] scanned {scanned} panel-gene variant rows")
    print(f"[filter] actionable kept={len(kept)}  {counts}")
    print(f"[filter] wrote {out_sqlite} and {out_json}")
    return kept


def _write_sqlite(path, kept, panel, patient):
    if os.path.exists(path):
        os.remove(path)
    conn = sqlite3.connect(path)
    cur = conn.cursor()
    scalar_keys = list(PULL_KEYS) + ["tier", "reason_codes"]
    coldefs = ", ".join(f'"{k}" TEXT' for k in scalar_keys)
    cur.execute(f"CREATE TABLE variant ({coldefs})")
    cur.executemany(
        f"INSERT INTO variant ({','.join(chr(34)+k+chr(34) for k in scalar_keys)}) "
        f"VALUES ({','.join('?'*len(scalar_keys))})",
        [[(str(k[c]) if k.get(c) is not None else None) for c in scalar_keys] for k in kept],
    )
    cur.execute("CREATE TABLE panel_gene (hugo TEXT, support INT, hpo TEXT, go TEXT)")
    cur.executemany(
        "INSERT INTO panel_gene VALUES (?,?,?,?)",
        [(g, v["support"], ";".join(v["hpo"]), ";".join(v["go"])) for g, v in panel.items()],
    )
    cur.execute("CREATE TABLE report_meta (patient TEXT, kept INT)")
    cur.execute("INSERT INTO report_meta VALUES (?,?)", (patient, len(kept)))
    conn.commit()
    conn.close()


def _write_json(path, kept, panel, patient, scanned, cfg):
    records = []
    for k in kept:
        rec = {kk: k.get(kk) for kk in PULL_KEYS}
        rec["tier"] = k["tier"]
        rec["reason_codes"] = k["reason_codes"].split(";") if k["reason_codes"] else []
        rec["evidence"] = k["evidence"]
        records.append(rec)
    out = {
        "patient": patient,
        "domain": cfg.get("domain"),
        "report_title": cfg.get("report_title", "Ontology-Driven Actionable Report"),
        "panel_gene_count": len(panel),
        "scanned_panel_variants": scanned,
        "actionable_count": len(records),
        "tier_counts": {t: sum(1 for r in records if r["tier"] == t)
                        for t in ("Tier1", "Tier2", "Tier3")},
        "records": records,
    }
    with open(path, "w") as f:
        json.dump(out, f, indent=2)


def main():
    ap = argparse.ArgumentParser(description="Ontology-driven actionable variant filter")
    ap.add_argument("--raw-db", required=True)
    ap.add_argument("--panel", required=True)
    ap.add_argument("--schema", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--out-sqlite", required=True)
    ap.add_argument("--out-json", required=True)
    ap.add_argument("--patient", default="Patient")
    args = ap.parse_args()
    run(args.raw_db, args.panel, args.schema, args.config,
        args.out_sqlite, args.out_json, args.patient)


if __name__ == "__main__":
    main()
