#!/usr/bin/env python3
"""
Offline validation suite for ontology_report (domain-agnostic engine).

Builds a tiny mock OpenCRAVAT-style SQLite (no WGS, no annotation needed) with
hand-crafted variants covering each decision path, then asserts that the
ontology filter tiers and reason codes come out as expected. Runs in <1s.

Usage:  python3 tests/test_flow.py
"""
import json
import os
import sqlite3
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "lib"))

import ontology_filter as of  # noqa: E402
import schema_probe as sp  # noqa: E402

CONFIG = os.path.join(ROOT, "config", "cardiology.yaml")

# Mock variant columns (subset of a real OC variant table).
VCOLS = [
    "base__uid", "base__hugo", "base__so", "base__coding", "base__achange",
    "base__cchange", "base__transcript", "base__chrom", "base__pos",
    "base__ref_base", "base__alt_base",
    "gnomad4__af", "allofus250k__gvs_all_af",
    "clinvar__sig", "clinvar__id", "clinvar__disease_names", "clinvar__rev_stat",
    "alphamissense__am_pathogenicity", "alphamissense__am_class", "revel__score",
    "cardioboost__cardiomyopathy", "cardioboost__cardiomyopathy1",
    "bayesdel__bayesdel_addAF_score", "metarnn__score", "esm1b__score",
    "varity_r__varity_r",
    "spliceai__ds_ag", "spliceai__ds_al", "spliceai__ds_dg", "spliceai__ds_dl",
    "cadd__phred", "linsight__value", "ncer__score", "regulomedb__ra",
    "arrvars__lqt", "arrvars__function",
]
GCOLS = ["base__hugo", "hpo__id", "hpo__term", "go__bpo_name", "go__mfo_name",
         "go__cco_name", "clingen__classification", "clingen__disease"]


def _mk_variant(uid, hugo, so, **kw):
    row = {c: None for c in VCOLS}
    row["base__uid"] = uid
    row["base__hugo"] = hugo
    row["base__so"] = so
    row["base__chrom"] = "chr1"
    row["base__pos"] = 1000 + uid
    row["base__ref_base"] = "A"
    row["base__alt_base"] = "G"
    for k, v in kw.items():
        col = {
            "af": "gnomad4__af", "aou": "allofus250k__gvs_all_af",
            "clinvar": "clinvar__sig", "clinvar_id": "clinvar__id",
            "am": "alphamissense__am_pathogenicity", "revel": "revel__score",
            "metarnn": "metarnn__score", "esm1b": "esm1b__score",
            "varity": "varity_r__varity_r", "bayesdel": "bayesdel__bayesdel_addAF_score",
            "ds_ag": "spliceai__ds_ag", "cadd": "cadd__phred",
            "linsight": "linsight__value", "ncer": "ncer__score", "achange": "base__achange",
            "cardioboost": "cardioboost__cardiomyopathy", "arrvars": "arrvars__lqt",
        }.get(k, k)
        row[col] = v
    return [row[c] for c in VCOLS]


def build_mock_db(path):
    conn = sqlite3.connect(path)
    cur = conn.cursor()
    vdef = ", ".join('"%s"' % c for c in VCOLS)
    gdef = ", ".join('"%s"' % c for c in GCOLS)
    cur.execute("CREATE TABLE variant (%s)" % vdef)
    cur.execute("CREATE TABLE gene (%s)" % gdef)

    variants = [
        # 1: ClinVar Pathogenic missense in panel gene, rare -> Tier1
        _mk_variant(1, "MYH7", "MIS", clinvar="Pathogenic", clinvar_id="111",
                    af=0.0, revel=0.9, am=0.9, achange="p.R100C"),
        # 2: LoF (stop-gained) in haploinsufficient LMNA, rare -> Tier1 (PVS1)
        _mk_variant(2, "LMNA", "STG", af=None, achange="p.Q50*"),
        # 3: high SpliceAI splice variant, rare -> Tier1
        _mk_variant(3, "KCNQ1", "SPL", af=1e-5, ds_ag=0.8),
        # 4: multi-predictor consensus missense, rare -> Tier1
        _mk_variant(4, "SCN5A", "MIS", af=1e-5, revel=0.9, am=0.9, metarnn=0.9,
                    esm1b=-9.0, varity=0.9),
        # 5: ClinVar VUS -> Tier2
        _mk_variant(5, "TTN", "MIS", clinvar="Uncertain significance",
                    clinvar_id="222", af=1e-4),
        # 6: single predictor + rare missense -> Tier2
        _mk_variant(6, "RYR2", "MIS", af=1e-5, revel=0.8),
        # 7: common benign missense (AF 20%) not ClinVar P/LP -> dropped
        _mk_variant(7, "MYH7", "MIS", af=0.20, revel=0.9),
        # 8: non-panel gene coding -> dropped (not in panel)
        _mk_variant(8, "BRCA1", "STG", af=None),
        # 9: intronic no signal in panel gene -> dropped (no evidence)
        _mk_variant(9, "MYH7", "INT", af=1e-5),
        # 10: ClinVar Pathogenic but COMMON -> kept, demoted Tier3 + flag
        _mk_variant(10, "MYBPC3", "MIS", clinvar="Pathogenic", clinvar_id="333",
                    af=0.30),
        # 11: non-coding regulatory double-signal, rare -> Tier3
        _mk_variant(11, "DSP", "UT3", af=1e-5, cadd=25, ncer=0.95),
        # 12: intronic, no core signal, but a config-driven domain-evidence
        #     column (ArrVars) is present -> actionable via domain_evidence
        _mk_variant(12, "SCN5A", "INT", af=1e-5, arrvars="LQT3 pathogenic"),
    ]
    cur.executemany(
        f"INSERT INTO variant VALUES ({','.join('?'*len(VCOLS))})", variants
    )
    genes = [
        ("MYH7", "HP:0001639", "Hypertrophic cardiomyopathy", "cardiac muscle contraction", "", "", None, None),
        ("LMNA", "HP:0001644", "Dilated cardiomyopathy", "", "", "", "Definitive", "DCM"),
        ("KCNQ1", "HP:0001657", "Prolonged QT interval", "", "potassium channel activity", "", None, None),
        ("SCN5A", "HP:0011675", "Arrhythmia", "cardiac conduction", "sodium channel activity", "", None, None),
        ("TTN", "HP:0001644", "Cardiomyopathy", "cardiac muscle contraction", "", "", None, None),
        ("RYR2", "HP:0004308", "Ventricular arrhythmia", "", "calcium channel activity", "", None, None),
        ("MYBPC3", "HP:0001639", "Hypertrophic cardiomyopathy", "", "", "", None, None),
        ("DSP", "HP:0011663", "Right ventricular cardiomyopathy", "", "", "", "Definitive", "ARVC"),
    ]
    cur.executemany(f"INSERT INTO gene VALUES ({','.join('?'*len(GCOLS))})", genes)
    conn.commit()
    conn.close()


def expect(cond, msg, failures):
    status = "PASS" if cond else "FAIL"
    print(f"  [{status}] {msg}")
    if not cond:
        failures.append(msg)


def main():
    failures = []
    tmp = tempfile.mkdtemp(prefix="ontotest_")
    db = os.path.join(tmp, "mock.sqlite")
    panel_path = os.path.join(tmp, "panel.json")
    schema_path = os.path.join(tmp, "schema.json")
    out_sqlite = os.path.join(tmp, "act.sqlite")
    out_json = os.path.join(tmp, "act.json")

    build_mock_db(db)

    # Panel: a fixed set that mimics the ontology output (so the test is
    # independent of the installed ontology DB).
    panel = {"genes": {g: {"support": 2, "forced": False, "hpo": ["x"], "go": ["y"]}
                       for g in ["MYH7", "LMNA", "KCNQ1", "SCN5A", "TTN", "RYR2",
                                 "MYBPC3", "DSP"]}}
    json.dump(panel, open(panel_path, "w"))

    schema = sp.probe(db)
    json.dump(schema, open(schema_path, "w"))

    kept = of.run(db, panel_path, schema_path, CONFIG, out_sqlite, out_json, "MOCK")
    by_uid = {int(k["uid"]): k for k in kept}

    print("Assertions:")
    expect(1 in by_uid and by_uid[1]["tier"] == "Tier1",
           "ClinVar Pathogenic rare -> Tier1", failures)
    expect(2 in by_uid and by_uid[2]["tier"] == "Tier1"
           and "PVS1_HAPLOINSUFFICIENT" in by_uid[2]["reason_codes"],
           "LoF in haploinsufficient gene -> Tier1 + PVS1", failures)
    expect(3 in by_uid and by_uid[3]["tier"] == "Tier1"
           and "SPLICEAI_HIGH" in by_uid[3]["reason_codes"],
           "High SpliceAI -> Tier1", failures)
    expect(4 in by_uid and by_uid[4]["tier"] == "Tier1"
           and "PP3_CONSENSUS" in by_uid[4]["reason_codes"],
           "Multi-predictor consensus rare -> Tier1", failures)
    expect(5 in by_uid and by_uid[5]["tier"] == "Tier2",
           "ClinVar VUS -> Tier2", failures)
    expect(6 in by_uid and by_uid[6]["tier"] == "Tier2",
           "Single predictor + rare missense -> Tier2", failures)
    expect(7 not in by_uid, "Common benign missense -> dropped", failures)
    expect(8 not in by_uid, "Non-panel gene -> dropped", failures)
    expect(9 not in by_uid, "Intronic no-signal panel variant -> dropped", failures)
    expect(10 in by_uid and by_uid[10]["tier"] == "Tier3"
           and "COMMON_AF_FLAG" in by_uid[10]["reason_codes"],
           "ClinVar Pathogenic but common -> kept, Tier3 + flag", failures)
    expect(11 in by_uid and by_uid[11]["tier"] == "Tier3",
           "Non-coding double regulatory rare -> Tier3", failures)
    expect(12 in by_uid and "ARRVARS_KNOWN" in by_uid[12]["reason_codes"],
           "Config-driven domain evidence (ArrVars) -> actionable", failures)

    n_checks = 12
    print()
    if failures:
        print(f"RESULT: {len(failures)} FAILURE(S)")
        for f in failures:
            print(f"   - {f}")
        sys.exit(1)
    print(f"RESULT: ALL {n_checks} CHECKS PASSED")


if __name__ == "__main__":
    main()
