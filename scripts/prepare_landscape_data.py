#!/usr/bin/env python3
"""
prepare_landscape_data.py
Precomputes high-density clinical genomics chart aggregations and variant datasets
for the Landscape Clinical Dashboard (dual-cohort Daniel Ehrle & Melinda Ehrle).
"""

import json
import os
import sys
from collections import defaultdict, Counter

SCRIPT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(SCRIPT_DIR, "data")
REPORTS_DIR = os.path.join(SCRIPT_DIR, "reports")

DE_JSON = os.path.join(REPORTS_DIR, "Daniel_Ehrle-09-10-2026", "Daniel_Ehrle_ontology_pharma_alphagenome.json")
ME_JSON = os.path.join(REPORTS_DIR, "Melinda_Ehrle-09-10-2026", "Melinda_Ehrle_ontology_pharma_alphagenome.json")

# Therapeutic drug class keyword matcher
DRUG_CLASSES = [
    ("Cardiovascular & Lipids", ["statin", "simvastatin", "atorvastatin", "warfarin", "clopidogrel", "prasugrel", "ticagrelor", "metoprolol", "carvedilol", "digoxin", "lipid", "coagulat"]),
    ("Neuropsychiatry", ["citalopram", "escitalopram", "sertraline", "fluoxetine", "paroxetine", "venlafaxine", "amitriptyline", "nortriptyline", "aripiprazole", "haloperidol", "antidepress", "psych"]),
    ("Pain & Analgesia", ["codeine", "tramadol", "hydrocodone", "oxycodone", "morphine", "fentanyl", "ibuprofen", "celecoxib", "meloxicam", "opioid", "analges"]),
    ("Oncology & Immune", ["tamoxifen", "fluorouracil", "capecitabine", "irinotecan", "mercaptopurine", "thioguanine", "azathioprine", "tacrolimus", "methotrexate", "cancer"]),
    ("Gastroenterology", ["omeprazole", "pantoprazole", "lansoprazole", "esomeprazole", "rabeprazole", "antiemetic", "ondansetron"]),
    ("Infectious Disease & Other", ["efavirenz", "abacavir", "voriconazole", "peginterferon", "antibiotic", "antiviral", "contraceptive", "estrogen"])
]

def map_drug_to_class(drug_name):
    d_lower = (drug_name or "").lower()
    for class_name, keywords in DRUG_CLASSES:
        if any(kw in d_lower for kw in keywords):
            return class_name
    return "Infectious Disease & Other"

def process_patient(json_path, patient_id, patient_name):
    print(f"Processing {patient_name} ({patient_id}) from {json_path}...")
    with open(json_path, "r", encoding="utf-8") as f:
        data = json.load(f)

    genes_list = data.get("genes", [])
    if isinstance(genes_list, dict):
        genes_list = list(genes_list.values())

    # Build flat variant list
    flat_records = []
    for g in genes_list:
        for v in g.get("variants", []):
            flat_records.append((g, v))

    total_variants = len(flat_records)
    total_genes = len(genes_list)

    # -------------------------------------------------------------------------
    # ROW 1: ONTOLOGY (HPO Clinical Domains)
    # -------------------------------------------------------------------------
    # Define standardized clinical presentation categories matching reference image
    hpo_domain_defs = [
        ("Neurology", ["nervous", "seizure", "epilepsy", "encephal", "brain", "neuropath", "ataxia", "dementia", "cns", "intellectual", "cntn1", "c19orf12"]),
        ("Metabolic", ["metabol", "inborn", "acid", "peroxis", "mitochon", "glycogen", "storage", "cblif", "hfe", "mthfr", "diabetes"]),
        ("Cardiovascular", ["cardio", "vascular", "heart", "arrhythm", "channelopath", "cardiomyopath", "scn5a", "apob", "ldlr", "long qt", "brugada"]),
        ("GI / Liver", ["gastro", "liver", "hepat", "digest", "intest", "bowel", "colon", "pancrea"]),
        ("Heme / Immune", ["immune", "autoimmun", "hematolog", "anemia", "coagulat", "leuk", "lymph", "platelet", "thromb", "hla", "ptpn22", "f5"]),
        ("Morphogenesis", ["skelet", "bone", "cranio", "dysmorphism", "syndrom", "connect", "collagen", "ehlers", "marfan"]),
        ("Neoplasm", ["neoplasm", "cancer", "tumor", "carcinoma", "lynch", "brca", "pms2", "rad51"]),
        ("Renal / Urinary", ["urin", "renal", "kidney", "nephr", "glomerul"]),
        ("Other / Multisystem", [])
    ]

    def get_variant_hpo_domain(gene_sym, gene_name, v):
        text = f"{gene_sym} {gene_name} {' '.join(v.get('consequence', []))} {v.get('cchange', '')} {v.get('achange', '')} {v.get('clinvar', '')}".lower()
        for domain, kws in hpo_domain_defs[:-1]:
            if any(kw in text for kw in kws):
                return domain
        return "Other / Multisystem"

    hpo_domain_counts = defaultdict(lambda: {"total": 0, "actionable": 0})
    for g, v in flat_records:
        domain = get_variant_hpo_domain(g.get("symbol", ""), g.get("name", ""), v)
        is_actionable = (
            v.get("tier") in ["Tier1", "Tier2"] or
            v.get("category") in ["concern", "protective"] or
            "pathogenic" in (v.get("clinvar") or "").lower() or
            (v.get("aviPhred") is not None and float(v.get("aviPhred")) >= 20.0) or
            (v.get("cadd") is not None and float(v.get("cadd")) >= 25.0)
        )
        hpo_domain_counts[domain]["total"] += 1
        if is_actionable:
            hpo_domain_counts[domain]["actionable"] += 1

    # Format Row 1 Stacked Bar
    ordered_domains = ["Neurology", "Metabolic", "Cardiovascular", "GI / Liver", "Heme / Immune", "Morphogenesis", "Neoplasm", "Renal / Urinary", "Other / Multisystem"]
    row1_bar = []
    for d in ordered_domains:
        tot = hpo_domain_counts[d]["total"]
        act = hpo_domain_counts[d]["actionable"]
        sec = tot - act
        rate = round((act / tot * 100), 1) if tot > 0 else 0.0
        row1_bar.append({
            "category": d,
            "actionable": act,
            "secondary": sec,
            "total": tot,
            "fraction": f"{act}/{tot}",
            "ratePct": rate
        })

    # Row 1 Pie B (Presentation Prior to RGS)
    row1_pie_prior = []
    colors_palette = ["#17222B", "#B3261E", "#475569", "#64748B", "#94A3B8", "#B5722A", "#1C5B5E", "#0284C7", "#CBD5E1"]
    tot_all_hpo = sum(hpo_domain_counts[d]["total"] for d in ordered_domains) or 1
    for i, d in enumerate(ordered_domains):
        cnt = hpo_domain_counts[d]["total"]
        pct = round(cnt / tot_all_hpo * 100, 1)
        row1_pie_prior.append({
            "label": d,
            "count": cnt,
            "total": tot_all_hpo,
            "pct": pct,
            "callout": f"{d} {pct}%, {cnt}/{tot_all_hpo}",
            "color": colors_palette[i % len(colors_palette)]
        })

    # Row 1 Pie C (Diagnosis Category After RGS - Confirmed Actionable Findings)
    tot_act_hpo = sum(hpo_domain_counts[d]["actionable"] for d in ordered_domains) or 1
    row1_pie_post = []
    diag_colors = ["#B3261E", "#17222B", "#1C5B5E", "#B5722A", "#0284C7", "#475569", "#047857", "#6D28D9", "#94A3B8"]
    for i, d in enumerate(ordered_domains):
        cnt = hpo_domain_counts[d]["actionable"]
        if cnt > 0:
            pct = round(cnt / tot_act_hpo * 100, 1)
            row1_pie_post.append({
                "label": d,
                "count": cnt,
                "total": tot_act_hpo,
                "pct": pct,
                "callout": f"{d} {pct}%, {cnt}/{tot_act_hpo}",
                "color": diag_colors[i % len(diag_colors)]
            })

    # -------------------------------------------------------------------------
    # ROW 2: ORGAN SYSTEMS (9 Anatomical Systems & Risk Matrix)
    # -------------------------------------------------------------------------
    organ_risk_matrix = data.get("organRiskMatrix", [])
    organ_names_standard = [
        "Heart & Cardiovascular",
        "Brain & Nervous System",
        "Digestive & Metabolism",
        "Immune & Autoimmunity",
        "Skeleton & Connective Tissue",
        "Lungs & Respiratory",
        "Kidneys & Urinary",
        "Cancer Predisposition",
        "Other / Multisystem"
    ]

    organ_counts = defaultdict(lambda: {"total": 0, "actionable": 0})
    for g, v in flat_records:
        org = g.get("organSystem") or "Other / Multisystem"
        matched_org = "Other / Multisystem"
        for std_o in organ_names_standard:
            if std_o.lower() in org.lower() or any(w in org.lower() for w in std_o.lower().split("&")):
                matched_org = std_o
                break
        
        is_actionable = (
            v.get("tier") in ["Tier1", "Tier2"] or
            v.get("category") in ["concern", "protective"] or
            "pathogenic" in (v.get("clinvar") or "").lower() or
            (v.get("aviPhred") is not None and float(v.get("aviPhred")) >= 20.0)
        )
        organ_counts[matched_org]["total"] += 1
        if is_actionable:
            organ_counts[matched_org]["actionable"] += 1

    row2_bar = []
    for org in organ_names_standard:
        tot = organ_counts[org]["total"]
        act = organ_counts[org]["actionable"]
        sec = tot - act
        rate = round((act / tot * 100), 1) if tot > 0 else 0.0
        row2_bar.append({
            "category": org.replace("& ", "").replace(" System", "").replace(" Tissue", ""),
            "fullName": org,
            "actionable": act,
            "secondary": sec,
            "total": tot,
            "fraction": f"{act}/{tot}",
            "ratePct": rate
        })

    tot_all_org = sum(organ_counts[o]["total"] for o in organ_names_standard) or 1
    row2_pie_prior = []
    for i, org in enumerate(organ_names_standard):
        cnt = organ_counts[org]["total"]
        pct = round(cnt / tot_all_org * 100, 1)
        short_lbl = org.replace("& ", "").replace(" System", "").replace(" Tissue", "")
        row2_pie_prior.append({
            "label": short_lbl,
            "fullName": org,
            "count": cnt,
            "total": tot_all_org,
            "pct": pct,
            "callout": f"{short_lbl} {pct}%, {cnt}/{tot_all_org}",
            "color": colors_palette[i % len(colors_palette)]
        })

    tot_act_org = sum(organ_counts[o]["actionable"] for o in organ_names_standard) or 1
    row2_pie_post = []
    for i, org in enumerate(organ_names_standard):
        cnt = organ_counts[org]["actionable"]
        if cnt > 0:
            pct = round(cnt / tot_act_org * 100, 1)
            short_lbl = org.replace("& ", "").replace(" System", "").replace(" Tissue", "")
            row2_pie_post.append({
                "label": short_lbl,
                "fullName": org,
                "count": cnt,
                "total": tot_act_org,
                "pct": pct,
                "callout": f"{short_lbl} {pct}%, {cnt}/{tot_act_org}",
                "color": diag_colors[i % len(diag_colors)]
            })

    # -------------------------------------------------------------------------
    # ROW 3: PHARMACOGENOMICS (CPIC / DPWG Clinical Actionability)
    # -------------------------------------------------------------------------
    pgx_data = data.get("pharmacogenomics", {})
    ext_pharma = pgx_data.get("externalPharmaReport", {})
    drug_cards = ext_pharma.get("drugCards", [])

    pharma_class_counts = defaultdict(lambda: {"total": 0, "actionable": 0})
    phenotype_counts = Counter()
    action_counts = Counter()

    for card in drug_cards:
        drug_name = card.get("drug", "")
        action = card.get("action", "STANDARD")
        pheno = card.get("phenotype", "Normal Metabolizer")
        
        cls = map_drug_to_class(drug_name)
        pharma_class_counts[cls]["total"] += 1
        
        is_act = action in ["AVOID", "DOSE_DOWN", "DOSE_UP", "CAUTION"]
        if is_act:
            pharma_class_counts[cls]["actionable"] += 1

        phenotype_counts[pheno] += 1
        action_counts[action] += 1

    standard_pharma_classes = [c[0] for c in DRUG_CLASSES]
    row3_bar = []
    for cls in standard_pharma_classes:
        tot = pharma_class_counts[cls]["total"]
        act = pharma_class_counts[cls]["actionable"]
        sec = tot - act
        rate = round((act / tot * 100), 1) if tot > 0 else 0.0
        row3_bar.append({
            "category": cls.replace("& ", "").replace(" Disease", ""),
            "fullName": cls,
            "actionable": act,
            "secondary": sec,
            "total": tot,
            "fraction": f"{act}/{tot}",
            "ratePct": rate
        })

    tot_phenos = sum(phenotype_counts.values()) or 1
    row3_pie_pheno = []
    pheno_colors = {
        "Normal Metabolizer": "#17222B",
        "Intermediate Metabolizer": "#B5722A",
        "Poor Metabolizer": "#B3261E",
        "Rapid Metabolizer": "#0284C7",
        "Ultrarapid Metabolizer": "#7C3AED",
        "Unknown": "#64748B",
        "Indeterminate": "#94A3B8"
    }
    for i, (pheno, cnt) in enumerate(phenotype_counts.most_common()):
        pct = round(cnt / tot_phenos * 100, 1)
        clr = pheno_colors.get(pheno, colors_palette[i % len(colors_palette)])
        row3_pie_pheno.append({
            "label": pheno,
            "count": cnt,
            "total": tot_phenos,
            "pct": pct,
            "callout": f"{pheno} {pct}%, {cnt}/{tot_phenos}",
            "color": clr
        })

    tot_actions = sum(action_counts.values()) or 1
    row3_pie_action = []
    action_colors = {
        "AVOID": "#B3261E",
        "DOSE_DOWN": "#DC2626",
        "DOSE_UP": "#EA580C",
        "CAUTION": "#D97706",
        "MONITOR": "#2563EB",
        "FAVOUR": "#059669",
        "STANDARD": "#17222B"
    }
    for i, (act, cnt) in enumerate(action_counts.most_common()):
        pct = round(cnt / tot_actions * 100, 1)
        clr = action_colors.get(act, colors_palette[i % len(colors_palette)])
        row3_pie_action.append({
            "label": act,
            "count": cnt,
            "total": tot_actions,
            "pct": pct,
            "callout": f"{act} {pct}%, {cnt}/{tot_actions}",
            "color": clr
        })

    # -------------------------------------------------------------------------
    # REPORT VIEW METRICS PRECOMPUTATION
    # -------------------------------------------------------------------------
    report_stats = {
        "visual_explorer": {
            "title": "Visual Explorer — Universal Genomic & Variant Catalog",
            "desc": "Comprehensive interactive catalog of all prioritized and annotated genomic variants across HPO, GO, and ClinVar classifications, with multidimensional impact scoring and parental phase attribution.",
            "totalGenes": len(genes_list),
            "totalVariants": total_variants,
            "topMetricLabel": "Priority Tier"
        },
        "rare": {
            "title": "Rare & Ultra-Rare Alleles (gnomAD AF ≤ 1%)",
            "desc": "Low-frequency and singleton alleles filtered by gnomAD population allele frequency (AF ≤ 1.0%), highlighting rare heterozygous or homozygous variants with potential high-penetrance phenotypic expression.",
            "totalGenes": len([g for g in genes_list if any(v.get("maf") is not None and v.get("maf") <= 0.01 for v in g.get("variants", []))]),
            "totalVariants": len([v for g, v in flat_records if v.get("maf") is not None and v.get("maf") <= 0.01]),
            "topMetricLabel": "gnomAD AF"
        },
        "clinvar": {
            "title": "ClinVar Clinical Pathogenicity & Actionability Report",
            "desc": "Curated variants in NCBI ClinVar with established clinical significance, including Pathogenic, Likely Pathogenic, Established Risk Factors, and Curated Protective alleles.",
            "totalGenes": len([g for g in genes_list if any(v.get("clinvar") and v.get("clinvar") not in ["None", "Not reviewed", ""] for v in g.get("variants", []))]),
            "totalVariants": len([v for g, v in flat_records if v.get("clinvar") and v.get("clinvar") not in ["None", "Not reviewed", ""]]),
            "topMetricLabel": "ClinVar Classification"
        },
        "cadd": {
            "title": "CADD High Functional Deleteriousness (Phred ≥ 15)",
            "desc": "Variants prioritized by Combined Annotation Dependent Depletion (CADD Phred ≥ 15, top 3% to 0.1% most deleterious mutations in the human genome), assessing both coding and regulatory disruptions.",
            "totalGenes": len([g for g in genes_list if any(v.get("cadd") is not None and v.get("cadd") >= 15.0 for v in g.get("variants", []))]),
            "totalVariants": len([v for g, v in flat_records if v.get("cadd") is not None and v.get("cadd") >= 15.0]),
            "topMetricLabel": "CADD Phred"
        },
        "spliceai": {
            "title": "SpliceAI Deep Learning Splicing Disruption Report (Δ ≥ 0.10)",
            "desc": "Deep neural network predictions of non-canonical and cryptic splice site disruption across Acceptor Gain (AG), Acceptor Loss (AL), Donor Gain (DG), and Donor Loss (DL) delta scores.",
            "totalGenes": len([g for g in genes_list if any(v.get("spliceai") is not None and v.get("spliceai") >= 0.10 for v in g.get("variants", []))]),
            "totalVariants": len([v for g, v in flat_records if v.get("spliceai") is not None and v.get("spliceai") >= 0.10]),
            "topMetricLabel": "SpliceAI Max Δ"
        },
        "atlas": {
            "title": "Google DeepMind AlphaGenome Atlas Variant Impact (AVI Score)",
            "desc": "Whole-genome regulatory and functional impact predictions scored by Google DeepMind AlphaGenome (Phred-scaled AVI score and genome-wide percentiles), capturing chromatin accessibility, histone marks, and transcription.",
            "totalGenes": len([g for g in genes_list if any(v.get("aviPhred") is not None for v in g.get("variants", []))]),
            "totalVariants": len([v for g, v in flat_records if v.get("aviPhred") is not None]),
            "topMetricLabel": "AlphaGenome AVI"
        },
        "alphamissense": {
            "title": "AlphaMissense Structural Missense Pathogenicity",
            "desc": "DeepMind AlphaFold-derived structural pathogenicity predictions for non-synonymous amino acid substitutions, categorizing variants as Likely Pathogenic, Ambiguous, or Likely Benign.",
            "totalGenes": len([g for g in genes_list if any(v.get("alphamissense") is not None for v in g.get("variants", []))]),
            "totalVariants": len([v for g, v in flat_records if v.get("alphamissense") is not None]),
            "topMetricLabel": "AlphaMissense Score"
        },
        "revel": {
            "title": "REVEL Ensembled Missense Pathogenicity Report (Score ≥ 0.5)",
            "desc": "Ensemble score integrating 13 individual pathogenicity prediction algorithms for missense variants (scores ≥ 0.50 indicate significant pathogenic probability; ≥ 0.75 indicates high confidence).",
            "totalGenes": len([g for g in genes_list if any(v.get("revel") is not None and v.get("revel") >= 0.50 for v in g.get("variants", []))]),
            "totalVariants": len([v for g, v in flat_records if v.get("revel") is not None and v.get("revel") >= 0.50]),
            "topMetricLabel": "REVEL Score"
        },
        "ai_deep": {
            "title": "Comprehensive AI Deep Research Genomic Dossier (Alt B)",
            "desc": "Multi-page publication-grade clinical genomics evidence synthesis correlating multi-engine AI scores with ClinVar, OMIM, GWAS, and PubMed literature, adhering to VSCP-DF clinical decision standards.",
            "htmlFile": f"reports/{os.path.basename(os.path.dirname(json_path))}/{patient_id}_deep_research_report.html"
        },
        "ai_clinical": {
            "title": "Executive AI Clinical Action Brief (Alt C 3-Page Brief)",
            "desc": "Actionable 3-page clinical summary designed for medical providers and genetic counselors, highlighting definitive diagnoses, pharmacogenomic contraindications, and multi-system risk mitigation.",
            "htmlFile": f"reports/{os.path.basename(os.path.dirname(json_path))}/{patient_id}_clinical_brief.html"
        }
    }

    return {
        "patientId": patient_id,
        "patientName": patient_name,
        "meta": {
            "sampleId": patient_id,
            "sampleName": patient_name,
            "totalGenes": total_genes,
            "totalVariants": total_variants,
            "date": data.get("metadata", {}).get("date") or "2026-10-09",
            "pipeline": "Genomic Ontology Engine v5.2 (Phased WGS + AlphaGenome AVI)",
            "deepHtml": f"reports/{os.path.basename(os.path.dirname(json_path))}/{patient_id}_deep_research_report.html",
            "clinicalHtml": f"reports/{os.path.basename(os.path.dirname(json_path))}/{patient_id}_clinical_brief.html"
        },
        "charts": {
            "row1_ontology": {
                "title": "Row 1 — Groupings by Ontology (HPO Disease Classification)",
                "bar": row1_bar,
                "piePrior": row1_pie_prior,
                "piePost": row1_pie_post
            },
            "row2_organ": {
                "title": "Row 2 — Groupings by Organ System (9 Anatomical Systems & Risk Tiers)",
                "bar": row2_bar,
                "piePrior": row2_pie_prior,
                "piePost": row2_pie_post
            },
            "row3_pharma": {
                "title": "Row 3 — Groupings by Pharmacogenomics (CPIC / DPWG Clinical Actionability)",
                "bar": row3_bar,
                "piePrior": row3_pie_pheno,
                "piePost": row3_pie_action
            }
        },
        "reportStats": report_stats,
        "genes": genes_list
    }

def main():
    print("==================================================================")
    print("PREPARING DATASET FOR LANDSCAPE CLINICAL DASHBOARD")
    print("==================================================================")
    os.makedirs(DATA_DIR, exist_ok=True)

    de_data = process_patient(DE_JSON, "Daniel_Ehrle", "Daniel Ehrle")
    me_data = process_patient(ME_JSON, "Melinda_Ehrle", "Melinda Ehrle")

    cohort_payload = {
        "Daniel_Ehrle": de_data,
        "Melinda_Ehrle": me_data
    }

    out_js = os.path.join(DATA_DIR, "cohort_landscape_data.js")
    print(f"\nWriting cohort JS bundle to: {out_js}...")
    with open(out_js, "w", encoding="utf-8") as f:
        f.write("/**\n * Cohort Landscape Data for Genomic Ontology Explorer\n * Precomputed for Daniel Ehrle & Melinda Ehrle\n */\n")
        f.write("window.COHORT_DATA = ")
        json.dump(cohort_payload, f, indent=None)
        f.write(";\n")

    size_mb = os.path.getsize(out_js) / (1024 * 1024)
    print(f"Cohort JS bundle generated successfully! ({size_mb:.2f} MB)")
    print("Done.\n")

if __name__ == "__main__":
    main()
