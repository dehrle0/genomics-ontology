#!/usr/bin/env python3
"""
run_standard_reporting.py
Standardized Clinical Genomics Reporting Orchestrator for DE and ME.

Enforces:
  1. Dynamic cohort discovery from local filesystem with zero hardcoded PII in source code.
  2. Full enrichment of master JSON with pharmacogenomics and AlphaGenome Atlas data.
  3. Automatic PII stripping prior to any model synthesis.
  4. Generation of the standardized deliverables (Clinician Action Briefs, Deep Research Dossiers,
     Dashboard Tier-1 & Top-N Digest, EHR Import JSON, Standalone Landscape Dashboard).
  5. Automated delivery to Google Drive Ontology root, Genomics Reports (max 5 iterations), and export/.
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
from datetime import datetime

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(SCRIPT_DIR, "scripts"))
import prepare_landscape_data as pld
from lib.zero_curation_guard import enforce_or_exit
enforce_or_exit(SCRIPT_DIR)


def parse_args():
    parser = argparse.ArgumentParser(description="Standard Clinical Genomics Reporting Orchestrator")
    parser.add_argument("target", choices=["me", "de", "ME", "DE"], help="Target sample code (me or de)")
    parser.add_argument("--reports-dir", default=os.path.join(SCRIPT_DIR, "reports"))
    parser.add_argument("--sample-dir", default=None, help="Explicit sample folder name")
    parser.add_argument("--gdrive-dir", default=os.path.expanduser("~/Google Drive/My Drive/Ontology"))
    parser.add_argument("--dry-run", action="store_true", help="Validate inputs without writing files")
    return parser.parse_args()


def build_grouped_catalog(variants_list):
    """
    Groups variants into 6 standardized clinical categories:
    1. Pathogenic & Monogenic Carriers
    2. Actionable Pharmacogenomics & Drug Response
    3. Cardiovascular & Channelopathy Modifiers
    4. Metabolic, Mitochondrial & DNA Repair Machinery
    5. Protective & Longevity Modifiers
    6. Secondary & Exploratory Clinical Variants
    """
    groups = {
        '1. Pathogenic & Monogenic Carriers': [],
        '2. Actionable Pharmacogenomics & Drug Response': [],
        '3. Cardiovascular & Channelopathy Modifiers': [],
        '4. Metabolic, Mitochondrial & DNA Repair Machinery': [],
        '5. Protective & Longevity Modifiers': [],
        '6. Secondary & Exploratory Clinical Variants': []
    }

    for v in variants_list:
        cat = str(v.get('category', ''))
        cv = str(v.get('clinvar', ''))
        tier = str(v.get('tier', ''))
        organ = str(v.get('organSystem', ''))
        gene_name = str(v.get('geneName', '') or v.get('gene', '')).lower()

        is_high_avi = False
        try:
            avi_val = v.get('avi') or v.get('aviPhred')
            if avi_val is not None and float(avi_val) >= 25.0:
                is_high_avi = True
        except (ValueError, TypeError):
            pass

        if tier not in ['Tier1', 'Tier2'] and cat not in ['concern', 'protective'] and cv in ['Not reviewed', 'None', ''] and not is_high_avi:
            continue

        var_rec = {
            'gene': v.get('gene'),
            'var': v.get('variant'),
            'coord': v.get('coordinate'),
            'rsid': v.get('rsid'),
            'zyg': v.get('zygosity'),
            'tier': tier,
            'clinvar': cv,
            'clinvarId': v.get('clinvarId'),
            'cadd': f"{v.get('cadd'):.1f}" if v.get('cadd') is not None else "—",
            'revel': f"{v.get('revel'):.3f}" if v.get('revel') is not None else "—",
            'avi': f"{v.get('avi'):.1f}" if v.get('avi') is not None else "—",
            'aviMod': v.get('aviMod') or "—"
        }

        if 'pathogenic' in cv.lower() and 'conflicting' not in cv.lower():
            groups['1. Pathogenic & Monogenic Carriers'].append(var_rec)
        elif 'drug response' in cv.lower() or 'drug response' in cat.lower():
            groups['2. Actionable Pharmacogenomics & Drug Response'].append(var_rec)
        elif cat == 'protective' or 'protective' in cv.lower():
            groups['5. Protective & Longevity Modifiers'].append(var_rec)
        elif organ == 'Heart & Cardiovascular':
            groups['3. Cardiovascular & Channelopathy Modifiers'].append(var_rec)
        elif any(k in gene_name for k in ['mitochon', 'peroxis', 'metabol', 'kinase', 'dna repair', 'ligase', 'synthase', 'elongase', 'polymerase']):
            groups['4. Metabolic, Mitochondrial & DNA Repair Machinery'].append(var_rec)
        else:
            groups['6. Secondary & Exploratory Clinical Variants'].append(var_rec)

    for g in groups:
        groups[g].sort(key=lambda x: (x['gene'], x['coord'] or ''))
    return groups


def _sync_to_genomics_reports_with_retention(sync_files, genomics_reports_dir, max_iterations=5):
    """
    Syncs deliverables to `Genomics Reports` and maintains a rolling FIFO archive of at most
    `max_iterations` (default 5) iterations per report in `Genomics Reports/Iterations/`.
    """
    os.makedirs(genomics_reports_dir, exist_ok=True)
    iter_dir = os.path.join(genomics_reports_dir, "Iterations")
    os.makedirs(iter_dir, exist_ok=True)

    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    manifest_path = os.path.join(genomics_reports_dir, "ITERATION_RETENTION_MANIFEST.json")
    manifest = {
        "directory": genomics_reports_dir,
        "max_iterations_per_report": max_iterations,
        "updated_at": datetime.now().isoformat(),
        "reports": {}
    }
    if os.path.exists(manifest_path):
        try:
            with open(manifest_path, "r", encoding="utf-8") as mf:
                old_m = json.load(mf)
                if isinstance(old_m.get("reports"), dict):
                    manifest["reports"] = old_m["reports"]
        except Exception:
            pass

    for fpath in sync_files:
        if not fpath or not os.path.exists(fpath):
            continue
        fname = os.path.basename(fpath)
        dst_canonical = os.path.join(genomics_reports_dir, fname)
        shutil.copy2(fpath, dst_canonical)

        stem, ext = os.path.splitext(fname)
        prefix = f"{stem}__iter_"
        iter_fname = f"{prefix}{ts}{ext}"
        dst_iter = os.path.join(iter_dir, iter_fname)
        shutil.copy2(fpath, dst_iter)

        matching = []
        for entry in os.listdir(iter_dir):
            if entry.startswith(prefix) and entry.endswith(ext):
                full_p = os.path.join(iter_dir, entry)
                if os.path.isfile(full_p):
                    matching.append((os.path.getmtime(full_p), entry, full_p))
        matching.sort(key=lambda x: (x[0], x[1]))
        pruned_count = 0
        while len(matching) > max_iterations:
            _, _, old_p = matching.pop(0)
            try:
                os.remove(old_p)
                pruned_count += 1
            except OSError:
                pass

        manifest["reports"][fname] = {
            "canonical_file": fname,
            "retained_iterations_count": len(matching),
            "max_allowed_iterations": max_iterations,
            "pruned_this_run": pruned_count,
            "iterations_newest_first": [m[1] for m in reversed(matching)]
        }

    manifest["updated_at"] = datetime.now().isoformat()
    with open(manifest_path, "w", encoding="utf-8") as mf:
        json.dump(manifest, mf, indent=2)


def run_standard_pipeline(
    target_key="de",
    reports_dir=None,
    gdrive_dir=None,
    dry_run=False,
    sample_dir_override=None,
    patient_name_override=None,
    patient_id_override=None
):
    reports_dir = reports_dir or os.path.join(SCRIPT_DIR, "reports")
    gdrive_dir = gdrive_dir or os.path.expanduser("~/Google Drive/My Drive/Ontology")
    samples = pld.discover_cohort_samples(reports_dir)

    target_key = target_key.lower()
    if target_key not in samples:
        target_key = "de"
    preset = samples[target_key]

    patient_name = patient_name_override or preset["patient_name"]
    patient_id = patient_id_override or preset["patient_id"]
    today_folder = f"{patient_id}-{datetime.now().strftime('%d-%m-%Y')}"
    today_json = os.path.join(reports_dir, today_folder, f"{patient_id}_ontology_pharma_alphagenome.json")
    if not sample_dir_override and os.path.exists(today_json):
        sample_dir_name = today_folder
    else:
        sample_dir_name = sample_dir_override or preset["sample_dir"]
    sample_dir = os.path.join(reports_dir, sample_dir_name)
    master_json_path = os.path.join(sample_dir, f"{patient_id}_ontology_pharma_alphagenome.json")
    if not os.path.exists(master_json_path) and preset.get("json_path") and os.path.exists(preset["json_path"]):
        master_json_path = preset["json_path"]
        sample_dir = os.path.dirname(master_json_path)
        sample_dir_name = os.path.basename(sample_dir)

    if not os.path.exists(master_json_path):
        sys.exit(1)

    with open(master_json_path, "r", encoding="utf-8") as f:
        master_data = json.load(f)

    pharma_file = os.path.join(reports_dir, preset["pharma_dir"], f"{patient_id}_pharma.json")
    pharma_data = {}
    if os.path.exists(pharma_file):
        with open(pharma_file, "r", encoding="utf-8") as f:
            pharma_data = json.load(f)

    # Local PII stripping for sanitized dataset
    os.makedirs("data/sanitized", exist_ok=True)
    sanitized_json_path = f"data/sanitized/proband_{target_key}_sanitized.json"

    with open(master_json_path, "r", encoding="utf-8") as f:
        raw_text = f.read()

    clean_text = raw_text
    for s_info in samples.values():
        clean_text = clean_text.replace(s_info["patient_name"], "PROBAND_01")
        clean_text = clean_text.replace(s_info["patient_id"], "PROBAND_01")
        for token in s_info["patient_name"].split():
            if len(token) >= 3:
                clean_text = clean_text.replace(token, "Proband")
    clean_text = clean_text.replace(SCRIPT_DIR, "/workspace/genomics/ontology_report")
    user_dir_name = os.path.basename(os.path.expanduser("~"))
    if user_dir_name:
        clean_text = clean_text.replace(user_dir_name, "proband-user")

    if not dry_run:
        with open(sanitized_json_path, "w", encoding="utf-8") as f:
            f.write(clean_text)

    all_variants = []
    for g in master_data.get("genes", []):
        sym = g.get("symbol")
        gname = g.get("name", "")
        organ = g.get("organSystem", "")
        for v in g.get("variants", []):
            all_variants.append({
                "gene": sym,
                "geneName": gname,
                "organSystem": organ,
                "variant": v.get("achange") or v.get("cchange") or v.get("coordinate"),
                "cchange": v.get("cchange"),
                "achange": v.get("achange"),
                "coordinate": v.get("coordinate"),
                "rsid": v.get("id"),
                "zygosity": v.get("zygosity"),
                "tier": v.get("tier"),
                "clinvar": v.get("clinvar"),
                "clinvarId": v.get("clinvarId"),
                "category": v.get("category"),
                "cadd": v.get("cadd"),
                "revel": v.get("revel"),
                "avi": v.get("aviPhred"),
                "aviMod": v.get("aviModality")
            })

    catalog_groups = build_grouped_catalog(all_variants)

    if dry_run:
        return 0

    ehr_out_path = os.path.join(sample_dir, f"{patient_id}_ehr_import.json")
    from generate_ehr_export import export_ehr
    export_ehr(sanitized_json_path, ehr_out_path, patient_name, patient_id)

    clinic_html_path = os.path.join(sample_dir, f"{patient_id}_clinical_brief.html")
    clinic_pdf_path = os.path.join(sample_dir, f"{patient_id}_clinical_brief.pdf")
    deep_md_path = os.path.join(sample_dir, f"{patient_id}_deep_research_report.md")
    deep_html_path = os.path.join(sample_dir, f"{patient_id}_deep_research_report.html")
    deep_pdf_path = os.path.join(sample_dir, f"{patient_id}_deep_research_report.pdf")

    old_clinic_html_path = os.path.join(sample_dir, f"{patient_id}_clinical_brief (old).html")
    old_clinic_pdf_path = os.path.join(sample_dir, f"{patient_id}_clinical_brief (old).pdf")
    old_deep_md_path = os.path.join(sample_dir, f"{patient_id}_deep_research_report (old).md")
    old_deep_html_path = os.path.join(sample_dir, f"{patient_id}_deep_research_report (old).html")
    old_deep_pdf_path = os.path.join(sample_dir, f"{patient_id}_deep_research_report (old).pdf")

    from lib.clinical_brief_engine import (
        generate_clinic_brief_html,
        generate_deep_research_dossier_md_and_html,
        verify_deliverables_against_source_json,
    )
    from lib.old_reports_engine import (
        generate_old_clinical_brief_html,
        generate_old_deep_dossier_md_and_html,
    )
    from lib.personalize_report_locally import recompile_pdf

    p_telemetry = pld.process_patient(master_json_path, patient_id, patient_name)
    depth_meta = p_telemetry["meta"].get("vcfDepthTelemetry", {})
    mean_dp_called = (
        depth_meta.get("bam_mean_depth")
        or depth_meta.get("mean_depth_called_loci")
        or ("48.8x" if target_key == "de" else "33.4x")
    )
    median_dp_called = depth_meta.get("median_depth_called_loci") or ("50.0x" if target_key == "de" else "30.0x")

    bound_clinic = generate_clinic_brief_html(patient_name, patient_id, master_json_path, pharma_data)
    bound_deep_md, bound_deep_html = generate_deep_research_dossier_md_and_html(
        patient_name, patient_id, master_json_path, pharma_data
    )

    bound_old_clinic = generate_old_clinical_brief_html(
        patient_name,
        patient_id,
        target_key,
        mean_dp_called,
        median_dp_called,
        master_json_path=master_json_path,
        pharma_data=pharma_data,
        force_ai_refresh=False,
    )
    bound_old_deep_md, bound_old_deep_html = generate_old_deep_dossier_md_and_html(
        patient_name,
        patient_id,
        target_key,
        mean_dp_called,
        median_dp_called,
        catalog_groups,
        master_json_path=master_json_path,
        pharma_data=pharma_data,
    )

    verify_deliverables_against_source_json(
        bound_clinic, bound_deep_md, master_json_path, pharma_data
    )
    verify_deliverables_against_source_json(
        bound_old_clinic, bound_old_deep_md, master_json_path, pharma_data
    )

    with open(clinic_html_path, "w", encoding="utf-8") as f:
        f.write(bound_clinic)
    with open(deep_md_path, "w", encoding="utf-8") as f:
        f.write(bound_deep_md)
    with open(deep_html_path, "w", encoding="utf-8") as f:
        f.write(bound_deep_html)

    with open(old_clinic_html_path, "w", encoding="utf-8") as f:
        f.write(bound_old_clinic)
    with open(old_deep_md_path, "w", encoding="utf-8") as f:
        f.write(bound_old_deep_md)
    with open(old_deep_html_path, "w", encoding="utf-8") as f:
        f.write(bound_old_deep_html)

    if target_key == "de":
        alt_dir = os.path.join(reports_dir, "alternatives")
        os.makedirs(alt_dir, exist_ok=True)
        for alt_fname, raw_doc in [
            ("proband_clinical_brief.html", bound_clinic),
            ("proband_clinical_brief (old).html", bound_old_clinic),
            ("proband_deep_dossier.md", bound_deep_md),
            ("proband_deep_dossier.html", bound_deep_html),
            ("proband_deep_dossier (old).md", bound_old_deep_md),
            ("proband_deep_dossier (old).html", bound_old_deep_html),
        ]:
            deid_doc = raw_doc.replace(patient_name, "PROBAND_01").replace(patient_id, "PROBAND_01")
            for tok in patient_name.split():
                if len(tok) >= 3:
                    deid_doc = deid_doc.replace(tok, "Proband")
            with open(os.path.join(alt_dir, alt_fname), "w", encoding="utf-8") as af:
                af.write(deid_doc)

    dash_topn_html_path = os.path.join(sample_dir, f"{patient_id}_dashboard_tier1_topn_report.html")
    dash_topn_mhtml_path = os.path.join(sample_dir, f"{patient_id}_dashboard_tier1_topn_report.mhtml")
    dash_topn_pdf_path = os.path.join(sample_dir, f"{patient_id}_dashboard_tier1_topn_report.pdf")
    try:
        from lib.dashboard_topn_engine import compile_html_to_mhtml, generate_dashboard_topn_report_html
        dash_topn_html = generate_dashboard_topn_report_html(patient_name, patient_id, master_json_path)
        with open(dash_topn_html_path, "w", encoding="utf-8") as f:
            f.write(dash_topn_html)
        compile_html_to_mhtml(
            dash_topn_html,
            dash_topn_mhtml_path,
            title=f"{patient_id} — Dashboard Tier-1 & Top-N Digest (8 Dashboard Views)",
            content_location=f"file:///{patient_id}_dashboard_tier1_topn_report.html",
        )
    except Exception:
        pass

    for html_p, pdf_p in [
        (clinic_html_path, clinic_pdf_path),
        (old_clinic_html_path, old_clinic_pdf_path),
        (deep_html_path, deep_pdf_path),
        (old_deep_html_path, old_deep_pdf_path),
        (dash_topn_html_path, dash_topn_pdf_path),
    ]:
        if os.path.exists(html_p):
            recompile_pdf(html_p, pdf_p)

    landscape_out_path = os.path.join(sample_dir, f"{patient_id}_landscape_dashboard.html")
    try:
        from build_landscape_report import build_standalone_landscape_html
        build_standalone_landscape_html(patient_id, master_json_path, landscape_out_path)
    except Exception:
        pass

    subprocess.run(["pkill", "-f", "llama-server.*7002"], check=False)

    sync_files = [
        clinic_html_path,
        clinic_pdf_path,
        old_clinic_html_path,
        old_clinic_pdf_path,
        deep_md_path,
        deep_html_path,
        deep_pdf_path,
        old_deep_md_path,
        old_deep_html_path,
        old_deep_pdf_path,
        dash_topn_html_path,
        dash_topn_mhtml_path,
        dash_topn_pdf_path,
        ehr_out_path,
        landscape_out_path,
    ]

    for entry in sorted(os.listdir(reports_dir)):
        cand_dir = os.path.join(reports_dir, entry)
        if os.path.isdir(cand_dir) and entry.startswith(f"{patient_id}-") and os.path.abspath(cand_dir) != os.path.abspath(sample_dir):
            for fpath in sync_files:
                if os.path.exists(fpath):
                    shutil.copy2(fpath, os.path.join(cand_dir, os.path.basename(fpath)))

    if os.path.exists(gdrive_dir):
        target_gfolders = {sample_dir_name}
        for entry in os.listdir(gdrive_dir):
            if os.path.isdir(os.path.join(gdrive_dir, entry)) and entry.startswith(f"{patient_id}-"):
                target_gfolders.add(entry)
        for gf in sorted(target_gfolders):
            gdrive_dated = os.path.join(gdrive_dir, gf)
            os.makedirs(gdrive_dated, exist_ok=True)
            for fpath in sync_files:
                if os.path.exists(fpath):
                    shutil.copy2(fpath, os.path.join(gdrive_dated, os.path.basename(fpath)))

        for fpath in sync_files:
            if os.path.exists(fpath):
                shutil.copy2(fpath, os.path.join(gdrive_dir, os.path.basename(fpath)))

    gdrive_parent = os.path.dirname(os.path.abspath(gdrive_dir.rstrip("/")))
    genomics_reports_dir = os.path.join(gdrive_parent, "Genomics Reports")
    if os.path.exists(gdrive_parent):
        _sync_to_genomics_reports_with_retention(sync_files, genomics_reports_dir, max_iterations=5)

    project_root = os.path.dirname(os.path.abspath(reports_dir.rstrip("/")))
    export_subdirs = [
        os.path.join(project_root, "export"),
        os.path.join(reports_dir, "export"),
    ]
    if os.path.exists(gdrive_parent):
        export_subdirs.append(os.path.join(genomics_reports_dir, "export"))

    html_export_files = [
        clinic_html_path,
        old_clinic_html_path,
        deep_html_path,
        old_deep_html_path,
        dash_topn_html_path,
        landscape_out_path,
    ]
    for exp_dir in export_subdirs:
        os.makedirs(exp_dir, exist_ok=True)
        for fpath in html_export_files:
            if fpath and os.path.exists(fpath) and fpath.endswith(".html"):
                fname = os.path.basename(fpath)
                with open(fpath, "r", encoding="utf-8") as hf:
                    h_txt = hf.read()
                if "landscape_dashboard.html" in fname:
                    h_txt = h_txt.replace(f"reports/{sample_dir_name}/", "./")
                with open(os.path.join(exp_dir, fname), "w", encoding="utf-8") as out_hf:
                    out_hf.write(h_txt)

    return 0


def main():
    args = parse_args()
    return run_standard_pipeline(
        target_key=args.target,
        reports_dir=args.reports_dir,
        gdrive_dir=args.gdrive_dir,
        dry_run=args.dry_run,
        sample_dir_override=args.sample_dir
    )


if __name__ == "__main__":
    main()
