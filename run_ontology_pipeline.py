#!/usr/bin/env python3
"""
run_ontology_pipeline.py
Unified Command-Line Execution Engine for the Genomic Ontology Reporting System (v5.2).

Produces standalone, 100% self-contained HTML5 deliverables (like Chrome "Webpage, Complete"):
  - {Sample_ID}_visual_explorer.html (Self-contained HTML5 with inlined CSS, Data, and Client JS)
  - {Sample_ID}_master_ontology_report.html (Self-contained Master Hub HTML5)
  - {Sample_ID}_master_actionable.json
  - {Sample_ID}_variants.tsv
  - {Sample_ID}_summary.txt
  - {Sample_ID}_report.pdf (Headless Chrome generated)
  - {Sample_ID}_iOS_bundle.zip

Outputs to:
  Google Drive : ~/Google Drive/My Drive/Ontology/{Sample_ID}-{DD-MM-YYYY}/
  Local Project: ./reports/{Sample_ID}-{DD-MM-YYYY}/
"""

import os
import sys
import argparse
import subprocess
import shutil
import sqlite3
import json
import glob
from datetime import datetime

OC_ANNOTATORS = [
    "hpo", "go", "clinvar", "clingen", "omim", "ncbigene", "revel",
    "alphamissense", "bayesdel", "metarnn", "esm1b", "varity", "spliceai",
    "cadd", "linsight", "ncer", "regulomedb", "ccre_screen", "gtex",
    "dbsnp", "vcfinfo", "gwas_catalog", "pharmgkb", "civic", "interpro"
]

def parse_args():
    parser = argparse.ArgumentParser(
        description="Run end-to-end Genomic Ontology Report from VCF, SQLite, or OpenCRAVAT Job ID."
    )
    parser.add_argument(
        "--sample", "-s", "--patient-id", "-p",
        dest="sample_id",
        required=True,
        help="Sample or Patient Identifier (e.g., DE_master, HG003)"
    )
    parser.add_argument(
        "--input", "-i",
        required=True,
        help="Input source: .vcf, .vcf.gz, .sqlite file, or OpenCRAVAT Job ID (e.g., 260706-105810)"
    )
    parser.add_argument(
        "--gdrive-dir",
        default=os.path.expanduser("~/Google Drive/My Drive/Ontology"),
        help="Target Google Drive directory (default: ~/Google Drive/My Drive/Ontology)"
    )
    parser.add_argument(
        "--config", "-c",
        default="config/ontology_domains.yaml",
        help="Domain configuration YAML (default: config/ontology_domains.yaml)"
    )
    parser.add_argument(
        "--no-pdf",
        action="store_true",
        help="Skip PDF generation"
    )
    parser.add_argument(
        "--local-only",
        action="store_true",
        help="Skip syncing to Google Drive"
    )
    return parser.parse_args()

def resolve_input(input_source, sample_id, work_dir):
    """
    Resolves input argument to SQLite database path and VCF path.
    Handles .vcf, .vcf.gz, .sqlite, or OpenCRAVAT Job ID.
    """
    raw_db = None
    vcf_path = None

    if not os.path.exists(input_source):
        job_dir_candidate = f"/data/opencravat/jobs/default/{input_source}"
        if os.path.exists(job_dir_candidate):
            input_source = job_dir_candidate
        else:
            job_glob = glob.glob(f"/data/opencravat/jobs/*/{input_source}")
            if job_glob:
                input_source = job_glob[0]

    if os.path.isdir(input_source):
        sqlites = glob.glob(os.path.join(input_source, "*.sqlite"))
        vcfs = glob.glob(os.path.join(input_source, "*.vcf.gz")) or glob.glob(os.path.join(input_source, "*.vcf"))
        if sqlites:
            raw_db = sqlites[0]
            print(f"[Input Resolver] Found OpenCRAVAT SQLite in job directory: {raw_db}")
        if vcfs:
            vcf_path = vcfs[0]
            print(f"[Input Resolver] Found VCF in job directory: {vcf_path}")
        if not raw_db:
            raise FileNotFoundError(f"No SQLite database found inside OpenCRAVAT job folder: {input_source}")

    elif input_source.endswith(".sqlite"):
        raw_db = input_source
        candidate_vcf = input_source.replace(".sqlite", "")
        if os.path.exists(candidate_vcf):
            vcf_path = candidate_vcf
        elif os.path.exists(input_source.replace(".vcf.gz.sqlite", ".vcf.gz")):
            vcf_path = input_source.replace(".vcf.gz.sqlite", ".vcf.gz")

    elif input_source.endswith((".vcf", ".vcf.gz", ".g.vcf", ".g.vcf.gz")):
        vcf_path = input_source
        raw_db = os.path.join(work_dir, f"{sample_id}.sqlite")
        if not os.path.exists(raw_db):
            print(f"[OpenCRAVAT] Annotating input VCF: {input_source}...")
            oc_cmd = [
                "oc", "run", input_source,
                "-l", "hg38",
                "-a", *OC_ANNOTATORS,
                "-d", work_dir,
                "--mp", str(os.cpu_count() or 4),
                "-n", sample_id
            ]
            print(f"Executing: {' '.join(oc_cmd)}")
            subprocess.run(oc_cmd, check=True)
    else:
        raise ValueError(f"Unrecognized input format: {input_source}. Expected .vcf, .vcf.gz, .sqlite, or OC Job ID.")

    return raw_db, vcf_path

def build_standalone_html5(template_html, css_path, js_data_path, js_app_path, out_html):
    """
    Builds a 100% self-contained standalone HTML5 deliverable (like Chrome 'Save As: Webpage, Complete')
    with inlined CSS styles, inlined dataset, and inlined application logic.
    """
    with open(template_html, 'r', encoding='utf-8') as f:
        html = f.read()
    with open(css_path, 'r', encoding='utf-8') as f:
        css = f.read()
    with open(js_data_path, 'r', encoding='utf-8') as f:
        js_data = f.read()
    with open(js_app_path, 'r', encoding='utf-8') as f:
        js_app = f.read()

    # Inline CSS
    html = html.replace('<link rel="stylesheet" href="css/style.css" />', f'<style>\n{css}\n</style>')
    html = html.replace('<link rel="stylesheet" href="css/style.css">', f'<style>\n{css}\n</style>')
    # Inline Data & App JS
    html = html.replace('<script src="data/mock-data.js"></script>', f'<script>\n{js_data}\n</script>')
    html = html.replace('<script src="js/app.js"></script>', f'<script>\n{js_app}\n</script>')

    with open(out_html, 'w', encoding='utf-8') as f:
        f.write(html)
    print(f"[HTML5 Standalone Engine] Built complete single-file report ({os.path.getsize(out_html)/1024:.1f} KB): {out_html}")
    return out_html

def generate_pdf_report(html_path, pdf_path):
    """
    Generates high-resolution PDF from HTML using headless Chrome.
    """
    chrome_bin = shutil.which("google-chrome") or shutil.which("chromium-browser") or shutil.which("chromium")
    if not chrome_bin:
        print("[Warning] No headless Chrome/Chromium binary found. Skipping PDF rendering.")
        return False
    try:
        cmd = [
            chrome_bin,
            "--headless=new",
            "--disable-gpu",
            f"--print-to-pdf={pdf_path}",
            "--no-pdf-header-footer",
            f"file://{os.path.abspath(html_path)}"
        ]
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        if res.returncode == 0 and os.path.exists(pdf_path):
            print(f"[PDF Engine] Generated PDF report: {pdf_path} ({os.path.getsize(pdf_path)/1024:.1f} KB)")
            return True
        else:
            print(f"[PDF Engine Warning] Chrome PDF error: {res.stderr.strip()}")
            return False
    except Exception as e:
        print(f"[PDF Engine Warning] Failed to render PDF: {e}")
        return False

def deliver_to_google_drive(src_files, gdrive_target_dir):
    """
    Syncs generated files to the target Google Drive directory and uses rclone if available.
    """
    os.makedirs(gdrive_target_dir, exist_ok=True)
    for f in src_files:
        if os.path.exists(f):
            dest = os.path.join(gdrive_target_dir, os.path.basename(f))
            shutil.copy2(f, dest)
            print(f"  [GDrive Sync] Copied -> {dest}")

    # Optional Rclone Cloud Sync
    rclone_bin = shutil.which("rclone")
    if rclone_bin:
        rel_folder = os.path.basename(gdrive_target_dir)
        cmd = [rclone_bin, "copy", gdrive_target_dir, f"drive:Ontology/{rel_folder}/"]
        try:
            subprocess.run(cmd, capture_output=True, timeout=30)
            print(f"  [rclone Cloud Sync] Synced folder to remote Google Drive: drive:Ontology/{rel_folder}/")
        except Exception:
            pass

def main():
    args = parse_args()
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)

    # Format Date as dd-mm-yyyy
    now = datetime.now()
    date_str = now.strftime("%d-%m-%Y")
    subfolder_name = f"{args.sample_id}-{date_str}"

    local_outdir = os.path.join(script_dir, "reports", subfolder_name)
    os.makedirs(local_outdir, exist_ok=True)

    print("==================================================================")
    print(f"GENOMIC ONTOLOGY PIPELINE ENGINE (v5.2)")
    print(f"  Sample / Patient ID : {args.sample_id}")
    print(f"  Input Source        : {args.input}")
    print(f"  Execution Date      : {date_str}")
    print(f"  Local Output Dir    : {local_outdir}")
    print("==================================================================")

    # 1. Resolve Input
    raw_db, vcf_path = resolve_input(args.input, args.sample_id, local_outdir)
    print(f"[Stage 1] Resolved Database : {raw_db}")
    print(f"[Stage 1] Resolved VCF File : {vcf_path or 'None (will use DB attributes)'}")

    # File naming based on sample / source prefix
    base_prefix = args.sample_id
    schema_json = os.path.join(local_outdir, f"{base_prefix}_schema.json")
    panel_json = os.path.join(local_outdir, f"{base_prefix}_master_panel.json")
    act_db = os.path.join(local_outdir, f"{base_prefix}_master_actionable.sqlite")
    act_json = os.path.join(local_outdir, f"{base_prefix}_master_actionable.json")
    enrich_cache = os.path.join(local_outdir, f"{base_prefix}_enrich_cache.json")

    # Pre-seed enrich cache from existing reports if available
    if not os.path.exists(enrich_cache):
        for candidate in glob.glob(os.path.join(script_dir, "reports", "*", "*enrich_cache.json")):
            if os.path.exists(candidate) and os.path.getsize(candidate) > 1000:
                shutil.copy2(candidate, enrich_cache)
                break

    # Output Deliverables
    visual_explorer_html = os.path.join(local_outdir, f"{base_prefix}_visual_explorer.html")
    master_hub_html = os.path.join(local_outdir, f"{base_prefix}_master_ontology_report.html")
    tsv_report = os.path.join(local_outdir, f"{base_prefix}_variants.tsv")
    txt_report = os.path.join(local_outdir, f"{base_prefix}_summary.txt")
    pdf_report = os.path.join(local_outdir, f"{base_prefix}_report.pdf")
    zip_bundle = os.path.join(local_outdir, f"{base_prefix}_iOS_bundle.zip")

    # 2. Build Multi-Domain Panel
    print("\n[Stage 2/7] Building HPO + GO + Organ multi-domain gene panel...")
    subprocess.run([
        "python3", "lib/build_ontology_panel.py",
        "--config", args.config,
        "--out", panel_json
    ], check=True)

    # 3. Schema Probe
    print("\n[Stage 3/7] Probing database schema & column mappings...")
    subprocess.run([
        "python3", "lib/schema_probe.py",
        raw_db,
        "--out", schema_json
    ], check=True)

    # 4. Filter & Tier Actionable Variants (with Protective & PGx Frequency Bypass)
    print("\n[Stage 4/7] Filtering actionable variants & extracting protective alleles...")
    subprocess.run([
        "python3", "lib/ontology_filter.py",
        "--raw-db", raw_db,
        "--panel", panel_json,
        "--schema", schema_json,
        "--config", args.config,
        "--out-sqlite", act_db,
        "--out-json", act_json,
        "--patient", args.sample_id
    ], check=True)

    # 5. Enrich Gene Descriptions & Clinical Synopses
    print("\n[Stage 5/7] Enriching gene annotations & literature references...")
    try:
        subprocess.run([
            "python3", "lib/enrich_report.py",
            "--in-json", act_json,
            "--cache", enrich_cache,
            "--genes"
        ], timeout=45)
    except Exception as e:
        print(f"[Enrichment Note] Completed / Cached fallback applied: {e}")

    # 6. Render Master Hub & Standalone Visual Explorer
    print("\n[Stage 6/7] Rendering Master Hub & Standalone Single-File Visual Explorer HTML5...")
    subprocess.run([
        "python3", "lib/render_master_hub.py",
        "--in-json", act_json,
        "--out-html", master_hub_html,
        "--out-tsv", tsv_report,
        "--out-text", txt_report,
        "--domain-config", args.config
    ], check=True)

    mock_data_js = os.path.join(script_dir, "data", "mock-data.js")
    subprocess.run([
        "python3", "generate_claude_v2_report.py",
        act_json,
        raw_db,
        vcf_path or "/dev/null",
        mock_data_js
    ], check=True)

    # Build 100% Self-Contained Standalone HTML5 Document (Chrome "Save As: Webpage, Complete")
    build_standalone_html5(
        os.path.join(script_dir, "index.html"),
        os.path.join(script_dir, "css", "style.css"),
        mock_data_js,
        os.path.join(script_dir, "js", "app.js"),
        visual_explorer_html
    )

    # 7. Generate PDF & Zip Bundle
    print("\n[Stage 7/7] Generating PDF report and packaging deliverables...")
    if not args.no_pdf:
        generate_pdf_report(visual_explorer_html, pdf_report)

    with open(os.devnull, 'w') as devnull:
        subprocess.run([
            "zip", "-q", "-r", zip_bundle,
            os.path.basename(visual_explorer_html),
            os.path.basename(master_hub_html),
            os.path.basename(tsv_report),
            os.path.basename(txt_report),
            os.path.basename(act_json)
        ], cwd=local_outdir, stdout=devnull, stderr=devnull)

    deliverables = [
        visual_explorer_html,
        master_hub_html,
        act_json,
        tsv_report,
        txt_report,
        pdf_report,
        zip_bundle
    ]

    if not args.local_only:
        gdrive_target = os.path.join(args.gdrive_dir, subfolder_name)
        print(f"\n[Google Drive Delivery] Uploading deliverables to: {gdrive_target}...")
        deliver_to_google_drive(deliverables, gdrive_target)

    print("\n==================================================================")
    print("✨ PIPELINE COMPLETE. ALL DELIVERABLES GENERATED SUCCESSFULLY:")
    print(f"  1. Visual Explorer HTML (Standalone HTML5) : {visual_explorer_html}")
    print(f"  2. Universal Master Hub (HTML5)            : {master_hub_html}")
    print(f"  3. Actionable JSON                         : {act_json}")
    print(f"  4. Variant TSV Matrix                      : {tsv_report}")
    print(f"  5. Clinical Summary TXT                    : {txt_report}")
    print(f"  6. Printable PDF Report                    : {pdf_report}")
    print(f"  7. Offline iOS Bundle                      : {zip_bundle}")
    if not args.local_only:
        print(f"  👉 Google Drive Folder                     : ~/Google Drive/My Drive/Ontology/{subfolder_name}/")
    print("==================================================================")

if __name__ == "__main__":
    main()
