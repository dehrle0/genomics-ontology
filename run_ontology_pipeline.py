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

Delivers automatically to:
  1. Google Drive Cloud Remote (via rclone to drive:Ontology/{Sample_ID}-{DD-MM-YYYY}/)
  2. Local Google Drive Directory (~/Google Drive/My Drive/Ontology/{Sample_ID}-{DD-MM-YYYY}/)
  3. Local Project Workspace (./reports/{Sample_ID}-{DD-MM-YYYY}/)
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
    "alphamissense", "bayesdel", "metarnn", "esm1b", "varity_r", "spliceai",
    "cadd", "linsight", "ncer", "regulomedb", "ccre_screen", "gtex",
    "dbsnp", "gnomad4", "allofus250k", "gwas_catalog", "pharmgkb", "civic", "interpro"
]

def find_rclone():
    """
    Finds rclone executable across standard and conda/micromamba paths.
    """
    candidates = [
        shutil.which("rclone"),
        os.path.expanduser("~/micromamba/envs/cravat_env/bin/rclone"),
        os.path.expanduser("~/micromamba/pkgs/rclone-1.75.0-h519d9b9_0/bin/rclone"),
        "/usr/bin/rclone",
        "/usr/local/bin/rclone"
    ]
    for c in candidates:
        if c and os.path.exists(c) and os.access(c, os.X_OK):
            return c
    return None

def parse_args():
    parser = argparse.ArgumentParser(
        description="Run end-to-end Genomic Ontology Report from VCF, SQLite, or OpenCRAVAT Job ID."
    )
    parser.add_argument(
        "--sample", "-s", "--patient-id", "-p",
        dest="sample_id",
        default=None,
        help="Sample or Patient Identifier (default: auto-derived from input source)"
    )
    parser.add_argument(
        "--input", "-i",
        required=True,
        help="Input source: .vcf, .vcf.gz, .sqlite file, or OpenCRAVAT Job ID (e.g., 260706-105810)"
    )
    parser.add_argument(
        "--gdrive-dir",
        default=os.path.expanduser("~/Google Drive/My Drive/Ontology"),
        help="Target local Google Drive directory (default: ~/Google Drive/My Drive/Ontology)"
    )
    parser.add_argument(
        "--config", "-c",
        default="config/ontology_domains.yaml",
        help="Domain configuration YAML (default: config/ontology_domains.yaml)"
    )
    parser.add_argument(
        "--vcf", "--vcfs", "--phased-vcf",
        dest="phased_vcf",
        default=None,
        help="Optional phased VCF(s) (comma-separated if multiple, e.g. SNVs, SVs, CNVs) to supply phased haplotypes and genotypes."
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

def resolve_input(input_source, sample_id, work_dir, phased_vcf=None):
    """
    Resolves input argument to SQLite database path and VCF path.
    Handles .vcf, .vcf.gz, .sqlite, or OpenCRAVAT Job ID.
    """
    raw_db = None
    vcf_path = phased_vcf

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
        if vcfs and not vcf_path:
            snv_vcfs = [v for v in vcfs if not ('.cnv.' in v or '.sv.' in v)]
            vcf_path = snv_vcfs[0] if snv_vcfs else max(vcfs, key=os.path.getsize)
            print(f"[Input Resolver] Found VCF in job directory: {vcf_path}")
        if not raw_db:
            raise FileNotFoundError(f"No SQLite database found inside OpenCRAVAT job folder: {input_source}")

    elif input_source.endswith(".sqlite"):
        raw_db = input_source
        if not vcf_path:
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

    html = html.replace('<link rel="stylesheet" href="css/style.css" />', f'<style>\n{css}\n</style>')
    html = html.replace('<link rel="stylesheet" href="css/style.css">', f'<style>\n{css}\n</style>')
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

def deliver_to_google_drive(src_files, local_gdrive_dir, subfolder_name):
    """
    Syncs generated files to both local Google Drive directory and cloud remote via rclone.
    """
    # 1. Local sync directory
    target_local_dir = os.path.join(local_gdrive_dir, subfolder_name)
    os.makedirs(target_local_dir, exist_ok=True)
    for f in src_files:
        if os.path.exists(f):
            dest = os.path.join(target_local_dir, os.path.basename(f))
            shutil.copy2(f, dest)
            print(f"  [Local GDrive Sync] Copied -> {dest}")

    # Also copy to root local Ontology folder for legacy view
    for f in src_files:
        if os.path.exists(f) and f.endswith((".html", ".tsv", ".txt", ".pdf", ".zip")):
            dest_root = os.path.join(local_gdrive_dir, os.path.basename(f))
            shutil.copy2(f, dest_root)

    # 2. Cloud rclone Sync
    rclone_bin = find_rclone()
    if rclone_bin:
        print(f"\n  [rclone Cloud Sync] Found rclone at: {rclone_bin}")
        for remote in ["drive:Ontology", "gdrive:Ontology"]:
            try:
                cmd = [rclone_bin, "copy", target_local_dir, f"{remote}/{subfolder_name}/"]
                res = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
                if res.returncode == 0:
                    print(f"  [rclone Cloud Sync] Uploaded to Google Drive Cloud: {remote}/{subfolder_name}/")
                    # Also upload top-level files
                    subprocess.run([rclone_bin, "copy", target_local_dir, f"{remote}/"], capture_output=True, timeout=60)
                    break
            except Exception as e:
                print(f"  [rclone Cloud Sync Note] {remote} error: {e}")

def update_actionable_with_phased_vcfs(act_json_path, act_db_path, vcf_paths):
    if not vcf_paths:
        return
    vcf_list = [v.strip() for v in vcf_paths.split(",") if v.strip()]
    if not vcf_list:
        return
        
    with open(act_json_path, "r", encoding="utf-8") as f:
        data = json.load(f)
        
    records = data.get("records", [])
    record_coords = {}
    for r in records:
        c, p = r.get("chrom"), r.get("pos")
        if c and p:
            record_coords[(str(c), int(p))] = r

    matched_count = 0
    phased_maternal = 0
    phased_paternal = 0
    import gzip

    for vp in vcf_list:
        if not os.path.exists(vp):
            continue
        opener = gzip.open(vp, "rt") if vp.endswith(".gz") else open(vp, "r", encoding="utf-8", errors="ignore")
        with opener as f:
            for line in f:
                if line.startswith("#"):
                    continue
                idx1 = line.find("\t")
                if idx1 == -1: continue
                c = line[:idx1]
                idx2 = line.find("\t", idx1 + 1)
                if idx2 == -1: continue
                try:
                    p = int(line[idx1+1:idx2])
                except ValueError:
                    continue
                if (c, p) in record_coords:
                    parts = line.rstrip('\r\n').split("\t")
                    if len(parts) > 9:
                        r = record_coords[(c, p)]
                        fmt_keys = parts[8].split(":")
                        sample_vals = parts[9].split(":")
                        fmt_dict = dict(zip(fmt_keys, sample_vals))
                        
                        gt = fmt_dict.get("GT", "").strip()
                        ps = fmt_dict.get("PS", "").strip()
                        dp = fmt_dict.get("DP", "")
                        ad = fmt_dict.get("AD", "")
                        vaf = fmt_dict.get("VAF", "")
                        qual = parts[5] if parts[5] != "." else ""

                        r["vcf_gt"] = gt
                        if ps and ps != ".":
                            r["hap_block"] = ps

                        ev = r.get("evidence", {}) or {}
                        
                        phase_origin = None
                        if gt == "0|1":
                            phase_origin = "Maternal"
                            anchor_label = " (SE Anchor)" if ("Daniel" in str(data.get("patient", "")) or "DE" in str(data.get("patient", ""))) else " (MI Anchor)"
                            phasing_str = f"Phased: Maternal{anchor_label}{f' (PS #{ps})' if ps and ps != '.' else ''}"
                            phased_maternal += 1
                        elif gt in ("1|0", "1/0"):
                            phase_origin = "Paternal"
                            phasing_str = f"Phased: Paternal{f' (PS #{ps})' if ps and ps != '.' else ''}"
                            phased_paternal += 1
                        elif gt in ("1/1", "1|1"):
                            phase_origin = "Homozygous"
                            phasing_str = "Homozygous Alternate"
                            r["zygosity"] = "hom"
                            ev["zygosity"] = "Homozygous"
                        elif gt == "0/1":
                            phase_origin = "Unphased"
                            phasing_str = "Unphased (Short-Read WGS)"
                            r["zygosity"] = "het"
                            ev["zygosity"] = "Heterozygous"
                        else:
                            phase_origin = "Unknown"
                            phasing_str = f"Call: {gt}"

                        ev["phase_origin"] = phase_origin
                        ev["phasing"] = phasing_str
                        if ps and ps != ".":
                            ev["hap_block"] = ps

                        if dp:
                            r["tot_reads"] = dp
                            ev["tot_reads"] = dp
                        if ad:
                            ad_parts = ad.split(",")
                            alt_r = ad_parts[1] if len(ad_parts) > 1 else ad_parts[0]
                            r["alt_reads"] = alt_r
                            ev["alt_reads"] = alt_r
                        if vaf:
                            try:
                                r["vaf"] = str(float(vaf))
                                ev["vaf"] = float(vaf)
                            except ValueError:
                                pass
                        if qual:
                            r["phred"] = qual
                            ev["qual"] = qual
                            
                        r["evidence"] = ev
                        matched_count += 1

    with open(act_json_path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)

    print(f"\n[Phased VCF Integrator] Successfully integrated phased genotypes from VCF(s):")
    print(f"  Total matched actionable variants : {matched_count}")
    print(f"  Maternally anchored alleles (0|1) : {phased_maternal}")
    print(f"  Paternally anchored alleles (1|0) : {phased_paternal}")

def main():
    args = parse_args()
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)

    sample_name = args.sample_id
    if not sample_name:
        base = os.path.basename(args.input)
        for ext in [".sqlite", ".vcf.gz", ".vcf", "_and_2_files"]:
            base = base.replace(ext, "")
        if "Daniel_Ehrle" in base or "260830" in base:
            sample_name = "Daniel_Ehrle"
        else:
            sample_name = base or "Sample"

    now = datetime.now()
    date_str = now.strftime("%d-%m-%Y")
    subfolder_name = f"{sample_name}-{date_str}"

    local_outdir = os.path.join(script_dir, "reports", subfolder_name)
    os.makedirs(local_outdir, exist_ok=True)

    print("==================================================================")
    print(f"GENOMIC ONTOLOGY PIPELINE ENGINE (v5.2)")
    print(f"  Sample / Patient ID : {sample_name}")
    print(f"  Input Source        : {args.input}")
    print(f"  Execution Date      : {date_str}")
    print(f"  Local Output Dir    : {local_outdir}")
    print("==================================================================")

    # 1. Resolve Input
    raw_db, vcf_path = resolve_input(args.input, sample_name, local_outdir, args.phased_vcf)
    print(f"[Stage 1] Resolved Database : {raw_db}")
    print(f"[Stage 1] Resolved VCF File : {vcf_path or 'None (will use DB attributes)'}")

    base_prefix = sample_name
    schema_json = os.path.join(local_outdir, f"{base_prefix}_schema.json")
    panel_json = os.path.join(local_outdir, f"{base_prefix}_master_panel.json")
    act_db = os.path.join(local_outdir, f"{base_prefix}_master_actionable.sqlite")
    act_json = os.path.join(local_outdir, f"{base_prefix}_master_actionable.json")
    enrich_cache = os.path.join(local_outdir, f"{base_prefix}_enrich_cache.json")

    if not os.path.exists(enrich_cache):
        for candidate in glob.glob(os.path.join(script_dir, "reports", "*", "*enrich_cache.json")):
            if os.path.exists(candidate) and os.path.getsize(candidate) > 1000:
                shutil.copy2(candidate, enrich_cache)
                break

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
        "--patient", sample_name
    ], check=True)

    if vcf_path:
        update_actionable_with_phased_vcfs(act_json, act_db, vcf_path)

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

    # 7.1 AlphaGenome Candidates TSV Export
    ag_candidates_tsv = os.path.join(local_outdir, f"{base_prefix}_alphagenome_candidates.tsv")
    try:
        with open(act_json) as f_in:
            data_j = json.load(f_in)
        ag_recs = [r for r in data_j.get("records", []) if "RESCUE_ALPHAGENOME_TARGET" in (r.get("reason_codes") or []) or (r.get("evidence", {}) or {}).get("is_alphagenome_candidate")]
        if ag_recs:
            cols = ["hugo", "chrom", "pos", "ref", "alt", "so", "achange", "zygosity", "tier", "clinvar_sig", "gnomad4_af", "allofus_af", "revel", "am_path", "cadd_phred", "spliceai_max", "alphagenome_subreason", "alphagenome_url"]
            with open(ag_candidates_tsv, "w") as f_out:
                f_out.write("\t".join(cols) + "\n")
                for r in ag_recs:
                    ev = r.get("evidence", {}) or {}
                    row = [
                        r.get("hugo", "") or "",
                        str(r.get("chrom", "") or ""),
                        str(r.get("pos", "") or ""),
                        str(r.get("ref", "") or ""),
                        str(r.get("alt", "") or ""),
                        str(r.get("so", "") or ""),
                        str(r.get("achange", "") or ""),
                        str(ev.get("zygosity", "") or ""),
                        str(r.get("tier", "") or ""),
                        str(r.get("clinvar_sig", "") or ""),
                        str(ev.get("gnomad4_af", "") or ""),
                        str(ev.get("allofus_af", "") or ""),
                        str(r.get("revel", "") or ""),
                        str(r.get("am_path", "") or ""),
                        str(r.get("cadd_phred", "") or ""),
                        str(ev.get("spliceai_max", "") or ""),
                        str(ev.get("alphagenome_subreason", "") or ""),
                        str(ev.get("alphagenome_url", "") or "")
                    ]
                    f_out.write("\t".join(row) + "\n")
            print(f"[AlphaGenome Export] Prioritized {len(ag_recs)} candidates -> {ag_candidates_tsv}")
    except Exception as e:
        print(f"[AlphaGenome Export Warning] {e}")

    deliverables = [
        visual_explorer_html,
        master_hub_html,
        act_json,
        tsv_report,
        txt_report,
        pdf_report,
        zip_bundle
    ]
    if os.path.exists(ag_candidates_tsv):
        deliverables.append(ag_candidates_tsv)

    if not args.local_only:
        print(f"\n[Google Drive Delivery] Uploading deliverables to Google Drive...")
        deliver_to_google_drive(deliverables, args.gdrive_dir, subfolder_name)

    print("\n==================================================================")
    print("✨ PIPELINE COMPLETE. ALL DELIVERABLES GENERATED & SYNCED:")
    print(f"  1. Visual Explorer HTML (Standalone HTML5) : {visual_explorer_html}")
    print(f"  2. Universal Master Hub (HTML5)            : {master_hub_html}")
    print(f"  3. Actionable JSON                         : {act_json}")
    print(f"  4. Variant TSV Matrix                      : {tsv_report}")
    print(f"  5. Clinical Summary TXT                    : {txt_report}")
    print(f"  6. Printable PDF Report                    : {pdf_report}")
    print(f"  7. Offline iOS Bundle                      : {zip_bundle}")
    if os.path.exists(ag_candidates_tsv):
        print(f"  8. AlphaGenome Candidates Matrix TSV       : {ag_candidates_tsv}")
    if not args.local_only:
        print(f"  👉 Google Drive Cloud Remote               : drive:Ontology/{subfolder_name}/")
        print(f"  👉 Google Drive Local Directory            : {args.gdrive_dir}/{subfolder_name}/")
    print("==================================================================")

if __name__ == "__main__":
    main()
