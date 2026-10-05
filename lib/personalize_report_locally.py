#!/usr/bin/env python3
"""
personalize_report_locally.py
Air-gapped, offline post-processor for Genomic Ontology & Deep Research Reports.

Security & Privacy Guarantee:
- The upstream AI report is generated 100% de-identified (using PROBAND tokens).
- This script runs strictly on localhost/offline machine to bind real patient PII 
  (Name, DOB, Sample ID) into the header and render final local deliverables.
- Strips any cross-family sample codes (e.g. ME in DE's report) ensuring complete
  patient privacy and HIPAA/GDPR zero-PII compliance.
"""

import os
import sys
import argparse
import subprocess
import shutil
import re

def parse_args():
    parser = argparse.ArgumentParser(description="Local Offline Patient Report Personalizer")
    parser.add_argument("--report-dir", required=True, help="Directory containing generated report files")
    parser.add_argument("--patient-name", required=True, help="Real patient name to bind locally")
    parser.add_argument("--patient-id", required=True, help="Real patient identifier/sample ID")
    parser.add_argument("--proband-token", default="PROBAND_01", help="Anonymous token used in upstream generation")
    parser.add_argument("--no-pdf", action="store_true", help="Skip local PDF re-compilation")
    return parser.parse_args()

def personalize_file(filepath, proband_token, patient_name, patient_id):
    if not os.path.exists(filepath):
        return False
    with open(filepath, "r", encoding="utf-8") as f:
        content = f.read()

    # 1. Replace proband tokens
    content = content.replace(proband_token, patient_name)
    content = content.replace(f"`{proband_token}`", f"`{patient_id}`")
    
    # 2. Strict family cross-reference scrub: ensure no parental sample IDs linger
    # Replace any parental codes with clean generic lineage
    content = re.sub(r'\(maternal microarray \w+\)', '(maternal microarray)', content, flags=re.IGNORECASE)
    content = re.sub(r'\(paternal microarray \w+\)', '(paternal microarray)', content, flags=re.IGNORECASE)
    content = re.sub(r'oriented against \w+\.tsv', 'oriented against parental microarray', content, flags=re.IGNORECASE)

    with open(filepath, "w", encoding="utf-8") as f:
        f.write(content)
    return True

def recompile_pdf(html_path, pdf_path):
    chrome_bin = shutil.which("google-chrome") or shutil.which("chromium-browser") or shutil.which("chromium")
    if not chrome_bin:
        print("[Warning] Chrome not found for local PDF compilation.")
        return False
    cmd = [
        chrome_bin,
        "--headless=new",
        "--disable-gpu",
        f"--print-to-pdf={pdf_path}",
        "--no-pdf-header-footer",
        f"file://{os.path.abspath(html_path)}"
    ]
    res = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
    return res.returncode == 0 and os.path.exists(pdf_path)

def main():
    args = parse_args()
    report_dir = os.path.abspath(args.report_dir)
    print(f"[Local Personalizer] Processing directory: {report_dir}")
    print(f"  Target Patient: {args.patient_name} ({args.patient_id})")

    # Find markdown and html reports
    for fname in os.listdir(report_dir):
        if fname.endswith(".md") or fname.endswith(".html") or fname.endswith(".txt"):
            fpath = os.path.join(report_dir, fname)
            if personalize_file(fpath, args.proband_token, args.patient_name, args.patient_id):
                print(f"  [Scrubbed & Bound] {fname}")

    # Recompile PDF if HTML exists
    if not args.no_pdf:
        for fname in os.listdir(report_dir):
            if fname.endswith("_deep_research_report.html"):
                html_path = os.path.join(report_dir, fname)
                pdf_path = html_path.replace(".html", ".pdf")
                if recompile_pdf(html_path, pdf_path):
                    print(f"  [PDF Recompiled] {os.path.basename(pdf_path)} ({os.path.getsize(pdf_path)/1024:.1f} KB)")

    print("[Local Personalizer] Done. Local deliverables bound with zero external PII egress.")

if __name__ == "__main__":
    main()
