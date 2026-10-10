#!/usr/bin/env python3
"""
build_landscape_report.py
Builds standalone, 100% self-contained single-file HTML deliverables for the
Landscape Clinical Dashboard with inlined styles, inlined cohort dataset, and logic.
"""

import argparse
import json
import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(SCRIPT_DIR, "scripts"))
import prepare_landscape_data as pld

def build_standalone_landscape_html(patient_id, json_path, out_html_path):
    print("=" * 70)
    print(f"BUILDING STANDALONE LANDSCAPE DASHBOARD: {patient_id}")
    print(f"  Input Master JSON : {json_path}")
    print(f"  Target Deliverable: {out_html_path}")
    print("=" * 70)

    if not os.path.exists(json_path):
        raise FileNotFoundError(f"Master JSON not found at {json_path}")

    # Process patient dataset
    patient_name = patient_id.replace("_", " ")
    p_data = pld.process_patient(json_path, patient_id, patient_name)

    # Read template, CSS, and JS
    template_path = os.path.join(SCRIPT_DIR, "landscape_clinical_dashboard.html")
    css_path = os.path.join(SCRIPT_DIR, "css", "landscape.css")
    js_path = os.path.join(SCRIPT_DIR, "js", "landscape.js")

    with open(template_path, "r", encoding="utf-8") as f:
        html = f.read()
    with open(css_path, "r", encoding="utf-8") as f:
        css = f.read()
    with open(js_path, "r", encoding="utf-8") as f:
        js = f.read()

    # Create embedded cohort payload with active patient
    cohort_payload = {
        patient_id: p_data
    }
    data_js_inline = f"window.COHORT_DATA = {json.dumps(cohort_payload)};\n"

    # Replace external links with inlined content
    html = html.replace(
        '<link rel="stylesheet" href="css/landscape.css" />',
        f"<style>\n{css}\n</style>"
    )
    html = html.replace(
        '<script src="data/cohort_landscape_data.js"></script>',
        f"<script>\n{data_js_inline}\n</script>"
    )
    html = html.replace(
        '<script src="js/landscape.js"></script>',
        f"<script>\n{js}\n</script>"
    )

    # Adjust cohort buttons in header for standalone single-patient view
    if patient_id == "Melinda_Ehrle":
        html = html.replace('id="cohort-de" class="cohort-btn active"', 'id="cohort-de" class="cohort-btn"')
        html = html.replace('id="cohort-me" class="cohort-btn"', 'id="cohort-me" class="cohort-btn active"')

    os.makedirs(os.path.dirname(os.path.abspath(out_html_path)), exist_ok=True)
    with open(out_html_path, "w", encoding="utf-8") as f:
        f.write(html)

    size_mb = os.path.getsize(out_html_path) / (1024 * 1024)
    print(f"Successfully compiled standalone Landscape Dashboard ({size_mb:.2f} MB): {out_html_path}")
    return out_html_path

def build_standalone_dual_cohort_html(de_json, me_json, out_html_path):
    print("=" * 70)
    print(f"BUILDING STANDALONE DUAL-COHORT LANDSCAPE PORTAL")
    print(f"  Target Deliverable: {out_html_path}")
    print("=" * 70)

    de_data = pld.process_patient(de_json, "Daniel_Ehrle", "Daniel Ehrle")
    me_data = pld.process_patient(me_json, "Melinda_Ehrle", "Melinda Ehrle")

    template_path = os.path.join(SCRIPT_DIR, "landscape_clinical_dashboard.html")
    css_path = os.path.join(SCRIPT_DIR, "css", "landscape.css")
    js_path = os.path.join(SCRIPT_DIR, "js", "landscape.js")

    with open(template_path, "r", encoding="utf-8") as f:
        html = f.read()
    with open(css_path, "r", encoding="utf-8") as f:
        css = f.read()
    with open(js_path, "r", encoding="utf-8") as f:
        js = f.read()

    cohort_payload = {
        "Daniel_Ehrle": de_data,
        "Melinda_Ehrle": me_data
    }
    data_js_inline = f"window.COHORT_DATA = {json.dumps(cohort_payload)};\n"

    html = html.replace(
        '<link rel="stylesheet" href="css/landscape.css" />',
        f"<style>\n{css}\n</style>"
    )
    html = html.replace(
        '<script src="data/cohort_landscape_data.js"></script>',
        f"<script>\n{data_js_inline}\n</script>"
    )
    html = html.replace(
        '<script src="js/landscape.js"></script>',
        f"<script>\n{js}\n</script>"
    )

    os.makedirs(os.path.dirname(os.path.abspath(out_html_path)), exist_ok=True)
    with open(out_html_path, "w", encoding="utf-8") as f:
        f.write(html)

    size_mb = os.path.getsize(out_html_path) / (1024 * 1024)
    print(f"Successfully compiled dual-cohort standalone portal ({size_mb:.2f} MB): {out_html_path}")
    return out_html_path

def main():
    parser = argparse.ArgumentParser(description="Compile standalone Landscape Clinical Dashboard deliverable.")
    parser.add_argument("--patient", default="Daniel_Ehrle", help="Patient ID (Daniel_Ehrle, Melinda_Ehrle, or dual)")
    parser.add_argument("--json", default=None, help="Path to input ontology_pharma_alphagenome.json")
    parser.add_argument("--out", default=None, help="Path to output standalone HTML file")
    args = parser.parse_args()

    patient_id = args.patient
    de_json = os.path.join(SCRIPT_DIR, "reports", "Daniel_Ehrle-09-10-2026", "Daniel_Ehrle_ontology_pharma_alphagenome.json")
    me_json = os.path.join(SCRIPT_DIR, "reports", "Melinda_Ehrle-09-10-2026", "Melinda_Ehrle_ontology_pharma_alphagenome.json")

    if patient_id.lower() in ["dual", "both", "all"]:
        out_path = args.out or os.path.join(SCRIPT_DIR, "reports", "landscape_clinical_dashboard_standalone.html")
        build_standalone_dual_cohort_html(de_json, me_json, out_path)
        return

    if not args.json:
        date_folder = "Daniel_Ehrle-09-10-2026" if "Daniel" in patient_id else "Melinda_Ehrle-09-10-2026"
        args.json = os.path.join(SCRIPT_DIR, "reports", date_folder, f"{patient_id}_ontology_pharma_alphagenome.json")

    if not args.out:
        date_folder = "Daniel_Ehrle-09-10-2026" if "Daniel" in patient_id else "Melinda_Ehrle-09-10-2026"
        args.out = os.path.join(SCRIPT_DIR, "reports", date_folder, f"{patient_id}_landscape_dashboard.html")

    build_standalone_landscape_html(patient_id, args.json, args.out)

if __name__ == "__main__":
    main()
