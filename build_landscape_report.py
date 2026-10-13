#!/usr/bin/env python3
"""
build_landscape_report.py
Builds standalone, 100% self-contained single-file HTML deliverables for the
Landscape Clinical Dashboard with inlined styles, inlined cohort dataset, and logic.
Zero hardcoded PII or static variant identifiers.
"""

import argparse
import json
import os
import shutil
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(SCRIPT_DIR, "scripts"))
import prepare_landscape_data as pld
from lib.zero_curation_guard import enforce_or_exit
enforce_or_exit(SCRIPT_DIR)


def build_standalone_landscape_html(patient_id, json_path, out_html_path):
    if not os.path.exists(json_path):
        raise FileNotFoundError(f"Master JSON not found at {json_path}")

    patient_name = patient_id.replace("_", " ")
    p_data = pld.process_patient(json_path, patient_id, patient_name)

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
        patient_id: p_data
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

    parts = [p for p in str(patient_id or "").split("_") if p]
    code = "".join(p[0] for p in parts).lower() if len(parts) >= 2 else str(patient_id or "")[:2].lower()
    if code == "me":
        html = html.replace('id="cohort-de" class="cohort-btn active"', 'id="cohort-de" class="cohort-btn"')
        html = html.replace('id="cohort-me" class="cohort-btn"', 'id="cohort-me" class="cohort-btn active"')

    os.makedirs(os.path.dirname(os.path.abspath(out_html_path)), exist_ok=True)
    with open(out_html_path, "w", encoding="utf-8") as f:
        f.write(html)

    return out_html_path


def build_standalone_dual_cohort_html(de_info, me_info, out_html_path):
    de_data = pld.process_patient(de_info["json_path"], de_info["patient_id"], de_info["patient_name"])
    me_data = pld.process_patient(me_info["json_path"], me_info["patient_id"], me_info["patient_name"])
    shared_data = pld.build_shared_concordance(de_data, me_data)

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
        de_info["patient_id"]: de_data,
        me_info["patient_id"]: me_data,
        "Shared_Concordance": shared_data
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

    return out_html_path


def main():
    parser = argparse.ArgumentParser(description="Compile standalone Landscape Clinical Dashboard deliverable.")
    parser.add_argument("--patient", default="de", help="Sample code (de, me, or dual)")
    parser.add_argument("--json", default=None, help="Path to input ontology_pharma_alphagenome.json")
    parser.add_argument("--out", default=None, help="Path to output standalone HTML file")
    args = parser.parse_args()

    samples = pld.discover_cohort_samples()
    de_info = samples["de"]
    me_info = samples["me"]

    target = args.patient.lower()
    if target in ["dual", "both", "all"]:
        out_path = args.out or os.path.join(SCRIPT_DIR, "reports", "landscape_clinical_dashboard_standalone.html")
        build_standalone_dual_cohort_html(de_info, me_info, out_path)
        gdrive_ont = os.path.expanduser("~/Google Drive/My Drive/Ontology")
        gdrive_gen = os.path.expanduser("~/Google Drive/My Drive/Genomics Reports")
        if os.path.exists(gdrive_ont):
            shutil.copy2(out_path, os.path.join(gdrive_ont, os.path.basename(out_path)))
        if os.path.exists(os.path.dirname(gdrive_gen)):
            from run_standard_reporting import _sync_to_genomics_reports_with_retention
            _sync_to_genomics_reports_with_retention([out_path], gdrive_gen, max_iterations=5)

        with open(out_path, "r", encoding="utf-8") as f:
            dual_html = f.read()
        dual_html_local = (
            dual_html.replace(f"reports/{de_info['sample_dir']}/", "./")
            .replace(f"reports/{me_info['sample_dir']}/", "./")
        )
        for exp_d in [
            os.path.join(SCRIPT_DIR, "export"),
            os.path.join(SCRIPT_DIR, "reports", "export"),
            os.path.join(gdrive_gen, "export") if os.path.exists(os.path.dirname(gdrive_gen)) else None,
        ]:
            if exp_d:
                os.makedirs(exp_d, exist_ok=True)
                with open(os.path.join(exp_d, os.path.basename(out_path)), "w", encoding="utf-8") as ef:
                    ef.write(dual_html_local)
        return

    info = samples.get(target)
    if not info:
        for k, v in samples.items():
            if v["patient_id"].lower() == target:
                info = v
                break
    if not info:
        info = de_info

    json_path = args.json or info["json_path"]
    out_path = args.out or os.path.join(SCRIPT_DIR, "reports", info["sample_dir"], f"{info['patient_id']}_landscape_dashboard.html")
    build_standalone_landscape_html(info["patient_id"], json_path, out_path)


if __name__ == "__main__":
    main()
