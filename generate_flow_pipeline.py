#!/usr/bin/env python3
"""
generate_flow_pipeline.py
Compatibility CLI for the End-to-End Clinical Genomics Report Pipeline (v5.3).
Delegates to run_standard_pipeline in run_standard_reporting.py for multi-cohort execution.

Flow Architecture:
  Stage 1: Local ingestion of master JSON (enriched with AlphaGenome & Pharmacogenomics).
  Stage 2: Automatic local PII stripping generating sanitized zero-PII artifacts.
  Stage 3: Evidence structuring across 6 clinical domains according to deep research skills.
  Stage 4: Synthesis of standard deliverables:
           - 3-Page Clinician Action Brief (Alt C format)
           - Deep Research Dossier (Alt B format)
           - Patient Self-Reported EHR Import JSON (FHIR bundle format)
  Stage 5: Local model memory unload & RAM reclamation.
  Stage 6: Automated sync to Google Drive dated run folders and main Ontology root.
"""

import os
import sys
import argparse
from run_standard_reporting import run_standard_pipeline

def parse_args():
    parser = argparse.ArgumentParser(description="End-to-End Clinical Genomics Report Pipeline")
    parser.add_argument("--sample-id", default="Daniel_Ehrle-07-10-2026", help="Dated sample folder identifier")
    parser.add_argument("--patient-name", default="Daniel Ehrle", help="Real patient name for local binding")
    parser.add_argument("--patient-id", default="Daniel_Ehrle", help="Real patient sample ID")
    parser.add_argument("--reports-dir", default="/home/daniel-ehrle/My-Projects/genomics/ontology_report/reports", help="Base reports directory")
    parser.add_argument("--gdrive-dir", default="/home/daniel-ehrle/Google Drive/My Drive/Ontology", help="Google Drive destination directory")
    parser.add_argument("--dry-run", action="store_true", help="Validate inputs without writing files")
    return parser.parse_args()

def main():
    args = parse_args()
    target_key = "me" if ("melinda" in args.patient_id.lower() or "melinda" in args.patient_name.lower()) else "de"
    sys.exit(run_standard_pipeline(
        target_key=target_key,
        reports_dir=args.reports_dir,
        gdrive_dir=args.gdrive_dir,
        dry_run=args.dry_run,
        sample_dir_override=args.sample_id,
        patient_name_override=args.patient_name,
        patient_id_override=args.patient_id
    ))

if __name__ == "__main__":
    main()
