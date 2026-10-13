#!/usr/bin/env python3
"""
generate_flow_pipeline.py
Compatibility CLI for the End-to-End Clinical Genomics Report Pipeline.
Delegates to run_standard_pipeline in run_standard_reporting.py for multi-cohort execution.
Zero hardcoded PII.
"""

import argparse
import os
import sys
from run_standard_reporting import run_standard_pipeline

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


def parse_args():
    parser = argparse.ArgumentParser(description="End-to-End Clinical Genomics Report Pipeline")
    parser.add_argument("--target", default="de", help="Sample code (de or me)")
    parser.add_argument("--sample-id", default=None, help="Dated sample folder identifier")
    parser.add_argument("--patient-name", default=None, help="Patient name for local binding")
    parser.add_argument("--patient-id", default=None, help="Patient sample ID")
    parser.add_argument("--reports-dir", default=os.path.join(SCRIPT_DIR, "reports"), help="Base reports directory")
    parser.add_argument("--gdrive-dir", default=os.path.expanduser("~/Google Drive/My Drive/Ontology"), help="Google Drive destination directory")
    parser.add_argument("--dry-run", action="store_true", help="Validate inputs without writing files")
    return parser.parse_args()


def main():
    args = parse_args()
    sys.exit(run_standard_pipeline(
        target_key=args.target,
        reports_dir=args.reports_dir,
        gdrive_dir=args.gdrive_dir,
        dry_run=args.dry_run,
        sample_dir_override=args.sample_id,
        patient_name_override=args.patient_name,
        patient_id_override=args.patient_id
    ))


if __name__ == "__main__":
    main()
