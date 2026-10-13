#!/usr/bin/env python3
"""
generate_me_reports.py
Compatibility CLI wrapper that delegates ME report generation to the dynamic,
100% ontology-driven `run_standard_pipeline(target_key="me")` in `run_standard_reporting.py`.
Zero static curations, zero hardcoded gene/variant IDs, and zero PII in source code.
"""

import sys
from run_standard_reporting import run_standard_pipeline


def main():
    return run_standard_pipeline(target_key="me")


if __name__ == "__main__":
    sys.exit(main())
