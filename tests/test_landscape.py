#!/usr/bin/env python3
"""
test_landscape.py
Unit tests for the Landscape Clinical Dashboard data processing and report compilation.
"""

import os
import sys
import unittest

SCRIPT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, SCRIPT_DIR)
sys.path.insert(0, os.path.join(SCRIPT_DIR, "scripts"))

import prepare_landscape_data as pld
import build_landscape_report as blr

class TestLandscapeDashboard(unittest.TestCase):

    def setUp(self):
        self.de_json = os.path.join(SCRIPT_DIR, "reports", "Daniel_Ehrle-09-10-2026", "Daniel_Ehrle_ontology_pharma_alphagenome.json")

    def test_01_patient_processing(self):
        """Verify data aggregation produces complete 3-row chart metrics and reports."""
        if not os.path.exists(self.de_json):
            self.skipTest("Sample JSON not present in reports directory")

        pdata = pld.process_patient(self.de_json, "Daniel_Ehrle", "Daniel Ehrle")
        self.assertEqual(pdata["patientId"], "Daniel_Ehrle")
        self.assertIn("charts", pdata)

        # Check Row 1: Ontology
        r1 = pdata["charts"]["row1_ontology"]
        self.assertTrue(len(r1["bar"]) > 0, "Row 1 Bar chart categories missing")
        self.assertTrue(len(r1["piePrior"]) > 0, "Row 1 Pie Prior slices missing")
        self.assertTrue(len(r1["piePost"]) > 0, "Row 1 Pie Post slices missing")

        # Check Row 2: Organ System
        r2 = pdata["charts"]["row2_organ"]
        self.assertTrue(len(r2["bar"]) >= 8, "Row 2 Organ system bars missing")
        self.assertTrue(len(r2["piePrior"]) >= 8, "Row 2 Organ system prior slices missing")

        # Check Row 3: Pharma
        r3 = pdata["charts"]["row3_pharma"]
        self.assertTrue(len(r3["bar"]) >= 5, "Row 3 Pharma classes missing")
        self.assertTrue(len(r3["piePrior"]) >= 3, "Row 3 Metabolizer phenotypes missing")
        self.assertTrue(len(r3["piePost"]) >= 3, "Row 3 Action categories missing")

        # Check Reports
        self.assertIn("reportStats", pdata)
        for r_key in ["visual_explorer", "rare", "clinvar", "cadd", "spliceai", "atlas", "alphamissense", "revel", "ai_deep", "ai_clinical"]:
            self.assertIn(r_key, pdata["reportStats"])

    def test_02_standalone_html_generation(self):
        """Verify build_standalone_landscape_html creates a self-contained file."""
        if not os.path.exists(self.de_json):
            self.skipTest("Sample JSON not present in reports directory")

        test_out = "/tmp/test_landscape_dashboard.html"
        blr.build_standalone_landscape_html("Daniel_Ehrle", self.de_json, test_out)

        self.assertTrue(os.path.exists(test_out), "Output HTML not generated")
        with open(test_out, "r", encoding="utf-8") as f:
            content = f.read()

        # Inlining checks
        self.assertIn("<style>", content)
        self.assertIn("window.COHORT_DATA =", content)
        self.assertIn("chart-row1-bar", content)
        self.assertIn("report-selector", content)
        os.remove(test_out)

if __name__ == "__main__":
    unittest.main()
