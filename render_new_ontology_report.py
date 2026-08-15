import json
import sys
import os
from typing import Dict, List
from genomics_ontology_io.models import VariantReport, MonogenicFinding, PolygenicRollup, PharmaRecommendation

def generate_html_report(report_data: dict, output_filepath: str):
    # Validate the data using Pydantic models to ensure complete conformance
    report = VariantReport(**report_data)
    
    # Calculate high-level summary metrics
    total_variants = len(report.monogenic_findings)
    tier1_pathogenic = sum(1 for f in report.monogenic_findings if f.clinvar_significance and "Pathogenic" in f.clinvar_significance)
    tier2_vus = sum(1 for f in report.monogenic_findings if f.clinvar_significance and "VUS" in f.clinvar_significance)
    high_polygenic = sum(1 for p in report.polygenic_findings if p.risk_category == "HIGH")
    actionable_pgx = sum(1 for ph in report.pharma_findings if ph.action_tier not in ["STANDARD", "FAVOUR"])

    # Pre-organize monogenic variants by HPO Level 1 Organ System
    organ_systems = {
        "HP:0001626": "Cardiovascular",
        "HP:0002715": "Immune",
        "HP:0000707": "Nervous",
        "HP:0000924": "Skeletal",
        "HP:0001939": "Metabolism",
        "HP:0002664": "Neoplasm / Oncology",
        "HP:0001871": "Blood / Tissues",
        "HP:0003011": "Musculature",
        "HP:0002086": "Respiratory",
        "HP:0000119": "Genitourinary",
        "HP:0000478": "Eye",
        "HP:0000818": "Endocrine",
        "HP:0025031": "Digestive"
    }

    # Generate HTML content
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>Gene Inspector Pro - Genomic & Polygenic Report</title>
    <style>
        :root {{
            --primary-color: #2c3e50;
            --secondary-color: #34495e;
            --accent-color: #3498db;
            --red-bg: #ffd2d2;
            --red-txt: #d8000c;
            --yellow-bg: #fff3cd;
            --yellow-txt: #856404;
            --green-bg: #dff2bf;
            --green-txt: #270000;
            --border-color: #e2e8f0;
        }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
            margin: 0;
            padding: 0;
            color: var(--primary-color);
            background-color: #f7fafc;
        }}
        header {{
            background: linear-gradient(135deg, var(--primary-color) 0%, var(--secondary-color) 100%);
            color: #ffffff;
            padding: 1.5rem 2rem;
            display: flex;
            justify-content: space-between;
            align-items: center;
            box-shadow: 0 4px 6px -1px rgba(0,0,0,0.1);
        }}
        header h1 {{
            margin: 0;
            font-size: 1.8rem;
            font-weight: 700;
            letter-spacing: -0.025em;
        }}
        header .meta {{
            font-size: 0.9rem;
            opacity: 0.9;
            text-align: right;
        }}
        .container {{
            max-width: 1400px;
            margin: 2rem auto;
            padding: 0 1.5rem;
        }}
        
        /* Dashboard Cards */
        .dashboard {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
            gap: 1.5rem;
            margin-bottom: 2rem;
        }}
        .card {{
            background: #ffffff;
            border-radius: 0.5rem;
            padding: 1.5rem;
            box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.05);
            border: 1px solid var(--border-color);
            text-align: center;
        }}
        .card h3 {{
            margin: 0 0 0.5rem 0;
            font-size: 0.875rem;
            text-transform: uppercase;
            letter-spacing: 0.05em;
            color: #718096;
        }}
        .card .value {{
            font-size: 2.25rem;
            font-weight: 800;
            margin: 0;
        }}
        .card.red {{ border-left: 5px solid var(--red-txt); }}
        .card.yellow {{ border-left: 5px solid #d89f00; }}
        .card.green {{ border-left: 5px solid #2d3748; }}
        
        /* Search & Filter Bar */
        .control-panel {{
            background: #ffffff;
            border-radius: 0.5rem;
            padding: 1rem 1.5rem;
            margin-bottom: 2rem;
            box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.05);
            border: 1px solid var(--border-color);
            display: flex;
            gap: 1rem;
            align-items: center;
        }}
        .search-bar {{
            flex: 1;
            padding: 0.75rem 1rem;
            border-radius: 0.375rem;
            border: 1px solid #cbd5e0;
            font-size: 1rem;
        }}
        .search-bar:focus {{
            outline: none;
            border-color: var(--accent-color);
            box-shadow: 0 0 0 3px rgba(66, 153, 225, 0.5);
        }}
        
        /* Section styling */
        .organ-section {{
            background: #ffffff;
            border-radius: 0.5rem;
            padding: 2rem;
            margin-bottom: 2.5rem;
            box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.05);
            border: 1px solid var(--border-color);
        }}
        .organ-header {{
            font-size: 1.5rem;
            border-bottom: 2px solid var(--border-color);
            padding-bottom: 0.5rem;
            margin-bottom: 1.5rem;
            display: flex;
            justify-content: space-between;
            align-items: center;
            color: var(--primary-color);
        }}
        .organ-header h2 {{ margin: 0; font-weight: 800; }}
        .curie-tag {{
            font-size: 0.8rem;
            background: #edf2f7;
            padding: 0.2rem 0.5rem;
            border-radius: 0.25rem;
            color: #4a5568;
            font-family: monospace;
        }}
        
        /* Gene Card and NCBI Header */
        .gene-card {{
            margin-bottom: 2rem;
            padding: 1rem;
            border-radius: 0.375rem;
            background-color: #fafbfd;
            border-left: 3px solid var(--accent-color);
        }}
        .gene-header {{
            font-size: 1.25rem;
            font-weight: 700;
            color: var(--secondary-color);
            margin: 0 0 0.5rem 0;
            display: flex;
            align-items: center;
            gap: 0.5rem;
        }}
        .gene-description {{
            font-size: 0.95rem;
            color: #4a5568;
            line-height: 1.5;
            margin: 0 0 1rem 0;
            font-style: italic;
        }}
        
        /* High-contrast variant tables */
        table.variant-table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 0.5rem;
            background: #ffffff;
        }}
        table.variant-table th, table.variant-table td {{
            padding: 0.75rem 1rem;
            text-align: left;
            border-bottom: 1px solid var(--border-color);
            font-size: 0.9rem;
        }}
        table.variant-table th {{
            background-color: #edf2f7;
            color: var(--secondary-color);
            font-weight: 700;
        }}
        table.variant-table tbody tr:hover {{
            background-color: #f8fafc;
        }}
        
        /* Highlights */
        .badge {{
            display: inline-block;
            padding: 0.25rem 0.5rem;
            font-size: 0.75rem;
            font-weight: 700;
            border-radius: 0.25rem;
            text-transform: uppercase;
        }}
        .badge.pathogenic {{ background-color: var(--red-bg); color: var(--red-txt); }}
        .badge.vus {{ background-color: var(--yellow-bg); color: var(--yellow-txt); }}
        .badge.benign {{ background-color: var(--green-bg); color: var(--green-txt); }}
        
        /* Polygenic Bar Graph */
        .prs-bar-container {{
            width: 100%;
            background-color: #e2e8f0;
            border-radius: 0.25rem;
            height: 1.25rem;
            position: relative;
            overflow: hidden;
            margin: 0.5rem 0;
        }}
        .prs-bar {{
            height: 100%;
            border-radius: 0.25rem;
        }}
        .prs-bar.high {{ background-color: #e53e3e; }}
        .prs-bar.moderate {{ background-color: #dd6b20; }}
        .prs-bar.protective {{ background-color: #319795; }}
        .prs-text {{
            position: absolute;
            top: 50%;
            left: 50%;
            transform: translate(-50%, -50%);
            font-size: 0.75rem;
            font-weight: bold;
            color: #ffffff;
            text-shadow: 0 1px 2px rgba(0,0,0,0.6);
        }}
        
        /* Phasing tag */
        .phase-tag {{
            font-size: 0.8rem;
            font-weight: 600;
            color: #4a5568;
            background: #e2e8f0;
            padding: 0.1rem 0.4rem;
            border-radius: 0.25rem;
        }}
        
        /* Page break rule for printing */
        @media print {{
            .tier-break {{
                page-break-before: always;
            }}
            header, .control-panel {{
                display: none !important;
            }}
        }}
    </style>
</head>
<body>

    <header>
        <div>
            <h1>Gene Inspector Pro - Clinical Report</h1>
            <div style="font-size: 1rem; margin-top: 0.3rem; opacity: 0.85;">Ontology & Polygenic Risk Integration Engine</div>
        </div>
        <div class="meta">
            Patient ID: <strong>{report.patient_id}</strong><br>
            Analysis Date: {report.run_date}<br>
            Validated by <strong>LinkML v1.4 / Pydantic v2</strong>
        </div>
    </header>

    <div class="container">
        
        <!-- Dashboard Summary -->
        <div class="dashboard">
            <div class="card green">
                <h3>Total Variants</h3>
                <div class="value">{total_variants}</div>
            </div>
            <div class="card red">
                <h3>Actionable Mendelian (Tier 1)</h3>
                <div class="value">{tier1_pathogenic}</div>
            </div>
            <div class="card yellow">
                <h3>VUS of Interest (Tier 2)</h3>
                <div class="value">{tier2_vus}</div>
            </div>
            <div class="card">
                <h3>High Polygenic Risks</h3>
                <div class="value">{high_polygenic}</div>
            </div>
            <div class="card">
                <h3>Actionable PGx Cards</h3>
                <div class="value">{actionable_pgx}</div>
            </div>
        </div>

        <!-- Search Bar -->
        <div class="control-panel">
            <span style="font-weight: bold; color: var(--secondary-color);">Dynamic Filter:</span>
            <input type="text" class="search-bar" id="searchBar" placeholder="Search by Gene, Variant, dbSNP rsID, HPO Phenotype, or Clinical Category..." onkeyup="filterReport()">
        </div>

        <!-- Section: Monogenic Clinical Findings (ACMG 3-Tier Classification) -->
        <h2 style="font-size: 1.75rem; margin-top: 3rem; margin-bottom: 1.5rem; font-weight: 800; border-bottom: 3px solid var(--primary-color); padding-bottom: 0.5rem;">
            🧬 Monogenic Findings (ACMG 3-Tier Classification)
        </h2>
"""
    
    # Process findings by Organ System HPO tags
    for curie, system_name in organ_systems.items():
        system_monogenic = [f for f in report.monogenic_findings if curie in f.associated_hpo_terms]
        system_polygenic = [p for p in report.polygenic_findings if p.hpo_level1_system == curie]
        
        if not system_monogenic and not system_polygenic:
            continue
            
        html += f"""
        <div class="organ-section target-section">
            <div class="organ-header">
                <h2>{system_name} Phenotypes</h2>
                <span class="curie-tag">{curie}</span>
            </div>
        """

        # 1. Display Polygenic Risk Roll-ups first under the organ node if they exist
        if system_polygenic:
            html += f"""
            <h3 style="font-size: 1.15rem; color: var(--secondary-color); margin-top: 0; margin-bottom: 1rem; border-bottom: 1px dashed var(--border-color); padding-bottom: 0.3rem;">
                📊 Polygenic Risk Roll-ups (Trait-Level Percentiles)
            </h3>
            <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 1rem; margin-bottom: 2rem;">
            """
            for prs in system_polygenic:
                bar_color_class = prs.risk_category.lower()
                html += f"""
                <div style="background: #f8fafc; padding: 1rem; border-radius: 0.375rem; border: 1px solid var(--border-color);">
                    <div style="display: flex; justify-content: space-between; font-weight: bold; font-size: 0.95rem; margin-bottom: 0.5rem;">
                        <span>{prs.trait_name}</span>
                        <span class="curie-tag">{prs.efo_trait_id}</span>
                    </div>
                    <div class="prs-bar-container">
                        <div class="prs-bar {bar_color_class}" style="width: {prs.percentile}%;"></div>
                        <div class="prs-text">{prs.risk_category} RISK ({prs.percentile:.1f}th Percentile)</div>
                    </div>
                    <div style="font-size: 0.8rem; color: #718096; text-align: right; margin-top: 0.25rem;">
                        PGS Catalog Source: {prs.pgs_catalog_id or 'N/A'} (Level 2: {prs.hpo_level2_subcategory})
                    </div>
                </div>
                """
            html += "</div>"

        # 2. Display Monogenic Variants grouped by Gene Symbol
        if system_monogenic:
            html += f"""
            <h3 style="font-size: 1.15rem; color: var(--secondary-color); margin-top: 1rem; margin-bottom: 1rem; border-bottom: 1px dashed var(--border-color); padding-bottom: 0.3rem;">
                Mendelian High-Penetrance Variants
            </h3>
            """
            # Group variants by gene
            variants_by_gene: Dict[str, List[MonogenicFinding]] = {}
            for f in system_monogenic:
                variants_by_gene.setdefault(f.gene_symbol, []).append(f)
                
            for gene_sym, findings in variants_by_gene.items():
                ncbi_desc = findings[0].ncbi_description or "No NCBI description available for this gene."
                html += f"""
                <div class="gene-card">
                    <h4 class="gene-header">🧬 {gene_sym}</h4>
                    <p class="gene-description">{ncbi_desc}</p>
                    <table class="variant-table">
                        <thead>
                            <tr>
                                <th>Variant / RSID</th>
                                <th>Consequence</th>
                                <th>Zygosity</th>
                                <th>Phase</th>
                                <th>ClinVar Records</th>
                                <th>REVEL</th>
                                <th>Coordinates</th>
                            </tr>
                        </thead>
                        <tbody>
                """
                for f in findings:
                    clinvar_badge = "vus"
                    if f.clinvar_significance and "Pathogenic" in f.clinvar_significance:
                        clinvar_badge = "pathogenic"
                    elif f.clinvar_significance and "Benign" in f.clinvar_significance:
                        clinvar_badge = "benign"
                        
                    html += f"""
                            <tr class="variant-row">
                                <td><strong>{f.rsid or 'Novel Variant'}</strong></td>
                                <td>{f.impact_consequence}</td>
                                <td>{f.zygosity}</td>
                                <td><span class="phase-tag">{f.phasing}</span></td>
                                <td><span class="badge {clinvar_badge}">{f.clinvar_significance or 'VUS'}</span></td>
                                <td style="font-weight: bold; color: {'#e53e3e' if (f.revel_score or 0) > 0.75 else '#4a5568'};">{f.revel_score or 'N/A'}</td>
                                <td style="font-size: 0.8rem; color: #718096; font-family: monospace;">{f.chromosome}:{f.position}</td>
                            </tr>
                    """
                html += """
                        </tbody>
                    </table>
                </div>
                """
        html += "</div>"

    # Section: Pharmacogenomics (PGx) Integrated Profile
    if report.pharma_findings:
        html += f"""
        <div class="tier-break"></div>
        <h2 style="font-size: 1.75rem; margin-top: 3rem; margin-bottom: 1.5rem; font-weight: 800; border-bottom: 3px solid var(--primary-color); padding-bottom: 0.5rem;">
            💊 Pharmacogenomics (PGx) Profile
        </h2>
        <div class="organ-section">
            <table class="variant-table" style="width:100%;">
                <thead>
                    <tr>
                        <th>Pharmacogene</th>
                        <th>Genotype (Star Alleles)</th>
                        <th>Predicted Phenotype</th>
                        <th>Drug Class / Affected Drug</th>
                        <th>prescribing Recommendations (CPIC / DPWG guidelines)</th>
                        <th>Action Actionability</th>
                    </tr>
                </thead>
                <tbody>
        """
        for ph in report.pharma_findings:
            pharma_badge = "vus" if ph.action_tier in ["CAUTION", "MONITOR"] else ("pathogenic" if ph.action_tier in ["AVOID", "DOSE_DOWN"] else "benign")
            html += f"""
                    <tr class="variant-row">
                        <td><strong>{ph.gene}</strong></td>
                        <td style="font-family: monospace; font-size: 1rem;">{ph.diplotype}</td>
                        <td>{ph.phenotype or 'N/A'}</td>
                        <td><strong>{ph.affected_drug}</strong></td>
                        <td style="font-size: 0.85rem; line-height: 1.4; color: #4a5568;">{ph.clinical_recommendation}</td>
                        <td><span class="badge {pharma_badge}">{ph.action_tier}</span></td>
                    </tr>
            """
        html += """
                </tbody>
            </table>
        </div>
        """

    # Add interactive search logic in vanilla JS
    html += """
    </div>

    <script>
        function filterReport() {
            var input = document.getElementById("searchBar");
            var filter = input.value.toLowerCase();
            
            // Filter target section divs (organ system cards)
            var sections = document.getElementsByClassName("target-section");
            
            for (var i = 0; i < sections.length; i++) {
                var section = sections[i];
                var geneCards = section.getElementsByClassName("gene-card");
                var sectionHasMatch = false;
                
                // Search inside each gene card
                for (var j = 0; j < geneCards.length; j++) {
                    var card = geneCards[j];
                    var textContent = card.textContent || card.innerText;
                    
                    if (textContent.toLowerCase().indexOf(filter) > -1) {
                        card.style.display = "";
                        sectionHasMatch = true;
                    } else {
                        card.style.display = "none";
                    }
                }
                
                // If there are polygenic risk roll-ups or if a gene card matched
                var prsContainer = section.querySelector("div[style*='grid-template-columns']");
                var prsMatch = false;
                if (prsContainer) {
                    var prsBlocks = prsContainer.children;
                    for (var k = 0; k < prsBlocks.length; k++) {
                        var block = prsBlocks[k];
                        if (block.textContent.toLowerCase().indexOf(filter) > -1) {
                            block.style.display = "";
                            prsMatch = true;
                            sectionHasMatch = true;
                        } else {
                            block.style.display = "none";
                        }
                    }
                }
                
                if (sectionHasMatch) {
                    section.style.display = "";
                } else {
                    section.style.display = "none";
                }
            }
        }
    </script>
</body>
</html>
"""

    with open(output_filepath, "w") as f:
        f.write(html)
    print(f"Successfully generated HTML report at: {output_filepath}")
