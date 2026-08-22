#!/usr/bin/env python3
"""
Visual Ontology Explorer - High-Fidelity Clinical Report Generator
Generates a clean, modern, interactive 2-pane clinical genomics report interface.
Features:
- Branded Clinical 2-Pane Architecture (Left Index Panel + Right Detailed Card Stream)
- Live Instant Search & Panel Filtering (REVEL, ClinVar, Phased, Flagged)
- Built-in Clinical HPO & Phenotype Dictionary Resolution
- Discrete Allele Genotype Badges, MAF Filters, and Maternal/Paternal Phasing
"""

import json
import os
import sys
import argparse
from pathlib import Path
from typing import Dict, List, Any

# Clinical HPO Phenotype Dictionary for human-readable resolution
HPO_DICTIONARY = {
    "HP:0001626": "Cardiovascular system abnormality",
    "HP:0001644": "Dilated cardiomyopathy",
    "HP:0001639": "Hypertrophic cardiomyopathy",
    "HP:0001635": "Cardiac arrhythmia / Conduction disease",
    "HP:0001678": "Atrioventricular block",
    "HP:0004756": "Long QT syndrome",
    "HP:0001662": "Brugada syndrome",
    "HP:0002664": "Neoplasm / Tumorigenesis",
    "HP:0003002": "Breast carcinoma susceptibility",
    "HP:0001250": "Seizure / Neurodevelopmental abnormality",
    "HP:0000707": "Abnormality of the nervous system",
    "HP:0000365": "Hearing impairment / Deafness",
    "HP:0001000": "Abnormality of skin morphology / Keratosis",
    "HP:0001427": "Mitochondrial respiratory chain deficiency",
    "HP:0000100": "Nephrotic syndrome / Renal abnormality",
    "HP:0003473": "Hyperkalemic periodic paralysis",
    "HP:0003470": "Paramyotonia congenita",
    "HP:0001249": "Intellectual disability",
    "HP:0001297": "Stroke-like episode",
    "HP:0000822": "Hypertension",
    "HP:0002099": "Asthma / Respiratory allergy",
    "HP:0002715": "Abnormality of the immune system",
    "HP:0001903": "Anemia",
    "HP:0001873": "Thrombocytopenia"
}

def resolve_hpo_term(hpo_id: str) -> str:
    """Resolves an HPO ID to a clinical phenotype name with fallback."""
    if not hpo_id:
        return ""
    clean_id = hpo_id.strip()
    if clean_id in HPO_DICTIONARY:
        return HPO_DICTIONARY[clean_id]
    return clean_id

def generate_report_html(report_data: dict, output_filepath: str) -> str:
    patient_id = report_data.get("patient_id", "PATIENT_WGS")
    run_date = report_data.get("run_date", "2026-08-21 18:30 UTC")
    monogenic = report_data.get("monogenic_findings", [])
    polygenic = report_data.get("polygenic_findings", [])
    pharma = report_data.get("pharma_findings", [])

    # Group findings by gene
    genes_dict = {}
    for item in monogenic:
        sym = item.get("gene_symbol", "Unknown")
        if sym not in genes_dict:
            # Resolve HPO terms
            hpo_list = item.get("associated_hpo_terms", [])
            resolved_hpos = []
            for h in hpo_list:
                term_name = resolve_hpo_term(h)
                resolved_hpos.append({"id": h, "name": term_name})

            genes_dict[sym] = {
                "symbol": sym,
                "name": item.get("gene_name", ""),
                "description": item.get("ncbi_description", ""),
                "omim_source": item.get("omim_source", "OMIM"),
                "pathologies": item.get("pathologies", []),
                "resolved_hpos": resolved_hpos,
                "variants": [],
                "max_revel": 0.0,
                "has_pathogenic": False,
                "has_phased": False
            }

        # Calculate max revel
        rev = item.get("revel_score")
        if rev is not None:
            try:
                rev_f = float(rev)
                if rev_f > genes_dict[sym]["max_revel"]:
                    genes_dict[sym]["max_revel"] = rev_f
            except (ValueError, TypeError):
                pass

        # Check clinvar
        sig = str(item.get("clinvar_significance", "")).lower()
        if "pathogenic" in sig:
            genes_dict[sym]["has_pathogenic"] = True

        # Check phasing
        ph = str(item.get("phasing", "")).lower()
        if ph in ["maternal", "paternal", "de_novo", "compound_het"]:
            genes_dict[sym]["has_phased"] = True

        genes_dict[sym]["variants"].append(item)

    # Sort genes by max REVEL score descending
    sorted_genes = sorted(genes_dict.values(), key=lambda g: g["max_revel"], reverse=True)
    total_variants_count = len(monogenic)

    # Serialize data for client-side JavaScript
    embedded_data_json = json.dumps({
        "patient_id": patient_id,
        "run_date": run_date,
        "genes": sorted_genes,
        "total_variants": total_variants_count,
        "polygenic": polygenic,
        "pharma": pharma
    }, indent=2)

    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Visual Ontology Explorer - Genomics Report</title>
    <!-- Modern typography & icons -->
    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700;800&family=JetBrains+Mono:wght@400;500;600&display=swap" rel="stylesheet">
    <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0/css/all.min.css">
    
    <style>
        :root {{
            --bg-page: #f8fafc;
            --bg-surface: #ffffff;
            --bg-sidebar: #ffffff;
            --border-color: #e2e8f0;
            --border-light: #edf2f7;
            --text-main: #0f172a;
            --text-secondary: #475569;
            --text-muted: #64748b;
            --text-light: #94a3b8;
            --primary: #2563eb;
            --primary-hover: #1d4ed8;
            --primary-light: #eff6ff;
            --danger: #dc2626;
            --danger-bg: #fee2e2;
            --warning: #d97706;
            --warning-bg: #fef3c7;
            --success: #16a34a;
            --success-bg: #dcfce7;
            --purple: #7c3aed;
            --purple-bg: #f3e8ff;
            --slate-badge: #f1f5f9;
        }}

        * {{
            box-sizing: border-box;
            margin: 0;
            padding: 0;
        }}

        body {{
            font-family: 'Inter', system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
            background-color: var(--bg-page);
            color: var(--text-main);
            height: 100vh;
            display: flex;
            flex-direction: column;
            overflow: hidden;
            -webkit-font-smoothing: antialiased;
        }}

        /* Header Navigation */
        .top-navbar {{
            background: #ffffff;
            border-bottom: 1px solid var(--border-color);
            padding: 10px 24px;
            display: flex;
            justify-content: space-between;
            align-items: center;
            height: 64px;
            shrink: 0;
            z-index: 30;
        }}

        .brand-section {{
            display: flex;
            align-items: center;
            gap: 12px;
        }}

        .brand-logo {{
            width: 38px;
            height: 38px;
            background: linear-gradient(135deg, #2563eb, #3b82f6);
            color: white;
            border-radius: 10px;
            display: flex;
            align-items: center;
            justify-content: center;
            font-size: 20px;
            box-shadow: 0 2px 6px rgba(37, 99, 235, 0.25);
        }}

        .brand-title {{
            font-size: 17px;
            font-weight: 700;
            color: var(--text-main);
            letter-spacing: -0.01em;
            display: flex;
            align-items: center;
            gap: 8px;
        }}

        .brand-subtitle {{
            font-size: 11px;
            color: var(--text-muted);
            font-weight: 500;
        }}

        /* Search Bar in Header */
        .search-container {{
            position: relative;
            width: 440px;
        }}

        .search-icon {{
            position: absolute;
            left: 14px;
            top: 50%;
            transform: translateY(-50%);
            color: #94a3b8;
            font-size: 14px;
        }}

        .search-input {{
            width: 100%;
            height: 40px;
            background: #f1f5f9;
            border: 1px solid #e2e8f0;
            border-radius: 20px;
            padding: 8px 16px 8px 38px;
            font-size: 13px;
            color: var(--text-main);
            outline: none;
            transition: all 0.2s ease;
        }}

        .search-input:focus {{
            background: #ffffff;
            border-color: var(--primary);
            box-shadow: 0 0 0 3px rgba(37, 99, 235, 0.15);
        }}

        .user-meta {{
            display: flex;
            align-items: center;
            gap: 16px;
        }}

        .patient-pill {{
            background: #f8fafc;
            border: 1px solid #e2e8f0;
            padding: 6px 14px;
            border-radius: 8px;
            font-size: 12px;
            text-align: right;
            line-height: 1.3;
        }}

        .avatar-circle {{
            width: 36px;
            height: 36px;
            border-radius: 50%;
            background: #0f172a;
            color: #ffffff;
            display: flex;
            align-items: center;
            justify-content: center;
            font-weight: 700;
            font-size: 13px;
        }}

        /* Secondary Panel Filter Tabs Bar */
        .panel-nav-bar {{
            background: #ffffff;
            border-bottom: 1px solid var(--border-color);
            padding: 6px 24px;
            display: flex;
            gap: 8px;
            align-items: center;
        }}

        .nav-pill {{
            padding: 6px 14px;
            border-radius: 6px;
            font-size: 12.5px;
            font-weight: 600;
            color: var(--text-muted);
            background: transparent;
            border: none;
            cursor: pointer;
            transition: all 0.15s ease;
        }}

        .nav-pill:hover {{
            color: var(--primary);
            background: #f1f5f9;
        }}

        .nav-pill.active {{
            color: var(--primary);
            background: var(--primary-light);
            font-weight: 700;
        }}

        /* Main Workspace Container (Two Pane) */
        .workspace {{
            display: flex;
            flex: 1;
            overflow: hidden;
            background: var(--bg-page);
        }}

        /* Left Sidebar: Gene Index */
        .sidebar-panel {{
            width: 360px;
            background: var(--bg-sidebar);
            border-right: 1px solid var(--border-color);
            display: flex;
            flex-direction: column;
            shrink: 0;
            height: 100%;
        }}

        .sidebar-header {{
            padding: 16px 18px;
            border-bottom: 1px solid var(--border-color);
        }}

        .sidebar-title-row {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 10px;
        }}

        .sidebar-title {{
            font-size: 14px;
            font-weight: 700;
            color: var(--text-main);
            display: flex;
            align-items: center;
            gap: 6px;
        }}

        .btn-report-badge {{
            background: #f1f5f9;
            border: 1px solid #cbd5e1;
            color: #1e293b;
            font-size: 11px;
            font-weight: 600;
            padding: 4px 10px;
            border-radius: 4px;
            cursor: pointer;
        }}

        .sidebar-search-box {{
            width: 100%;
            height: 34px;
            background: #ffffff;
            border: 1px solid #cbd5e1;
            border-radius: 6px;
            padding: 4px 10px;
            font-size: 12px;
            outline: none;
        }}

        .sidebar-search-box:focus {{
            border-color: var(--primary);
        }}

        .gene-list {{
            flex: 1;
            overflow-y: auto;
            list-style: none;
        }}

        .gene-item {{
            padding: 12px 18px;
            border-bottom: 1px solid var(--border-light);
            cursor: pointer;
            display: flex;
            justify-content: space-between;
            align-items: center;
            transition: all 0.15s ease;
            position: relative;
        }}

        .gene-item:hover {{
            background: #f8fafc;
        }}

        .gene-item.active {{
            background: #eff6ff;
            border-left: 3px solid var(--primary);
        }}

        .gene-item-info {{
            display: flex;
            align-items: flex-start;
            gap: 10px;
            flex: 1;
            overflow: hidden;
        }}

        .priority-dot {{
            width: 6px;
            height: 6px;
            border-radius: 50%;
            background: #ef4444;
            margin-top: 6px;
            shrink: 0;
        }}

        .priority-dot.amber {{
            background: #f59e0b;
        }}

        .gene-symbol-title {{
            font-size: 13.5px;
            font-weight: 700;
            color: var(--text-main);
            line-height: 1.2;
        }}

        .gene-desc-sub {{
            font-size: 11.5px;
            color: var(--text-muted);
            white-space: nowrap;
            overflow: hidden;
            text-overflow: ellipsis;
            max-width: 200px;
            margin-top: 2px;
        }}

        .gene-score-badge {{
            font-family: 'JetBrains Mono', monospace;
            font-size: 12px;
            font-weight: 600;
            color: var(--danger);
            margin-left: 8px;
            shrink: 0;
        }}

        .gene-score-badge.amber {{
            color: var(--warning);
        }}

        /* Right Panel: Content Area */
        .content-panel {{
            flex: 1;
            overflow-y: auto;
            padding: 24px 32px;
            height: 100%;
        }}

        .report-summary-bar {{
            background: #ffffff;
            border: 1px solid var(--border-color);
            border-radius: 8px;
            padding: 12px 20px;
            margin-bottom: 20px;
            display: flex;
            justify-content: space-between;
            align-items: center;
            box-shadow: 0 1px 2px rgba(0, 0, 0, 0.03);
        }}

        .summary-left {{
            display: flex;
            align-items: center;
            gap: 14px;
        }}

        .panel-report-title {{
            font-size: 14px;
            font-weight: 700;
            color: var(--text-main);
        }}

        .detected-variants-badge {{
            background: #64748b;
            color: #ffffff;
            font-size: 12px;
            font-weight: 600;
            padding: 4px 10px;
            border-radius: 6px;
        }}

        .summary-controls {{
            display: flex;
            align-items: center;
            gap: 20px;
            font-size: 12.5px;
            color: var(--text-secondary);
        }}

        .filter-select {{
            background: #ffffff;
            border: 1px solid #cbd5e1;
            border-radius: 6px;
            padding: 4px 8px;
            font-size: 12px;
            color: var(--text-main);
            outline: none;
        }}

        /* Gene Section Card */
        .gene-card {{
            background: var(--bg-surface);
            border: 1px solid var(--border-color);
            border-radius: 10px;
            padding: 22px 24px;
            margin-bottom: 24px;
            box-shadow: 0 1px 3px rgba(0, 0, 0, 0.04);
            scroll-margin-top: 16px;
            transition: all 0.2s ease;
        }}

        .gene-card:hover {{
            border-color: #cbd5e1;
            box-shadow: 0 4px 12px rgba(0, 0, 0, 0.05);
        }}

        .gene-card-header {{
            font-size: 16px;
            font-weight: 800;
            color: var(--text-main);
            padding-bottom: 12px;
            border-bottom: 1px solid var(--border-light);
            margin-bottom: 16px;
            display: flex;
            align-items: center;
            gap: 8px;
        }}

        .gene-full-name-label {{
            font-weight: 500;
            color: var(--text-secondary);
        }}

        /* Two Column Description & Pathology */
        .gene-details-grid {{
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 28px;
            margin-bottom: 20px;
            padding-bottom: 18px;
            border-bottom: 1px solid var(--border-light);
        }}

        .section-subhead {{
            font-size: 11px;
            font-weight: 700;
            text-transform: uppercase;
            letter-spacing: 0.05em;
            color: var(--text-muted);
            margin-bottom: 8px;
        }}

        .gene-description-text {{
            font-size: 12.5px;
            line-height: 1.6;
            color: var(--text-secondary);
        }}

        .source-link {{
            color: var(--primary);
            text-decoration: none;
            font-size: 11.5px;
            margin-top: 6px;
            display: inline-block;
            font-weight: 600;
        }}

        .source-link:hover {{
            text-decoration: underline;
        }}

        .pathology-list {{
            list-style: none;
            display: flex;
            flex-direction: column;
            gap: 8px;
        }}

        .pathology-item {{
            font-size: 12.5px;
            color: var(--text-main);
            display: flex;
            flex-direction: column;
            gap: 4px;
        }}

        .pathology-tags {{
            display: flex;
            gap: 6px;
            flex-wrap: wrap;
        }}

        .badge-inheritance {{
            font-size: 10.5px;
            font-weight: 600;
            padding: 2px 7px;
            border-radius: 4px;
            border: 1px solid #cbd5e1;
            background: #f8fafc;
            color: #334155;
        }}

        .badge-ad {{
            background: #fff7ed;
            border-color: #ffedd5;
            color: #c2410c;
        }}

        .badge-ar {{
            background: #f1f5f9;
            border-color: #e2e8f0;
            color: #475569;
        }}

        .badge-phase {{
            background: #faf5ff;
            border-color: #f3e8ff;
            color: #7e22ce;
        }}

        /* Clinical Variant Table */
        .variant-table-wrapper {{
            width: 100%;
            overflow-x: auto;
        }}

        .variant-table {{
            width: 100%;
            border-collapse: collapse;
            font-size: 12.5px;
            text-align: left;
        }}

        .variant-table th {{
            color: var(--text-muted);
            font-weight: 600;
            font-size: 11px;
            text-transform: uppercase;
            letter-spacing: 0.04em;
            padding: 10px 12px;
            border-bottom: 1px solid var(--border-color);
            background: #f8fafc;
        }}

        .variant-table td {{
            padding: 12px 12px;
            border-bottom: 1px solid var(--border-light);
            color: var(--text-main);
            vertical-align: middle;
        }}

        .variant-table tr:hover td {{
            background: #f8fafc;
        }}

        /* Table Badges & Cells */
        .rs-link {{
            color: var(--primary);
            text-decoration: none;
            font-family: 'JetBrains Mono', monospace;
            font-size: 12px;
            font-weight: 600;
        }}

        .rs-link:hover {{
            text-decoration: underline;
        }}

        .genotype-box-container {{
            display: inline-flex;
            gap: 4px;
        }}

        .allele-box {{
            width: 22px;
            height: 22px;
            border-radius: 4px;
            display: inline-flex;
            align-items: center;
            justify-content: center;
            font-family: 'JetBrains Mono', monospace;
            font-weight: 700;
            font-size: 12px;
        }}

        .allele-box.ref {{
            background: #475569;
            color: #ffffff;
        }}

        .allele-box.alt {{
            background: #d97706;
            color: #ffffff;
        }}

        .maf-badge {{
            font-family: 'JetBrains Mono', monospace;
            font-size: 11.5px;
            font-weight: 600;
            background: #fef9c3;
            color: #854d0e;
            padding: 3px 7px;
            border-radius: 4px;
            display: inline-block;
        }}

        .revel-score {{
            font-family: 'JetBrains Mono', monospace;
            font-size: 12.5px;
            font-weight: 700;
        }}

        .revel-score.high {{
            color: var(--danger);
        }}

        .revel-score.med {{
            color: var(--warning);
        }}

        .revel-score.low {{
            color: var(--success);
        }}

        .impact-pill {{
            font-size: 11.5px;
            font-weight: 600;
            padding: 3px 8px;
            border-radius: 4px;
            display: inline-block;
        }}

        .impact-pill.missense {{
            background: #fee2e2;
            color: #991b1b;
        }}

        .impact-pill.frameshift {{
            background: #fee2e2;
            color: #b91c1c;
        }}

        .impact-pill.nonsense {{
            background: #fee2e2;
            color: #b91c1c;
        }}

        .impact-pill.intron {{
            background: #f1f5f9;
            color: #475569;
            border: 1px solid #e2e8f0;
        }}

        .clinvar-badge {{
            font-size: 11px;
            font-weight: 700;
            padding: 3px 8px;
            border-radius: 4px;
            display: inline-block;
            text-transform: uppercase;
        }}

        .clinvar-badge.pathogenic {{
            background: #dc2626;
            color: #ffffff;
        }}

        .clinvar-badge.likely-pathogenic {{
            background: #ea580c;
            color: #ffffff;
        }}

        .clinvar-badge.vus {{
            background: #fef3c7;
            color: #92400e;
            border: 1px solid #fde68a;
        }}

        .clinvar-badge.likely-benign {{
            background: #dcfce7;
            color: #166534;
            border: 1px solid #bbf7d0;
        }}

        .clinvar-badge.benign {{
            background: #16a34a;
            color: #ffffff;
        }}

        .phase-badge {{
            font-size: 11px;
            font-weight: 600;
            padding: 3px 8px;
            border-radius: 4px;
            display: inline-flex;
            align-items: center;
            gap: 4px;
        }}

        .phase-badge.maternal {{
            background: #f3e8ff;
            color: #6b21a8;
            border: 1px solid #e9d5ff;
        }}

        .phase-badge.paternal {{
            background: #dbeafe;
            color: #1e40af;
            border: 1px solid #bfdbfe;
        }}

        .phase-badge.undetermined {{
            background: #f1f5f9;
            color: #64748b;
        }}

        /* Scrollbars */
        ::-webkit-scrollbar {{
            width: 6px;
            height: 6px;
        }}
        ::-webkit-scrollbar-track {{
            background: transparent;
        }}
        ::-webkit-scrollbar-thumb {{
            background: #cbd5e1;
            border-radius: 3px;
        }}
        ::-webkit-scrollbar-thumb:hover {{
            background: #94a3b8;
        }}

        .hidden {{
            display: none !important;
        }}
    </style>
</head>
<body>

    <!-- Header Navigation -->
    <header class="top-navbar">
        <div class="brand-section">
            <div class="brand-logo">
                <i class="fa-solid fa-dna"></i>
            </div>
            <div>
                <div class="brand-title">Visual Ontology Explorer</div>
                <div class="brand-subtitle">Genomics Report · Clinical Triage & Phased Variants</div>
            </div>
        </div>

        <!-- Global Search -->
        <div class="search-container">
            <i class="fa-solid fa-magnifying-glass search-icon"></i>
            <input type="text" id="mainSearch" class="search-input" placeholder="Search by gene's code, name, rsID, or phenotype..." oninput="onGlobalSearch(this.value)">
        </div>

        <!-- Patient Metadata -->
        <div class="user-meta">
            <div class="patient-pill">
                <div style="font-weight: 700; color: #0f172a;">Patient: {patient_id}</div>
                <div style="color: #64748b; font-size: 11px;">{run_date} | GRCh38</div>
            </div>
            <div class="avatar-circle">DE</div>
        </div>
    </header>

    <!-- Panel Filter Navigation Pills -->
    <div class="panel-nav-bar">
        <button class="nav-pill active" onclick="setPanelFilter('all', this)">All Panels</button>
        <button class="nav-pill" onclick="setPanelFilter('revel', this)">High REVEL (&ge; 0.7)</button>
        <button class="nav-pill" onclick="setPanelFilter('clinvar', this)">ClinVar (P/LP)</button>
        <button class="nav-pill" onclick="setPanelFilter('phased', this)">Phased Variants</button>
        <button class="nav-pill" onclick="setPanelFilter('flagged', this)">Flagged</button>
    </div>

    <!-- Main Workspace (Dual Pane Layout) -->
    <div class="workspace">
        
        <!-- Left Sidebar: Gene Index Panel -->
        <aside class="sidebar-panel">
            <div class="sidebar-header">
                <div class="sidebar-title-row">
                    <div class="sidebar-title" id="sidebarModeTitle">High REVEL score</div>
                    <button class="btn-report-badge" onclick="resetAllFilters()">Panel Report</button>
                </div>
                <input type="text" id="sidebarFilter" class="sidebar-search-box" placeholder="Filter by gene's code..." oninput="onSidebarFilter(this.value)">
            </div>

            <!-- Scrollable Gene List -->
            <ul class="gene-list" id="geneListContainer">
                <!-- Dynamically populated by JS -->
            </ul>
        </aside>

        <!-- Right Main Panel: Stream of Gene Cards -->
        <main class="content-panel" id="mainContentPanel">
            
            <!-- Summary Bar -->
            <div class="report-summary-bar">
                <div class="summary-left">
                    <div class="panel-report-title">Panel Report</div>
                    <div class="detected-variants-badge" id="detectedCountBadge">Detected {total_variants_count} variant(s)</div>
                </div>

                <div class="summary-controls">
                    <label style="display: flex; align-items: center; gap: 6px; cursor: pointer;">
                        <input type="checkbox" id="toggleDescriptions" checked onchange="toggleGeneDescriptions(this.checked)">
                        Gene Descriptions
                    </label>

                    <div style="display: flex; align-items: center; gap: 6px;">
                        <span>Filter by MAF:</span>
                        <select id="mafFilterSelect" class="filter-select" onchange="onMafFilterChange(this.value)">
                            <option value="all">Show all</option>
                            <option value="0.01">&lt; 0.01 (1%)</option>
                            <option value="0.001">&lt; 0.001 (0.1%)</option>
                            <option value="0.0001">&lt; 0.0001 (0.01%)</option>
                        </select>
                    </div>
                </div>
            </div>

            <!-- Gene Cards Container -->
            <div id="geneCardsContainer">
                <!-- Dynamically populated by JS -->
            </div>

        </main>
    </div>

    <!-- Client-Side Script & Data Store -->
    <script>
        const REPORT_DATA = {embedded_data_json};

        let activePanelFilter = 'all';
        let currentSearchQuery = '';
        let currentMafFilter = 'all';
        let showDescriptions = true;

        // Initialize UI
        document.addEventListener('DOMContentLoaded', () => {{
            renderGeneList();
            renderGeneCards();
        }});

        function setPanelFilter(panel, btnElement) {{
            activePanelFilter = panel;
            document.querySelectorAll('.nav-pill').forEach(el => el.classList.remove('active'));
            if (btnElement) btnElement.classList.add('active');

            const titleMap = {{
                'all': 'All Active Genes',
                'revel': 'High REVEL score',
                'clinvar': 'ClinVar Pathogenic / LP',
                'phased': 'Phased Variants (M/P)',
                'flagged': 'Flagged Genes'
            }};
            document.getElementById('sidebarModeTitle').textContent = titleMap[panel] || 'Genes';

            renderGeneList();
            renderGeneCards();
        }}

        function onGlobalSearch(query) {{
            currentSearchQuery = query.toLowerCase().trim();
            renderGeneList();
            renderGeneCards();
        }}

        function onSidebarFilter(query) {{
            currentSearchQuery = query.toLowerCase().trim();
            renderGeneList();
            renderGeneCards();
        }}

        function onMafFilterChange(value) {{
            currentMafFilter = value;
            renderGeneCards();
        }}

        function toggleGeneDescriptions(visible) {{
            showDescriptions = visible;
            document.querySelectorAll('.gene-details-grid').forEach(el => {{
                if (visible) {{
                    el.classList.remove('hidden');
                }} else {{
                    el.classList.add('hidden');
                }}
            }});
        }}

        function resetAllFilters() {{
            activePanelFilter = 'all';
            currentSearchQuery = '';
            currentMafFilter = 'all';
            document.getElementById('mainSearch').value = '';
            document.getElementById('sidebarFilter').value = '';
            document.getElementById('mafFilterSelect').value = 'all';
            document.querySelectorAll('.nav-pill').forEach(el => el.classList.remove('active'));
            document.querySelector('.nav-pill').classList.add('active');
            document.getElementById('sidebarModeTitle').textContent = 'All Active Genes';
            renderGeneList();
            renderGeneCards();
        }}

        function filterGenes() {{
            return REPORT_DATA.genes.filter(gene => {{
                // Panel Filter
                if (activePanelFilter === 'revel' && gene.max_revel < 0.7) return false;
                if (activePanelFilter === 'clinvar' && !gene.has_pathogenic) return false;
                if (activePanelFilter === 'phased' && !gene.has_phased) return false;
                if (activePanelFilter === 'flagged' && !gene.has_pathogenic && gene.max_revel < 0.9) return false;

                // Search Filter
                if (currentSearchQuery) {{
                    const sym = (gene.symbol || '').toLowerCase();
                    const name = (gene.name || '').toLowerCase();
                    const desc = (gene.description || '').toLowerCase();
                    const matchedVar = gene.variants.some(v => 
                        (v.rsid && v.rsid.toLowerCase().includes(currentSearchQuery)) ||
                        (v.impact_consequence && v.impact_consequence.toLowerCase().includes(currentSearchQuery))
                    );
                    const matchedHpo = gene.resolved_hpos.some(h => 
                        h.name.toLowerCase().includes(currentSearchQuery) || h.id.toLowerCase().includes(currentSearchQuery)
                    );
                    if (!sym.includes(currentSearchQuery) && !name.includes(currentSearchQuery) && !desc.includes(currentSearchQuery) && !matchedVar && !matchedHpo) {{
                        return false;
                    }}
                }}

                return true;
            }});
        }}

        function renderGeneList() {{
            const container = document.getElementById('geneListContainer');
            const filtered = filterGenes();
            container.innerHTML = '';

            if (filtered.length === 0) {{
                container.innerHTML = '<li style="padding: 20px; text-align: center; color: #94a3b8; font-size: 12px;">No matching genes found</li>';
                return;
            }}

            filtered.forEach((gene, index) => {{
                const li = document.createElement('li');
                li.className = 'gene-item' + (index === 0 ? ' active' : '');
                li.id = 'sidebar-gene-' + gene.symbol;
                li.onclick = () => scrollToGene(gene.symbol);

                const isHighRevel = gene.max_revel >= 0.7;
                const dotClass = isHighRevel ? 'priority-dot' : 'priority-dot amber';
                const scoreClass = isHighRevel ? 'gene-score-badge' : 'gene-score-badge amber';

                li.innerHTML = `
                    <div class="gene-item-info">
                        <div class="${{dotClass}}"></div>
                        <div style="overflow: hidden;">
                            <div class="gene-symbol-title">${{gene.symbol}}</div>
                            <div class="gene-desc-sub">${{gene.name || gene.description || ''}}</div>
                        </div>
                    </div>
                    <div class="${{scoreClass}}">${{gene.max_revel > 0 ? gene.max_revel.toFixed(3) : '—'}}</div>
                `;
                container.appendChild(li);
            }});
        }}

        function scrollToGene(symbol) {{
            document.querySelectorAll('.gene-item').forEach(el => el.classList.remove('active'));
            const sideEl = document.getElementById('sidebar-gene-' + symbol);
            if (sideEl) sideEl.classList.add('active');

            const card = document.getElementById('card-' + symbol);
            if (card) {{
                card.scrollIntoView({{ behavior: 'smooth', block: 'start' }});
            }}
        }}

        function renderGeneCards() {{
            const container = document.getElementById('geneCardsContainer');
            const filtered = filterGenes();
            container.innerHTML = '';

            let visibleVariantCount = 0;

            if (filtered.length === 0) {{
                container.innerHTML = '<div style="padding: 40px; text-align: center; color: #64748b; background: white; border-radius: 8px; border: 1px solid #e2e8f0;">No matching gene reports found for current filter criteria.</div>';
                document.getElementById('detectedCountBadge').textContent = 'Detected 0 variant(s)';
                return;
            }}

            filtered.forEach(gene => {{
                // Filter variants by MAF if selected
                const variants = gene.variants.filter(v => {{
                    if (currentMafFilter === 'all') return true;
                    const maxAf = parseFloat(currentMafFilter);
                    const af = v.gnomad_af !== undefined && v.gnomad_af !== null ? parseFloat(v.gnomad_af) : 0;
                    return af <= maxAf;
                }});

                visibleVariantCount += variants.length;

                const card = document.createElement('div');
                card.className = 'gene-card';
                card.id = 'card-' + gene.symbol;

                // Pathologies HTML
                let pathologyHtml = '';
                if (gene.pathologies && gene.pathologies.length > 0) {{
                    const items = gene.pathologies.map(p => {{
                        const badges = (p.inheritance || []).map(inh => {{
                            let cls = 'badge-inheritance';
                            if (inh.includes('Dominant')) cls += ' badge-ad';
                            else if (inh.includes('Recessive')) cls += ' badge-ar';
                            return `<span class="${{cls}}">${{inh}}</span>`;
                        }}).join(' ');
                        return `
                            <li class="pathology-item">
                                <div>• ${{p.name}}</div>
                                <div class="pathology-tags">${{badges}}</div>
                            </li>
                        `;
                    }}).join('');
                    pathologyHtml = `<ul class="pathology-list">${{items}}</ul>`;
                }} else if (gene.resolved_hpos && gene.resolved_hpos.length > 0) {{
                    const items = gene.resolved_hpos.map(h => `
                        <li class="pathology-item">
                            <div>• ${{h.name}}</div>
                            <div class="pathology-tags"><span class="badge-inheritance badge-ar">${{h.id}}</span></div>
                        </li>
                    `).join('');
                    pathologyHtml = `<ul class="pathology-list">${{items}}</ul>`;
                }} else {{
                    pathologyHtml = '<div style="font-size: 12px; color: #94a3b8;">No specific OMIM pathology annotations recorded.</div>';
                }}

                // Variant Rows HTML
                let variantRowsHtml = '';
                variants.forEach(v => {{
                    // Genotype boxes
                    const gt = v.genotype || 'N/A';
                    let gtBoxes = '';
                    if (gt.includes('/')) {{
                        const alleles = gt.split('/');
                        gtBoxes = `
                            <div class="genotype-box-container">
                                <span class="allele-box ref">${{alleles[0]}}</span>
                                <span class="allele-box alt">${{alleles[1]}}</span>
                            </div>
                        `;
                    }} else {{
                        gtBoxes = `<span class="allele-box alt">${{gt}}</span>`;
                    }}

                    // MAF
                    const afVal = v.gnomad_af !== undefined && v.gnomad_af !== null ? Number(v.gnomad_af).toFixed(6).replace(/0+$/, '').replace(/\\.$/, '') : '—';
                    const mafHtml = afVal !== '—' ? `<span class="maf-badge">${{afVal}}</span>` : '<span style="color:#94a3b8;">—</span>';

                    // REVEL
                    const rev = v.revel_score !== undefined && v.revel_score !== null ? parseFloat(v.revel_score).toFixed(3) : '—';
                    let revCls = 'revel-score';
                    if (rev !== '—') {{
                        const num = parseFloat(rev);
                        if (num >= 0.7) revCls += ' high';
                        else if (num >= 0.5) revCls += ' med';
                        else revCls += ' low';
                    }}

                    // Impact
                    const impact = v.impact_consequence || '—';
                    let impactCls = 'impact-pill';
                    if (impact.toLowerCase().includes('missense')) impactCls += ' missense';
                    else if (impact.toLowerCase().includes('frameshift')) impactCls += ' frameshift';
                    else if (impact.toLowerCase().includes('nonsense')) impactCls += ' nonsense';
                    else impactCls += ' intron';

                    // ClinVar
                    const clinSig = v.clinvar_significance || 'VUS';
                    let clinCls = 'clinvar-badge vus';
                    const sigLow = clinSig.toLowerCase();
                    if (sigLow.includes('likely pathogenic')) clinCls = 'clinvar-badge likely-pathogenic';
                    else if (sigLow.includes('pathogenic')) clinCls = 'clinvar-badge pathogenic';
                    else if (sigLow.includes('likely benign')) clinCls = 'clinvar-badge likely-benign';
                    else if (sigLow.includes('benign')) clinCls = 'clinvar-badge benign';

                    // Phase
                    const phase = (v.phasing || 'undetermined').toLowerCase();
                    let phaseBadge = '<span class="phase-badge undetermined">Unknown</span>';
                    if (phase === 'maternal') {{
                        phaseBadge = '<span class="phase-badge maternal"><i class="fa-solid fa-venus"></i> Maternal</span>';
                    }} else if (phase === 'paternal') {{
                        phaseBadge = '<span class="phase-badge paternal"><i class="fa-solid fa-mars"></i> Paternal</span>';
                    }}

                    // rsID link
                    const rsLink = v.rsid ? `<a href="https://www.ncbi.nlm.nih.gov/snp/${{v.rsid}}" target="_blank" class="rs-link">${{v.rsid}}</a>` : `<span style="color:#64748b; font-family:monospace;">${{v.chromosome}}:${{v.position}}</span>`;

                    variantRowsHtml += `
                        <tr>
                            <td>${{rsLink}}</td>
                            <td>${{gtBoxes}}</td>
                            <td>${{mafHtml}}</td>
                            <td><span class="${{revCls}}">${{rev}}</span></td>
                            <td><span class="${{impactCls}}">${{impact}}</span></td>
                            <td><span class="${{clinCls}}">${{clinSig}}</span></td>
                            <td>${{phaseBadge}}</td>
                        </tr>
                    `;
                }});

                card.innerHTML = `
                    <div class="gene-card-header">
                        <span>${{gene.symbol}}</span>
                        <span class="gene-full-name-label">— ${{gene.name || ''}}</span>
                    </div>

                    <div class="gene-details-grid ${{showDescriptions ? '' : 'hidden'}}">
                        <div>
                            <div class="section-subhead">Gene Description</div>
                            <div class="gene-description-text">${{gene.description || 'No description available.'}}</div>
                            <div>
                                <a href="https://www.omim.org/entry/${{gene.omim_source ? gene.omim_source.replace('OMIM:', '') : ''}}" target="_blank" class="source-link">(Source: ${{gene.omim_source || 'OMIM'}})</a>
                            </div>
                        </div>

                        <div>
                            <div class="section-subhead">Associated Pathology</div>
                            ${{pathologyHtml}}
                        </div>
                    </div>

                    <div class="variant-table-wrapper">
                        <table class="variant-table">
                            <thead>
                                <tr>
                                    <th>rs ID</th>
                                    <th>Genotype</th>
                                    <th>MAF</th>
                                    <th>REVEL</th>
                                    <th>Impact of the variant</th>
                                    <th>ClinVar records</th>
                                    <th>Phasing</th>
                                </tr>
                            </thead>
                            <tbody>
                                ${{variantRowsHtml}}
                            </tbody>
                        </table>
                    </div>
                `;

                container.appendChild(card);
            }});

            document.getElementById('detectedCountBadge').textContent = `Detected ${{visibleVariantCount}} variant(s)`;
        }}
    </script>
</body>
</html>"""

    with open(output_filepath, "w", encoding="utf-8") as f:
        f.write(html_content)
    
    print(f"Successfully generated Visual Ontology Explorer report at: {output_filepath}")
    return output_filepath

def main():
    parser = argparse.ArgumentParser(description="Visual Ontology Explorer HTML Report Generator")
    parser.add_argument("-i", "--input", help="Path to input JSON file containing variant report data")
    parser.add_argument("-o", "--output", default="reports/visual_ontology_explorer.html", help="Path to output HTML file")
    parser.add_argument("-d", "--demo", "--mock", action="store_true", help="Generate report using demo sample dataset")
    args = parser.parse_args()

    out_dir = os.path.dirname(args.output)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    if args.demo:
        demo_path = Path(__file__).parent / "samples" / "demo_variant_report.json"
        if demo_path.exists():
            with open(demo_path, "r", encoding="utf-8") as f:
                data = json.load(f)
        else:
            data = {
                "patient_id": "DE_WGS_2026",
                "run_date": "2026-08-21 18:30 UTC",
                "monogenic_findings": []
            }
        generate_report_html(data, args.output)
    else:
        if not args.input:
            parser.print_help()
            sys.exit(1)
        if not os.path.exists(args.input):
            print(f"Error: Input file does not exist: {args.input}")
            sys.exit(1)
        with open(args.input, "r", encoding="utf-8") as f:
            data = json.load(f)
        generate_report_html(data, args.output)

if __name__ == "__main__":
    main()
