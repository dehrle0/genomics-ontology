#!/usr/bin/env python3
"""
render_new_ontology_report.py
Render a modern Genomic & Polygenic Ontology Dashboard report with:
- Top metrics summary (Screened, Tier 1, Tier 2, High PGS traits)
- Polygenic trait risk overview bars
- Interactive live search & dynamic ontology category filter
- Gene cards with NCBI summary & detailed SNP tables
- Full TSV, JSON, and text exports
"""
import json
import os
import sys
import html
import argparse
from typing import Dict, List, Any

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import render_report as rr

SO_NAME = rr.SO_NAME

HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>Genomic Ontology & Polygenic Report - {sample_id}</title>
    <style>
        :root {{
            --bg-color: #f8fafc;
            --card-bg: #ffffff;
            --text-main: #0f172a;
            --text-muted: #64748b;
            --border-color: #e2e8f0;
            --red-bg: #fee2e2; --red-text: #991b1b; --red-border: #fca5a5;
            --yellow-bg: #fef3c7; --yellow-text: #92400e; --yellow-border: #fcd34d;
            --green-bg: #dcfce7; --green-text: #166534; --green-border: #86efac;
            --blue-accent: #2563eb;
        }}

        * {{ box-sizing: border-box; margin: 0; padding: 0; }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
            background-color: var(--bg-color);
            color: var(--text-main);
            padding: 24px;
            line-height: 1.5;
        }}

        .container {{ max-width: 1300px; margin: 0 auto; }}

        /* Header & Dashboard */
        .header-card {{
            background: var(--card-bg);
            border: 1px solid var(--border-color);
            border-radius: 12px;
            padding: 28px;
            margin-bottom: 24px;
            box-shadow: 0 1px 3px rgba(0,0,0,0.05);
        }}
        .header-card h2 {{ font-size: 1.6rem; font-weight: 800; color: var(--text-dark); }}
        .header-card p {{ color: var(--text-muted); font-size: 0.9rem; margin-top: 4px; }}

        .metrics-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
            gap: 16px;
            margin-top: 20px;
        }}

        .metric-box {{
            background: #f1f5f9;
            padding: 16px 20px;
            border-radius: 8px;
            text-align: center;
        }}
        .metric-box .number {{ font-size: 1.8rem; font-weight: 800; color: var(--blue-accent); }}
        .metric-box .label {{ font-size: 0.75rem; color: var(--text-muted); text-transform: uppercase; font-weight: 600; letter-spacing: 0.05em; margin-top: 4px; }}

        /* Polygenic Risk Summary Bar Charts */
        .pgs-section {{
            margin-top: 24px;
            padding-top: 18px;
            border-top: 1px solid var(--border-color);
        }}
        .pgs-section h3 {{ font-size: 1.1rem; font-weight: 700; margin-bottom: 12px; }}
        .pgs-row {{
            display: flex;
            align-items: center;
            margin-bottom: 10px;
        }}
        .pgs-trait {{ width: 240px; font-weight: 600; font-size: 0.85rem; }}
        .pgs-bar-container {{ flex-grow: 1; background: #e2e8f0; height: 18px; border-radius: 9px; overflow: hidden; position: relative; margin: 0 12px; }}
        .pgs-bar {{ height: 100%; border-radius: 9px; transition: width 0.4s ease; }}
        .pgs-value {{ width: 110px; font-size: 0.85rem; font-weight: 700; text-align: right; }}

        /* Interactive Filter Bar */
        .search-bar-container {{
            position: sticky;
            top: 12px;
            z-index: 100;
            background: var(--card-bg);
            padding: 16px 20px;
            border: 1px solid var(--border-color);
            border-radius: 10px;
            margin-bottom: 24px;
            box-shadow: 0 4px 12px rgba(0,0,0,0.06);
            display: flex;
            gap: 12px;
            flex-wrap: wrap;
        }}

        .search-input, .select-input {{
            padding: 10px 14px;
            border: 1px solid var(--border-color);
            border-radius: 6px;
            font-size: 0.9rem;
            outline: none;
        }}
        .search-input {{ flex-grow: 1; min-width: 260px; }}
        .search-input:focus, .select-input:focus {{ border-color: var(--blue-accent); }}

        /* Gene Card & Hierarchy */
        .gene-card {{
            background: var(--card-bg);
            border: 1px solid var(--border-color);
            border-radius: 12px;
            margin-bottom: 24px;
            overflow: hidden;
            box-shadow: 0 1px 3px rgba(0,0,0,0.04);
        }}

        .gene-header {{
            background: #f8fafc;
            padding: 18px 24px;
            border-bottom: 1px solid var(--border-color);
        }}

        .ncbi-gene-title {{
            font-size: 1.2rem;
            font-weight: 700;
            margin: 0 0 8px 0;
            color: var(--text-main);
            display: flex;
            align-items: center;
            flex-wrap: wrap;
            gap: 8px;
        }}

        .ncbi-description {{
            font-size: 0.95rem; 
            line-height: 1.5;
            color: #334155;
            margin-top: 6px;
        }}

        .ontology-badge {{
            display: inline-block;
            background: #e0f2fe;
            color: #0369a1;
            font-size: 0.75rem;
            padding: 3px 8px;
            border-radius: 6px;
            font-weight: 600;
            border: 1px solid #bae6fd;
        }}

        /* Standardized SNP Table */
        .snp-table {{
            width: 100%;
            border-collapse: collapse;
            font-size: 0.85rem;
        }}

        .snp-table th {{
            background: #f8fafc;
            text-align: left;
            padding: 12px 16px;
            border-bottom: 1px solid var(--border-color);
            color: var(--text-muted);
            font-weight: 600;
            text-transform: uppercase;
            font-size: 0.75rem;
            letter-spacing: 0.04em;
        }}

        .snp-table td {{
            padding: 12px 16px;
            border-bottom: 1px solid #f1f5f9;
            vertical-align: middle;
        }}

        /* Semantic Color Highlights */
        .status-red {{ background-color: var(--red-bg) !important; color: var(--red-text); font-weight: 700; border-radius: 4px; padding: 2px 6px; }}
        .status-yellow {{ background-color: var(--yellow-bg) !important; color: var(--yellow-text); font-weight: 700; border-radius: 4px; padding: 2px 6px; }}
        .status-green {{ background-color: var(--green-bg) !important; color: var(--green-text); font-weight: 700; border-radius: 4px; padding: 2px 6px; }}

        @media print {{
            .search-bar-container {{ display: none; }}
            body {{ background: white; padding: 0; }}
            .gene-card {{ border: 1px solid #ccc; break-inside: avoid; }}
        }}
    </style>
</head>
<body>

<div class="container">
    <div class="header-card">
        <h2>Genomic & Polygenic Ontology Summary</h2>
        <p><strong>Sample ID:</strong> {sample_id} | <strong>Domain:</strong> {domain}</p>
        
        <div class="metrics-grid">
            <div class="metric-box">
                <div class="number">{metrics_analyzed:,}</div>
                <div class="label">Variants Screened</div>
            </div>
            <div class="metric-box">
                <div class="number" style="color: #dc2626;">{metrics_tier1}</div>
                <div class="label">Tier 1 Pathogenic</div>
            </div>
            <div class="metric-box">
                <div class="number" style="color: #d97706;">{metrics_tier2}</div>
                <div class="label">Tier 2 VUS</div>
            </div>
            <div class="metric-box" style="display: {pgs_display};">
                <div class="number" style="color: #2563eb;">{metrics_high_pgs}</div>
                <div class="label">High Polygenic Traits</div>
            </div>
        </div>

        <div class="pgs-section" style="display: {pgs_display};">
            <h3>Polygenic Trait Risk Overview</h3>
            {pgs_rows_html}
        </div>
    </div>

    <div class="search-bar-container">
        <input type="text" id="searchInput" class="search-input" onkeyup="filterReport()" placeholder="Search by Gene (e.g. LDLR), RSID (rs121908025), Reason Code, or Keyword...">
        <select id="ontologyFilter" class="select-input" onchange="filterReport()">
            <option value="">All Ontology Categories</option>
            {ontology_options_html}
        </select>
    </div>

    <div id="genesContainer">
        {genes_html}
    </div>
</div>

<script>
function filterReport() {{
    let searchVal = document.getElementById('searchInput').value.toLowerCase();
    let ontologyVal = document.getElementById('ontologyFilter').value.toLowerCase();
    let geneCards = document.getElementsByClassName('gene-card');

    for (let card of geneCards) {{
        let textContent = card.innerText.toLowerCase();
        let matchesSearch = textContent.includes(searchVal);
        let matchesOntology = ontologyVal === "" || textContent.includes(ontologyVal);

        if (matchesSearch && matchesOntology) {{
            card.style.display = "";
        }} else {{
            card.style.display = "none";
        }}
    }}
}}
</script>

</body>
</html>
"""

def build_pgs_rows(pgs_list: List[Dict[str, Any]]) -> str:
    if not pgs_list:
        return ""
    html_out = ""
    for item in pgs_list:
        perc = item.get("percentile", 50)
        color = "#dc2626" if perc >= 80 else ("#16a34a" if perc <= 10 else "#2563eb")
        html_out += f"""
        <div class="pgs-row">
            <div class="pgs-trait">{html.escape(item.get('trait', 'Unknown'))} <br><small style="color:#64748b;">{html.escape(item.get('pgs_id', ''))}</small></div>
            <div class="pgs-bar-container">
                <div class="pgs-bar" style="width: {perc}%; background-color: {color};"></div>
            </div>
            <div class="pgs-value">{perc:.1f}th percentile</div>
        </div>
        """
    return html_out

def build_genes_sections(genes_list: List[Dict[str, Any]]) -> str:
    html_out = ""
    for gene in genes_list:
        hugo = html.escape(gene['gene_symbol'])
        ncbi_link = f"<a href='https://www.ncbi.nlm.nih.gov/gene/{gene['ncbi_gene_id']}' target='_blank' style='color:#2563eb;'>{gene['ncbi_gene_id']}</a>" if gene.get('ncbi_gene_id') else ""
        ncbi_id_span = f"<span style='font-size: 0.85rem; font-weight: normal; color:#64748b;'>(NCBI ID: {ncbi_link})</span>" if ncbi_link else ""
        
        html_out += f"""
        <div class="gene-card">
            <div class="gene-header">
                <div class="ncbi-gene-title">
                    Gene: <strong>{hugo}</strong> 
                    {ncbi_id_span}
                    <span class="ontology-badge">{html.escape(gene.get('ontology_category', ''))}</span>
                    <span class="ontology-badge" style="background: #f1f5f9; color: #475569; border-color:#cbd5e1;">{html.escape(gene['tier'])}</span>
                </div>
                <div class="ncbi-description">
                    <strong>NCBI Gene Summary:</strong> {html.escape(gene['ncbi_summary'])}
                </div>
            </div>
            
            <table class="snp-table">
                <thead>
                    <tr>
                        <th>Variant (RSID / Pos)</th>
                        <th>Consequence</th>
                        <th>Zygosity</th>
                        <th>Clinical Effect</th>
                        <th>OpenCRAVAT Predictors</th>
                        <th>Reason Codes</th>
                    </tr>
                </thead>
                <tbody>
        """
        for snp in gene.get("snps", []):
            color_class = f"status-{snp['status_color']}"
            rsid_txt = html.escape(str(snp['rsid']))
            if str(snp['rsid']).startswith("rs"):
                rsid_txt = f"<a href='https://www.ncbi.nlm.nih.gov/snp/{snp['rsid']}' target='_blank' style='color:#2563eb;'>{rsid_txt}</a>"
                
            html_out += f"""
                    <tr>
                        <td><strong>{rsid_txt}</strong><br><small style="color:#64748b; font-family:monospace;">{html.escape(snp['chrom_pos'])}</small></td>
                        <td><span style="background:#f1f5f9; padding:2px 6px; border-radius:4px; font-size:0.75rem;">{html.escape(snp['consequence'])}</span></td>
                        <td><strong>{html.escape(str(snp['zygosity']))}</strong></td>
                        <td><span class="{color_class}">{html.escape(str(snp['effect']))}</span></td>
                        <td>{html.escape(str(snp['opencravat_score']))}</td>
                        <td style="font-size:0.75rem; color:#64748b;">{html.escape(str(snp['reasons']))}</td>
                    </tr>
            """
        html_out += """
                </tbody>
            </table>
        </div>
        """
    return html_out

def adapt_pipeline_data(data: Dict[str, Any]) -> Dict[str, Any]:
    """Convert pipeline JSON into the format expected by the dashboard renderer."""
    adapted = {
        "sample_id": data.get("patient", "UNKNOWN"),
        "domain": data.get("domain", "general"),
        "summary_metrics": {
            "total_variants_analyzed": data.get("scanned_panel_variants", 0),
            "tier_1_pathogenic": data.get("tier_counts", {}).get("Tier1", 0),
            "tier_2_vus": data.get("tier_counts", {}).get("Tier2", 0),
            "polygenic_high_risk_traits": 0
        },
        "polygenic_scores": [], 
        "genes": [],
        "categories": set()
    }
    
    genes_map = {}
    for r in data.get("records", []):
        hugo = r.get("hugo", "UNKNOWN")
        if hugo not in genes_map:
            gi = r.get("gene_info") or {}
            ev = r.get("evidence", {})
            hpo_ctx = ", ".join(ev.get("hpo_context", []) or [])
            onto_cat = hpo_ctx if hpo_ctx else adapted["domain"].capitalize()
            if onto_cat:
                adapted["categories"].add(onto_cat)
            
            genes_map[hugo] = {
                "gene_symbol": hugo,
                "ncbi_gene_id": gi.get("ncbi_gene_id", ""),
                "ncbi_summary": gi.get("description", "No summary available."),
                "ontology_category": onto_cat,
                "tier": "Tier 3 - Monitor" if r.get("tier") == "Tier3" else ("Tier 2 - VUS" if r.get("tier") == "Tier2" else "Tier 1 - Actionable Pathogenic"),
                "snps": []
            }
        
        tier_val = r.get("tier", "Tier3")
        if tier_val == "Tier1":
             genes_map[hugo]["tier"] = "Tier 1 - Actionable Pathogenic"
        elif tier_val == "Tier2" and "Tier 1" not in genes_map[hugo]["tier"]:
             genes_map[hugo]["tier"] = "Tier 2 - VUS"
             
        rsid = r.get("rsid") or "-"
        chrom_pos = f"{r.get('chrom')}:{r.get('pos')} {r.get('ref')}>{r.get('alt')}"
        so_code = r.get("so")
        consequence = SO_NAME.get(so_code, so_code or "-")
        
        ev = r.get("evidence", {})
        zyg = ev.get("zygosity") or "-"
        vaf = ev.get("vaf")
        if vaf is not None:
            zyg += f" (VAF: {vaf:.2f})"
        
        if tier_val == "Tier1":
            color = "red"
        elif tier_val == "Tier2":
            color = "yellow"
        else:
            color = "green"
            
        effect = r.get("clinvar_sig") or " | ".join(r.get("reason_codes", []))
        if not effect:
            effect = "Uncertain"
            
        revel = r.get("revel")
        am = r.get("am_path")
        cadd = r.get("cadd_phred")
        preds = []
        if revel: preds.append(f"REVEL: {revel}")
        if am: preds.append(f"AM: {am}")
        if cadd: preds.append(f"CADD: {cadd}")
        oc_score = " | ".join(preds) if preds else "-"
        
        genes_map[hugo]["snps"].append({
            "rsid": rsid,
            "chrom_pos": chrom_pos,
            "consequence": consequence,
            "zygosity": zyg,
            "effect": effect,
            "opencravat_score": oc_score,
            "reasons": " | ".join(r.get("reason_codes", [])),
            "status_color": color
        })
    
    adapted["genes"] = sorted(genes_map.values(), key=lambda g: g["tier"])
    return adapted

def generate_report(data: Dict[str, Any], output_filename: str):
    pgs_html = build_pgs_rows(data.get("polygenic_scores", []))
    genes_html = build_genes_sections(data.get("genes", []))
    
    metrics = data.get("summary_metrics", {})
    pgs_display = "none" if not data.get("polygenic_scores", []) else "block"
    
    categories = sorted(list(data.get("categories", set())))
    options_html = "".join(f"<option value='{html.escape(c)}'>{html.escape(c)}</option>" for c in categories)
    
    full_html = HTML_TEMPLATE.format(
        sample_id=html.escape(data.get("sample_id", "N/A")),
        domain=html.escape(data.get("domain", "General")),
        metrics_analyzed=metrics.get("total_variants_analyzed", 0),
        metrics_tier1=metrics.get("tier_1_pathogenic", 0),
        metrics_tier2=metrics.get("tier_2_vus", 0),
        metrics_high_pgs=metrics.get("polygenic_high_risk_traits", 0),
        pgs_display=pgs_display,
        pgs_rows_html=pgs_html,
        ontology_options_html=options_html,
        genes_html=genes_html
    )
    
    with open(output_filename, "w", encoding="utf-8") as f:
        f.write(full_html)
    print(f"[render] Dashboard HTML -> {os.path.abspath(output_filename)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Render interactive ontology dashboard report")
    parser.add_argument("--in-json", required=True)
    parser.add_argument("--out-html", required=True)
    parser.add_argument("--out-tsv", required=True)
    parser.add_argument("--out-text", required=True)
    parser.add_argument("--split-reports", action="store_true", help="Emit separate monogenic and polygenic reports")
    args = parser.parse_args()
    
    raw_data = json.load(open(args.in_json))
    adapted_data = adapt_pipeline_data(raw_data)
    
    if args.split_reports:
        mono_html = args.out_html.replace(".html", "_monogenic.html")
        adapted_data["polygenic_scores"] = []
        generate_report(adapted_data, mono_html)
        
        poly_html = args.out_html.replace(".html", "_polygenic.html")
        generate_report(adapted_data, poly_html)
    else:
        generate_report(adapted_data, args.out_html)
        
    rr.write_tsv(raw_data.get("records", []), args.out_tsv)
    rr.write_text(raw_data, args.out_text)
    print(f"[render] TSV -> {os.path.abspath(args.out_tsv)}")
    print(f"[render] text -> {os.path.abspath(args.out_text)}")

