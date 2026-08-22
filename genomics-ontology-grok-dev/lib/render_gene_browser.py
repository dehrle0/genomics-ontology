#!/usr/bin/env python3
"""
render_gene_browser.py
Render actionable variants in a two-pane Genomic Variant Browser layout:
- Left sidebar: Flagged gene cards with tier-colored variant counters and quick filtering.
- Right pane: Rich gene detail pane with bold NCBI summary, OMIM Clinical Synopsis, phase status, clinical ontology links, and variant tables.
- Full TSV, JSON, and structured text summary exporters.
"""
import argparse
import html
import json
import os
import sys
from typing import Dict, Any, List

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import render_report as rr

SO_NAME = rr.SO_NAME

HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8"/>
    <meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no"/>
    <title>Genomic Variant Browser — {sample_id} ({domain})</title>
    <style>
        :root {{
            --bg-main: #f8fafc;
            --sidebar-bg: #ffffff;
            --text-dark: #0f172a;
            --text-muted: #64748b;
            --border-color: #e2e8f0;
            --primary: #2563eb;
            --primary-hover: #1d4ed8;
            --tier1-bg: #fee2e2; --tier1-txt: #991b1b; --tier1-border: #fca5a5;
            --tier2-bg: #fef3c7; --tier2-txt: #92400e; --tier2-border: #fcd34d;
            --tier3-bg: #e0f2fe; --tier3-txt: #075985; --tier3-border: #7dd3fc;
        }}
        * {{ box-sizing: border-box; margin: 0; padding: 0; }}
        body, html {{ height: 100%; font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif; background: var(--bg-main); color: var(--text-dark); overflow: hidden; }}
        
        .top-navbar {{
            height: 56px;
            background: linear-gradient(135deg, #1e293b 0%, #0f172a 100%);
            color: #ffffff;
            display: flex;
            align-items: center;
            justify-content: space-between;
            padding: 0 24px;
            border-bottom: 1px solid #334155;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .top-navbar .brand {{ font-size: 1.15rem; font-weight: 700; display: flex; align-items: center; gap: 8px; }}
        .top-navbar .meta {{ font-size: 0.85rem; color: #cbd5e1; }}
        .top-navbar .meta strong {{ color: #ffffff; }}
        
        .app-layout {{ display: flex; height: calc(100% - 56px); width: 100%; }}
        
        /* Left Sidebar */
        .sidebar {{
            width: 320px;
            min-width: 320px;
            background: var(--sidebar-bg);
            border-right: 1px solid var(--border-color);
            display: flex;
            flex-direction: column;
            height: 100%;
        }}
        .sidebar-header {{
            padding: 16px;
            border-bottom: 1px solid var(--border-color);
            background: #ffffff;
        }}
        .sidebar-header .title {{ font-size: 0.95rem; font-weight: 700; color: var(--text-dark); }}
        .sidebar-header .subtitle {{ font-size: 0.75rem; color: var(--text-muted); margin-top: 2px; }}

        .sidebar-search {{ padding: 10px 16px; border-bottom: 1px solid var(--border-color); background: #f8fafc; }}
        .sidebar-search input {{
            width: 100%;
            padding: 8px 12px;
            border: 1px solid var(--border-color);
            border-radius: 6px;
            font-size: 0.85rem;
            outline: none;
        }}
        .sidebar-search input:focus {{ border-color: var(--primary); }}

        .sidebar-list {{
            flex: 1;
            overflow-y: auto;
            padding: 8px;
        }}
        .gene-item {{
            padding: 10px 14px;
            border-radius: 8px;
            margin-bottom: 4px;
            cursor: pointer;
            border: 1px solid transparent;
            transition: all 0.15s ease;
            display: flex;
            flex-direction: column;
            gap: 4px;
        }}
        .gene-item:hover {{ background: #f1f5f9; }}
        .gene-item.active {{ background: #eff6ff; border-color: #bfdbfe; }}
        .gene-item-head {{ display: flex; justify-content: space-between; align-items: center; }}
        .gene-item-name {{ font-weight: 800; font-size: 0.95rem; color: var(--text-dark); }}
        .gene-item.active .gene-item-name {{ color: var(--primary); }}
        .gene-item-desc {{ font-size: 0.75rem; color: var(--text-muted); white-space: nowrap; overflow: hidden; text-overflow: ellipsis; }}
        
        .badge-group {{ display: flex; gap: 4px; }}
        .gene-item-badge {{
            font-size: 0.7rem;
            font-weight: 700;
            padding: 2px 6px;
            border-radius: 10px;
        }}
        .badge-tier1 {{ background: var(--tier1-bg); color: var(--tier1-txt); border: 1px solid var(--tier1-border); }}
        .badge-tier2 {{ background: var(--tier2-bg); color: var(--tier2-txt); border: 1px solid var(--tier2-border); }}
        .badge-tier3 {{ background: var(--tier3-bg); color: var(--tier3-txt); border: 1px solid var(--tier3-border); }}
        
        /* Right Content Area */
        .main-pane {{
            flex: 1;
            height: 100%;
            overflow-y: auto;
            padding: 28px 36px;
            background: var(--bg-main);
        }}
        .gene-content-pane {{ display: none; }}
        .gene-content-pane.active {{ display: block; }}
        
        .gene-hero {{
            background: #ffffff;
            border: 1px solid var(--border-color);
            border-radius: 12px;
            padding: 22px 26px;
            margin-bottom: 20px;
            box-shadow: 0 1px 3px rgba(0,0,0,0.03);
        }}
        .gene-hero-title {{ font-size: 1.6rem; font-weight: 800; color: var(--text-dark); display: flex; align-items: center; gap: 12px; }}
        .gene-hero-summary {{ font-size: 0.95rem; font-weight: 600; color: #1e293b; margin: 10px 0; line-height: 1.5; }}
        
        .omim-browser-box {{
            font-size: 0.85rem;
            color: #475569;
            background: #faf5ff;
            border-left: 4px solid #7c3aed;
            padding: 8px 12px;
            border-radius: 6px;
            margin-bottom: 12px;
        }}

        .gene-links {{ display: flex; gap: 10px; flex-wrap: wrap; margin-top: 10px; }}
        .gene-links a {{
            font-size: 0.8rem;
            color: var(--primary);
            text-decoration: none;
            font-weight: 600;
            background: #f1f5f9;
            padding: 4px 10px;
            border-radius: 6px;
        }}
        .gene-links a:hover {{ background: #e2e8f0; }}
        
        /* SNP Table */
        .table-wrap {{
            background: #ffffff;
            border: 1px solid var(--border-color);
            border-radius: 12px;
            overflow: hidden;
            box-shadow: 0 1px 3px rgba(0,0,0,0.03);
            margin-top: 16px;
        }}
        table.variant-table {{ width: 100%; border-collapse: collapse; text-align: left; font-size: 0.85rem; }}
        table.variant-table th {{
            background: #f1f5f9;
            color: #475569;
            font-weight: 700;
            padding: 10px 14px;
            border-bottom: 1px solid var(--border-color);
            text-transform: uppercase;
            font-size: 0.72rem;
            letter-spacing: 0.05em;
        }}
        table.variant-table td {{
            padding: 12px 14px;
            border-bottom: 1px solid #f1f5f9;
            color: #334155;
            vertical-align: middle;
        }}
        table.variant-table tr:hover td {{ background: #f8fafc; }}
        
        .badge-so {{ background: #f1f5f9; padding: 2px 6px; border-radius: 4px; font-size: 0.75rem; font-weight: 600; }}
        .badge-zyg {{ padding: 2px 6px; border-radius: 4px; font-size: 0.75rem; font-weight: 600; }}
        .zyg-het {{ background: #fef3c7; color: #92400e; }}
        .zyg-hom {{ background: #fee2e2; color: #991b1b; }}
        .zyg-hemi {{ background: #f3e8ff; color: #6b21a8; }}
        
        .reason-tag {{ display: inline-block; font-size: 0.7rem; padding: 1px 6px; border-radius: 4px; margin: 1px; background: #e0f2fe; color: #0369a1; }}
        .reason-tag.strong {{ background: #fee2e2; color: #991b1b; font-weight: 700; }}
        .reason-tag.pheno {{ background: #dcfce7; color: #166534; }}
        
        @media print {{
            body, html {{ overflow: visible; height: auto; }}
            .top-navbar, .sidebar, .sidebar-search {{ display: none !important; }}
            .app-layout {{ height: auto; width: 100%; }}
            .main-pane {{ padding: 0; overflow: visible; }}
            .gene-content-pane {{ display: block !important; page-break-after: always; }}
        }}
    </style>
</head>
<body>

<div class="top-navbar">
    <div class="brand">🧬 Genomic Variant Browser</div>
    <div class="meta">
        Sample: <strong>{sample_id}</strong> &nbsp;|&nbsp;
        Domain: <strong>{domain}</strong> &nbsp;|&nbsp;
        Actionable Genes: <strong>{total_genes}</strong> &nbsp;|&nbsp;
        Actionable Variants: <strong>{total_variants}</strong>
    </div>
</div>

<div class="app-layout">
    <!-- Left Sidebar -->
    <div class="sidebar">
        <div class="sidebar-header">
            <div class="title">Flagged Genes ({total_genes})</div>
            <div class="subtitle">{total_variants} actionable variants across tiers</div>
        </div>
        <div class="sidebar-search">
            <input type="text" id="geneSearch" placeholder="Filter genes (e.g. MYH7)..." oninput="filterGeneSidebar(this.value)">
        </div>
        <div class="sidebar-list" id="sidebarList">
            {sidebar_html}
        </div>
    </div>

    <!-- Right Detail Pane -->
    <div class="main-pane" id="mainPane">
        {main_html}
    </div>
</div>

<script>
    function showGene(hugo) {{
        document.querySelectorAll('.gene-content-pane').forEach(function(p) {{ p.classList.remove('active'); }});
        document.querySelectorAll('.gene-item').forEach(function(i) {{ i.classList.remove('active'); }});
        
        var pane = document.getElementById('pane-' + hugo);
        if (pane) pane.classList.add('active');
        
        var nav = document.getElementById('nav-' + hugo);
        if (nav) nav.classList.add('active');
    }}

    function filterGeneSidebar(query) {{
        query = query.toLowerCase().trim();
        var items = document.querySelectorAll('.gene-item');
        items.forEach(function(item) {{
            var txt = (item.dataset.hugo + ' ' + item.dataset.desc).toLowerCase();
            item.style.display = (!query || txt.indexOf(query) >= 0) ? '' : 'none';
        }});
    }}

    window.onload = function() {{
        var first = document.querySelector('.gene-item');
        if (first) first.click();
    }};
</script>
</body>
</html>
"""

def adapt_to_genes(data: Dict[str, Any]) -> List[Dict[str, Any]]:
    genes_map = {}
    for r in data.get("records", []):
        hugo = r.get("hugo", "UNKNOWN")
        if hugo not in genes_map:
            gi = r.get("gene_info") or {}
            genes_map[hugo] = {
                "gene_symbol": hugo,
                "ncbi_id": gi.get("ncbi_gene_id"),
                "ncbi_summary": gi.get("description", "No summary available"),
                "map_location": gi.get("map_location"),
                "omim_id": r.get("omim_id"),
                "clinvar_disease": r.get("clinvar_disease"),
                "tier": r.get("tier", "Tier3"),
                "snps": []
            }
        
        curr_tier = genes_map[hugo]["tier"]
        new_tier = r.get("tier", "Tier3")
        if new_tier == "Tier1" or (new_tier == "Tier2" and curr_tier == "Tier3"):
            genes_map[hugo]["tier"] = new_tier

        rsid = r.get("rsid") or "-"
        chrom_pos = f"{r.get('chrom')}:{r.get('pos')} {r.get('ref')}>{r.get('alt')}"
        so_code = r.get("so")
        consequence = SO_NAME.get(so_code, so_code or "-")
        ev = r.get("evidence", {})
        zyg = ev.get("zygosity") or "-"
        phasing = ev.get("phasing") or ""
        vaf = ev.get("vaf")
        revel = r.get("revel")
        am = r.get("am_path")
        cadd = r.get("cadd_phred")
        clinvar = r.get("clinvar_sig") or "-"
        cid = r.get("clinvar_id")
        reasons = r.get("reason_codes", [])
        
        genes_map[hugo]["snps"].append({
            "rsid": rsid,
            "chrom_pos": chrom_pos,
            "consequence": consequence,
            "zygosity": zyg,
            "phasing": phasing,
            "vaf": vaf,
            "clinvar": clinvar,
            "clinvar_id": cid,
            "revel": revel,
            "am": am,
            "cadd": cadd,
            "reasons": reasons,
            "tier": new_tier
        })
        
    return sorted(genes_map.values(), key=lambda g: (g["tier"], g["gene_symbol"]))


def render_gene_browser(data: Dict[str, Any], output_path: str):
    genes = adapt_to_genes(data)
    sidebar_html = ""
    main_html = ""
    
    for g in genes:
        hugo = g["gene_symbol"]
        var_count = len(g["snps"])
        title = html.escape(g["ncbi_summary"])
        badge_cls = f"badge-{g['tier'].lower()}"
        
        sidebar_html += f"""
        <div id="nav-{hugo}" class="gene-item" data-hugo="{hugo}" data-desc="{title}" onclick="showGene('{hugo}')">
            <div class="gene-item-head">
                <span class="gene-item-name">{hugo}</span>
                <span class="gene-item-badge {badge_cls}">{g['tier']} ({var_count})</span>
            </div>
            <div class="gene-item-desc">{title}</div>
        </div>
        """
        
        snps_html = ""
        for s in g["snps"]:
            z_cls = "zyg-het" if s["zygosity"] == "Heterozygous" else ("zyg-hom" if s["zygosity"] == "Homozygous" else "zyg-hemi")
            
            clinvar_txt = html.escape(s['clinvar'])
            if s.get("clinvar_id"):
                clinvar_txt = f"<a href='https://www.ncbi.nlm.nih.gov/clinvar/variation/{s['clinvar_id']}/' target='_blank' style='color:#2563eb; font-weight:600;'>{clinvar_txt}</a>"
                
            rsid_txt = html.escape(str(s['rsid']))
            if str(s['rsid']).startswith("rs"):
                rsid_txt = f"<a href='https://www.ncbi.nlm.nih.gov/snp/{s['rsid']}' target='_blank' style='color:#2563eb; font-family:monospace;'>{rsid_txt}</a>"
                
            reasons_html = ""
            for rc in s["reasons"]:
                r_cls = "reason-tag"
                if rc in ("PVS1_HAPLOINSUFFICIENT", "PP3_CONSENSUS", "SPLICEAI_HIGH", "PM2_RARE"):
                    r_cls += " strong"
                elif rc.startswith("HPO_") or rc.startswith("GO_") or rc.startswith("CLIN"):
                    r_cls += " pheno"
                reasons_html += f"<span class='{r_cls}'>{html.escape(rc)}</span> "

            vaf_str = f" (VAF: {s['vaf']:.2f})" if s.get('vaf') is not None else ""
            ph_str = f"<br><small style='color:#7c3aed; font-weight:600;'>{html.escape(s['phasing'])}</small>" if s.get('phasing') and s.get('phasing') != 'Unphased (Short-Read WGS)' else ""

            snps_html += f"""
            <tr>
                <td><strong>{rsid_txt}</strong><br><span style="color:#64748b; font-family:monospace; font-size:0.75rem;">{html.escape(s['chrom_pos'])}</span></td>
                <td><span class="badge-so">{html.escape(s['consequence'])}</span></td>
                <td><span class="badge-zyg {z_cls}">{html.escape(s['zygosity'])}{vaf_str}</span>{ph_str}</td>
                <td>{clinvar_txt}</td>
                <td>
                    <div style="font-size:0.75rem; color:#475569;">
                        {f"REVEL: {s['revel']}<br>" if s['revel'] else ""}
                        {f"AM: {s['am']}<br>" if s['am'] else ""}
                        {f"CADD: {s['cadd']}" if s['cadd'] else ""}
                    </div>
                </td>
                <td><span class="gene-item-badge badge-{s['tier'].lower()}">{s['tier']}</span></td>
                <td style="max-width:320px;">{reasons_html}</td>
            </tr>
            """
            
        ncbi_link = f"<a href='https://www.ncbi.nlm.nih.gov/gene/{g['ncbi_id']}' target='_blank'>NCBI Gene ↗</a>" if g.get("ncbi_id") else ""
        hpo_link = f"<a href='https://hpo.jax.org/app/browse/search?q={hugo}' target='_blank'>HPO ↗</a>"
        gencc_link = f"<a href='https://search.thegencc.org/genes?q={hugo}' target='_blank'>GenCC ↗</a>"
        
        omim_box = ""
        if g.get("omim_id") or g.get("clinvar_disease"):
            omim_link = f"<a href='https://omim.org/entry/{g['omim_id']}' target='_blank' style='color:#7c3aed; font-weight:700;'>OMIM #{g['omim_id']} ↗</a>" if g.get("omim_id") else ""
            dis = html.escape(str(g.get('clinvar_disease') or ''))
            omim_box = f"<div class='omim-browser-box'><strong>OMIM Clinical Synopsis:</strong> {omim_link} {dis}</div>"

        loc_str = f" <span style='font-family:monospace; color:#64748b; font-size:0.85rem;'>[{html.escape(g['map_location'])}]</span>" if g.get("map_location") else ""

        main_html += f"""
        <div id="pane-{hugo}" class="gene-content-pane">
            <div class="gene-hero">
                <div class="gene-hero-title">
                    {hugo} {loc_str}
                    <span class="gene-item-badge badge-{g['tier'].lower()}">{g['tier']}</span>
                </div>
                <div class="gene-hero-summary">{title}</div>
                {omim_box}
                <div class="gene-links">
                    {ncbi_link}
                    {hpo_link}
                    {gencc_link}
                </div>
            </div>
            
            <h4 style="font-size:1rem; font-weight:700; color:#1e293b; margin-top:20px;">Actionable Variants ({var_count})</h4>
            <div class="table-wrap">
                <table class="variant-table">
                    <thead>
                        <tr>
                            <th>Variant</th>
                            <th>Consequence</th>
                            <th>Zygosity / Phase</th>
                            <th>ClinVar</th>
                            <th>Predictors</th>
                            <th>Tier</th>
                            <th>Reason Codes</th>
                        </tr>
                    </thead>
                    <tbody>
                        {snps_html}
                    </tbody>
                </table>
            </div>
        </div>
        """

    if not genes:
        main_html = "<div style='padding:24px; background:#f8fafc; border:1px solid #e2e8f0; border-radius:8px;'>No actionable variants found for this panel.</div>"

    final_html = HTML_TEMPLATE.format(
        sample_id=html.escape(data.get("patient", "UNKNOWN")),
        domain=html.escape(data.get("domain", "General")),
        total_genes=len(genes),
        total_variants=sum(len(g["snps"]) for g in genes),
        sidebar_html=sidebar_html,
        main_html=main_html
    )
    
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(final_html)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Render Two-Pane Variant Browser Layout")
    parser.add_argument("--in-json", required=True)
    parser.add_argument("--out-html", required=True)
    parser.add_argument("--out-tsv", required=True)
    parser.add_argument("--out-text", required=True)
    parser.add_argument("--split-reports", action="store_true", help="Ignored in browser mode")
    args = parser.parse_args()
    
    raw_data = json.load(open(args.in_json))
    render_gene_browser(raw_data, args.out_html)
    rr.write_tsv(raw_data.get("records", []), args.out_tsv)
    rr.write_text(raw_data, args.out_text)
    print(f"[render] Variant Browser HTML -> {os.path.abspath(args.out_html)}")
    print(f"[render] TSV -> {os.path.abspath(args.out_tsv)}")
    print(f"[render] text -> {os.path.abspath(args.out_text)}")
