#!/usr/bin/env python3
"""
render_master_hub.py
Unified Master Genomics & Ontology Report Portal:
- Dynamic selection across ALL HPO Level 1 Organ Systems + Level 2 Subcategories
- Real-time client-side domain switching & search
- Bold NCBI Gene summaries & OMIM Clinical Synopsis
- Maternal vs Paternal phase tracking & chromosome phase blocks
- SVG Trait burden visualization & live GWAS evidence
- Monogenic & Polygenic integrated view
"""
import argparse
import html
import json
import math
import os
import sys
import yaml
from typing import Dict, List, Any

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import render_report as rr
import render_autoimmune as ra

SO_NAME = rr.SO_NAME
TIER_COLOR = rr.TIER_COLOR
TIER_LABEL = rr.TIER_LABEL


def load_domain_registry(config_path=None):
    if not config_path:
        config_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "config", "ontology_domains.yaml")
    if os.path.exists(config_path):
        try:
            return yaml.safe_load(open(config_path)).get("level1_systems", {})
        except Exception:
            pass
    return {}


def assign_domain_categories(records, domain_reg):
    for r in records:
        ev = r.get("evidence", {})
        hpo_ctx = " ".join(ev.get("hpo_context", []) or []).lower()
        go_ctx = " ".join(ev.get("go_context", []) or []).lower()
        reasons = " ".join(r.get("reason_codes", []) or []).lower()
        hugo = (r.get("hugo") or "").lower()
        text_corpus = f"{hugo} {hpo_ctx} {go_ctx} {reasons} {(r.get('clinvar_disease') or '').lower()}"

        matched_l1 = set()
        matched_l2 = set()

        for l1_key, l1_data in domain_reg.items():
            l1_title = l1_data.get("title", l1_key)
            l2_dict = l1_data.get("level2_subcategories", {})
            
            for l2_key, l2_data in l2_dict.items():
                l2_title = l2_data.get("title", l2_key)
                tokens = [l2_key, l1_key] + [t.lower() for t in l2_title.replace("(", " ").replace(")", " ").replace(",", " ").split() if len(t) > 3]
                if any(tok in text_corpus for tok in tokens):
                    matched_l1.add(l1_key)
                    matched_l2.add(f"{l1_key}:{l2_key}")

        if not matched_l1:
            matched_l1.add("cardiovascular" if "cardio" in text_corpus or "myh" in hugo else ("autoimmune_immune" if "immune" in text_corpus or "ptpn" in hugo else "other"))
            matched_l2.add(f"{list(matched_l1)[0]}:general")

        r["matched_level1"] = sorted(list(matched_l1))
        r["matched_level2"] = sorted(list(matched_l2))


def render_master_hub(data, output_path, domain_reg=None):
    if not domain_reg:
        domain_reg = load_domain_registry()

    records = data.get("records", [])
    assign_domain_categories(records, domain_reg)
    
    patient = data.get("patient", "Patient")
    title = data.get("report_title", "Universal Genomic & Ontology Master Report")
    tc = data.get("tier_counts", {"Tier1": 0, "Tier2": 0, "Tier3": 0})
    
    trait_rows = ra.collect_traits(records)
    chart_svg = ra.trait_chart_svg(trait_rows, top_n=12)

    import collections
    gene_groups = collections.OrderedDict()
    for r in records:
        hugo = r.get("hugo", "UNKNOWN")
        gene_groups.setdefault(hugo, []).append(r)

    l1_buttons_html = []
    l1_buttons_html.append(f'<button class="domain-pill active" data-l1="all" onclick="selectLevel1(\'all\')">All Systems <span class="count">{len(gene_groups)}</span></button>')
    
    for l1_key, l1_data in domain_reg.items():
        count = sum(1 for g_vars in gene_groups.values() if any(l1_key in r.get("matched_level1", []) for r in g_vars))
        icon = l1_data.get("icon", "")
        title_s = l1_data.get("title", l1_key)
        l1_buttons_html.append(f'<button class="domain-pill" data-l1="{l1_key}" onclick="selectLevel1(\'{l1_key}\')">{icon} {html.escape(title_s)} <span class="count">{count}</span></button>')

    l2_options_html = []
    l2_options_html.append('<option value="all">All Level 2 Subcategories</option>')
    for l1_key, l1_data in domain_reg.items():
        l1_title = l1_data.get("title", l1_key)
        for l2_key, l2_data in l1_data.get("level2_subcategories", {}).items():
            l2_title = l2_data.get("title", l2_key)
            val = f"{l1_key}:{l2_key}"
            l2_options_html.append(f'<option value="{val}" data-l1="{l1_key}">[{html.escape(l1_title)}] {html.escape(l2_title)}</option>')

    l1_btns_str = "\n".join(l1_buttons_html)
    l2_opts_str = "\n".join(l2_options_html)

    cards_html = []
    tier_order = {"Tier1": 0, "Tier2": 1, "Tier3": 2}
    for hugo, g_vars in gene_groups.items():
        best_tier = min((r.get("tier", "Tier3") for r in g_vars), key=lambda t: tier_order.get(t, 9))
        l1_set = set()
        l2_set = set()
        for r in g_vars:
            l1_set.update(r.get("matched_level1", []))
            l2_set.update(r.get("matched_level2", []))
        l1_classes = " ".join(sorted(l1_set))
        l2_classes = " ".join(sorted(l2_set))
        
        gene_card_str = ra._gene_card(hugo, g_vars)
        wrapped = f'<div class="hub-card-wrapper" data-tier="{best_tier}" data-gene="{hugo}" data-l1="{l1_classes}" data-l2="{l2_classes}">{gene_card_str}</div>'
        cards_html.append(wrapped)

    body_cards = "\n".join(cards_html)

    doc = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>{html.escape(title)} — {html.escape(patient)}</title>
    <script src="https://cdn.tailwindcss.com"></script>
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700;800&display=swap" rel="stylesheet">
    <style>
        body {{
            font-family: 'Inter', sans-serif;
            background: radial-gradient(circle at 10% 20%, #f8fafc, #edf2f7);
            color: #0f172a;
            min-height: 100vh;
        }}
        .domain-pill {{
            display: inline-flex;
            align-items: center;
            gap: 6px;
            padding: 8px 16px;
            background: #ffffff;
            border: 1px solid #cbd5e1;
            border-radius: 12px;
            font-size: 13px;
            font-weight: 600;
            color: #334155;
            cursor: pointer;
            transition: all 0.2s ease;
            box-shadow: 0 1px 2px rgba(0,0,0,0.04);
        }}
        .domain-pill:hover {{
            background: #f1f5f9;
            border-color: #94a3b8;
            transform: translateY(-1px);
        }}
        .domain-pill.active {{
            background: #1e293b;
            color: #ffffff;
            border-color: #0f172a;
            box-shadow: 0 4px 12px rgba(15, 23, 42, 0.15);
        }}
        .domain-pill .count {{
            background: rgba(0,0,0,0.06);
            padding: 1px 7px;
            border-radius: 8px;
            font-size: 11px;
        }}
        .domain-pill.active .count {{
            background: rgba(255,255,255,0.2);
            color: #ffffff;
        }}
        .card {{
            background: #ffffff;
            border-radius: 14px;
            padding: 20px 24px;
            margin: 16px 0;
            box-shadow: 0 2px 8px rgba(0,0,0,0.04);
            border: 1px solid #e2e8f0;
            transition: transform 0.15s ease, box-shadow 0.15s ease;
        }}
        .card:hover {{
            box-shadow: 0 6px 20px rgba(0,0,0,0.06);
            transform: translateY(-1px);
        }}
        .card-head {{
            display: flex;
            gap: 12px;
            align-items: baseline;
            flex-wrap: wrap;
            border-bottom: 1px solid #f1f5f9;
            padding-bottom: 12px;
        }}
        .card-head .gene {{
            font-weight: 800;
            font-size: 20px;
            color: #0f172a;
        }}
        .card-head .loc {{
            font-family: ui-monospace, monospace;
            color: #64748b;
            font-size: 13.5px;
        }}
        .card-head .so {{
            background: #f1f5f9;
            border-radius: 6px;
            padding: 2px 8px;
            font-size: 12px;
            font-weight: 600;
            color: #475569;
        }}
        .card-head .ach {{
            color: #4338ca;
            font-size: 14px;
            font-family: ui-monospace, monospace;
            font-weight: 600;
        }}
        .gene-desc-bold {{
            font-size: 14.5px;
            color: #0f172a;
            margin: 12px 0 6px;
            line-height: 1.55;
            font-weight: 600;
        }}
        .gene-desc-bold strong {{
            color: #0f172a;
            font-weight: 700;
        }}
        .gd-label {{
            font-size: 10px;
            text-transform: uppercase;
            letter-spacing: .05em;
            font-weight: 700;
            color: #fff;
            background: #2563eb;
            border-radius: 4px;
            padding: 2px 7px;
            margin-right: 8px;
            display: inline-block;
        }}
        .omim-block {{
            font-size: 13px;
            color: #475569;
            margin: 4px 0 8px;
            padding: 8px 14px;
            background: #faf5ff;
            border-left: 4px solid #7c3aed;
            border-radius: 6px;
        }}
        .omim-label {{
            font-size: 10.5px;
            text-transform: uppercase;
            letter-spacing: .05em;
            font-weight: 700;
            color: #7c3aed;
            margin-right: 8px;
        }}
        .grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(140px, 1fr));
            gap: 10px 18px;
            margin: 14px 0;
            background: #f8fafc;
            padding: 12px 16px;
            border-radius: 10px;
            border: 1px solid #f1f5f9;
        }}
        .grid label {{
            display: block;
            font-size: 10.5px;
            color: #64748b;
            text-transform: uppercase;
            letter-spacing: .04em;
            font-weight: 700;
            margin-bottom: 2px;
        }}
        .grid div {{
            font-size: 13.5px;
            color: #1e293b;
            font-weight: 500;
        }}
        .onto-box {{
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 14px;
            margin: 14px 0;
            padding: 14px 16px;
            background: #f8fafc;
            border: 1px solid #e2e8f0;
            border-radius: 10px;
        }}
        .onto-item label {{
            display: block;
            font-size: 10.5px;
            font-weight: 700;
            text-transform: uppercase;
            letter-spacing: .04em;
            color: #64748b;
            margin-bottom: 3px;
        }}
        .onto-text {{
            font-size: 13.5px;
            color: #1e293b;
            font-weight: 500;
        }}
        .onto-link {{
            font-size: 12px;
            color: #2563eb;
            font-weight: 700;
            text-decoration: none;
            margin-left: 4px;
        }}
        .onto-link:hover {{
            text-decoration: underline;
        }}
        .studies {{
            margin: 14px 0 10px;
        }}
        .studies-h {{
            font-size: 11.5px;
            text-transform: uppercase;
            letter-spacing: .05em;
            color: #7b2d5e;
            font-weight: 700;
            margin-bottom: 6px;
        }}
        .study-tbl {{
            width: 100%;
            border-collapse: collapse;
            font-size: 12.5px;
        }}
        .study-tbl th {{
            text-align: left;
            color: #64748b;
            font-weight: 700;
            border-bottom: 1px solid #e2e8f0;
            padding: 6px 8px;
            font-size: 11px;
            text-transform: uppercase;
        }}
        .study-tbl td {{
            padding: 6px 8px;
            border-bottom: 1px solid #f1f5f9;
            vertical-align: top;
        }}
        .study-tbl td.trait {{
            max-width: 320px;
            font-weight: 500;
        }}
        .study-tbl td.pv {{
            font-family: ui-monospace, monospace;
            color: #dc2626;
            font-weight: 700;
        }}
        .study-tbl td.ra {{
            font-family: ui-monospace, monospace;
        }}
        .reasons-wrap {{
            display: flex;
            flex-wrap: wrap;
            gap: 6px;
            margin-top: 12px;
        }}
        @media print {{
            .noprint {{ display: none !important; }}
            .card {{ page-break-inside: avoid; box-shadow: none !important; }}
        }}
    </style>
</head>
<body class="antialiased">

    <!-- Top Master Header -->
    <header class="relative overflow-hidden bg-slate-900 text-white px-8 py-10 shadow-xl">
        <div class="absolute top-0 right-0 -mr-20 -mt-20 w-96 h-96 rounded-full bg-indigo-500 blur-3xl opacity-20 pointer-events-none"></div>
        <div class="absolute bottom-0 left-0 -ml-20 -mb-20 w-80 h-80 rounded-full bg-blue-500 blur-3xl opacity-20 pointer-events-none"></div>
        
        <div class="relative z-10 max-w-7xl mx-auto flex justify-between items-center flex-wrap gap-4">
            <div>
                <h1 class="text-3xl font-extrabold tracking-tight mb-2">{html.escape(title)}</h1>
                <div class="text-slate-300 font-medium text-sm flex items-center flex-wrap gap-x-6 gap-y-2">
                    <span>Patient: <b class="text-white">{html.escape(patient)}</b></span>
                    <span class="text-slate-600">|</span>
                    <span>All Ontology Level 1 + Level 2 Combinations</span>
                    <span class="text-slate-600">|</span>
                    <span>Screened: <b class="text-white">{data.get('scanned_panel_variants', len(records))}</b> variants</span>
                </div>
            </div>
            <div class="noprint">
                <button onclick="window.print()" class="px-5 py-2.5 bg-white/10 hover:bg-white/20 border border-white/20 text-white rounded-xl shadow-sm font-semibold text-sm flex items-center gap-2 transition-all">
                    <svg class="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M17 17h2a2 2 0 002-2v-4a2 2 0 00-2-2H5a2 2 0 00-2 2v4a2 2 0 002 2h2m2 4h6a2 2 0 002-2v-4a2 2 0 00-2-2H9a2 2 0 00-2 2v4a2 2 0 002 2zm8-12V5a2 2 0 00-2-2H9a2 2 0 00-2 2v4h10z"></path></svg>
                    Print / Export PDF
                </button>
            </div>
        </div>
    </header>

    <div class="max-w-7xl mx-auto px-6 py-8">
        
        <!-- Summary Stats Cards -->
        <div class="grid grid-cols-2 md:grid-cols-5 gap-4 mb-8">
            <div class="bg-white border border-slate-200 shadow-sm rounded-2xl p-5 flex flex-col justify-center">
                <span class="text-3xl font-extrabold text-slate-800">{len(records)}</span>
                <span class="text-[11px] uppercase tracking-wider font-bold text-slate-500 mt-1">Actionable Variants</span>
            </div>
            <div class="bg-white border border-slate-200 shadow-sm rounded-2xl p-5 flex flex-col justify-center">
                <span class="text-3xl font-extrabold text-red-600">{tc.get('Tier1', 0)}</span>
                <span class="text-[11px] uppercase tracking-wider font-bold text-slate-500 mt-1">Tier 1 Pathogenic</span>
            </div>
            <div class="bg-white border border-slate-200 shadow-sm rounded-2xl p-5 flex flex-col justify-center">
                <span class="text-3xl font-extrabold text-amber-500">{tc.get('Tier2', 0)}</span>
                <span class="text-[11px] uppercase tracking-wider font-bold text-slate-500 mt-1">Tier 2 VUS</span>
            </div>
            <div class="bg-white border border-slate-200 shadow-sm rounded-2xl p-5 flex flex-col justify-center">
                <span class="text-3xl font-extrabold text-blue-600">{tc.get('Tier3', 0)}</span>
                <span class="text-[11px] uppercase tracking-wider font-bold text-slate-500 mt-1">Tier 3 Risk Alleles</span>
            </div>
            <div class="bg-white border border-slate-200 shadow-sm rounded-2xl p-5 flex flex-col justify-center">
                <span class="text-3xl font-extrabold text-purple-600">{len(set(r.get('hugo') for r in records))}</span>
                <span class="text-[11px] uppercase tracking-wider font-bold text-slate-500 mt-1">Flagged Genes</span>
            </div>
        </div>

        <!-- Trait Burden Overview Chart -->
        <div class="bg-white border border-slate-200 rounded-2xl p-6 shadow-sm mb-8">
            <h2 class="text-lg font-bold text-slate-900 mb-1">Global Trait Burden & Polygenic Evidence</h2>
            <p class="text-xs text-slate-500 mb-4">Catalogued GWAS risk-allele associations across your sequenced genome, coloured by association strength (-log10 p).</p>
            {chart_svg}
        </div>

        <!-- Master Interactive Domain Selector (Level 1 + Level 2) -->
        <div class="bg-white border border-slate-200 rounded-2xl p-4 md:p-6 shadow-sm mb-8 noprint relative z-10">
            <div class="mb-4">
                <label class="text-xs uppercase font-bold text-slate-500 tracking-wider block mb-2">Select Ontology Level 1 System:</label>
                <div class="flex flex-wrap gap-2" id="level1Container">
                    {l1_btns_str}
                </div>
            </div>

            <div class="grid grid-cols-1 md:grid-cols-2 gap-4 mt-4 pt-4 border-t border-slate-100">
                <div>
                    <label class="text-xs uppercase font-bold text-slate-500 tracking-wider block mb-1">Filter by Level 2 Subcategory:</label>
                    <select id="level2Select" class="w-full p-2.5 bg-slate-50 border border-slate-200 rounded-xl text-sm font-semibold text-slate-700 outline-none focus:ring-2 focus:ring-indigo-500" onchange="filterMasterHub()">
                        {l2_opts_str}
                    </select>
                </div>
                <div>
                    <label class="text-xs uppercase font-bold text-slate-500 tracking-wider block mb-1">Live Keyword / Gene / RSID Search:</label>
                    <input id="hubSearch" type="text" placeholder="Search gene (e.g. MYH7, PTPN22), RSID, or phenotype..." 
                        class="w-full p-2.5 bg-slate-50 border border-slate-200 rounded-xl text-sm font-medium text-slate-800 outline-none focus:ring-2 focus:ring-indigo-500"
                        oninput="filterMasterHub()">
                </div>
            </div>
        </div>

        <!-- Cards Container -->
        <div id="hubCardsContainer">
            {body_cards or '<div class="p-8 text-center text-slate-500 bg-white rounded-2xl border border-slate-200">No actionable variants found for the selected ontology filters.</div>'}
        </div>

    </div>

    <footer class="bg-slate-50 border-t border-slate-200 py-10 px-8 text-center text-slate-500 text-xs leading-relaxed mt-12">
        <div class="max-w-4xl mx-auto">
            Generated by <code class="bg-slate-100 text-slate-700 px-1 py-0.5 rounded font-mono">ontology_report</code> Master Hub.
            Supports dynamic multi-domain mapping across all HPO Level 1 Organ Systems and Level 2 Subcategories.
            <br><span class="font-bold text-slate-600 mt-2 block">Personal Screening & Exploratory Genomics — not an in-vitro clinical diagnostic device.</span>
        </div>
    </footer>

    <script>
    var currentL1 = 'all';

    function selectLevel1(l1Key) {{
        currentL1 = l1Key;
        document.querySelectorAll('#level1Container .domain-pill').forEach(function(btn) {{
            if (btn.dataset.l1 === l1Key) {{
                btn.classList.add('active');
            }} else {{
                btn.classList.remove('active');
            }}
        }});
        
        var l2Sel = document.getElementById('level2Select');
        for (var i = 0; i < l2Sel.options.length; i++) {{
            var opt = l2Sel.options[i];
            if (l1Key === 'all' || !opt.dataset.l1 || opt.dataset.l1 === l1Key) {{
                opt.style.display = '';
            }} else {{
                opt.style.display = 'none';
            }}
        }}
        l2Sel.value = 'all';
        filterMasterHub();
    }}

    function filterMasterHub() {{
        var q = document.getElementById('hubSearch').value.trim().toLowerCase();
        var l2Val = document.getElementById('level2Select').value;

        var wrappers = document.querySelectorAll('.hub-card-wrapper');
        wrappers.forEach(function(w) {{
            var l1Matches = (currentL1 === 'all') || (w.dataset.l1.indexOf(currentL1) >= 0);
            var l2Matches = (l2Val === 'all') || (w.dataset.l2.indexOf(l2Val) >= 0);
            var txtMatches = (!q) || (w.innerText.toLowerCase().indexOf(q) >= 0);

            if (l1Matches && l2Matches && txtMatches) {{
                w.style.display = '';
            }} else {{
                w.style.display = 'none';
            }}
        }});
    }}
    </script>
</body>
</html>"""
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(doc)
    print(f"[master-hub] Master Portal HTML -> {output_path}")


def main():
    ap = argparse.ArgumentParser(description="Render Unified Master Genomics & Ontology Hub")
    ap.add_argument("--in-json", required=True)
    ap.add_argument("--out-html", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-text", required=True)
    ap.add_argument("--domain-config", default=None)
    args = ap.parse_args()

    data = json.load(open(args.in_json))
    domain_reg = load_domain_registry(args.domain_config)
    render_master_hub(data, args.out_html, domain_reg)
    rr.write_tsv(data.get("records", []), args.out_tsv)
    rr.write_text(data, args.out_text)


if __name__ == "__main__":
    main()
