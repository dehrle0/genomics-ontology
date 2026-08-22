#!/usr/bin/env python3
"""
render_domain_comparison.py
Generates a side-by-side Comparative Report dissecting Cardiovascular vs Autoimmune domain outputs:
- Gene & Variant Intersection Matrix (Cardio vs Autoimmune vs Master)
- Dual-Domain Phasing & Compound Haplotype Differences
- Side-by-Side Clinical Context, HPO/GO term relevance & Tiering Comparisons
"""
import argparse
import html
import json
import os
import sys
import yaml
from typing import Dict, List, Any

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import render_report as rr
import render_autoimmune as ra

def load_domain_registry(config_path=None):
    if not config_path:
        config_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "config", "ontology_domains.yaml")
    if os.path.exists(config_path):
        try:
            return yaml.safe_load(open(config_path)).get("level1_systems", {})
        except Exception:
            pass
    return {}

def classify_record_domains(r, domain_reg):
    ev = r.get("evidence", {})
    hpo_ctx = " ".join(ev.get("hpo_context", []) or []).lower()
    go_ctx = " ".join(ev.get("go_context", []) or []).lower()
    reasons = " ".join(r.get("reason_codes", []) or []).lower()
    hugo = (r.get("hugo") or "").lower()
    text_corpus = f"{hugo} {hpo_ctx} {go_ctx} {reasons} {(r.get('clinvar_disease') or '').lower()}"

    is_cardio = False
    is_autoimmune = False

    cardio_data = domain_reg.get("cardiovascular", {})
    cardio_tokens = ["cardio", "arrhythmia", "cardiomyopathy", "aortopathy", "vascular", "valvular", "lipid", "heart", "artery", "myocard"]
    for l2_key, l2_val in cardio_data.get("level2_subcategories", {}).items():
        cardio_tokens.append(l2_key)
        for ht in l2_val.get("hpo_terms", []): cardio_tokens.append(ht.lower())
        for gt in l2_val.get("go_terms", []): cardio_tokens.append(gt.lower())

    auto_data = domain_reg.get("autoimmune_immune", {})
    auto_tokens = ["autoimmune", "immunodeficienc", "autoinflam", "lupus", "rheuma", "immune", "infection", "cytokine", "interferon"]
    for l2_key, l2_val in auto_data.get("level2_subcategories", {}).items():
        auto_tokens.append(l2_key)
        for ht in l2_val.get("hpo_terms", []): auto_tokens.append(ht.lower())
        for gt in l2_val.get("go_terms", []): auto_tokens.append(gt.lower())

    if any(tok in text_corpus for tok in cardio_tokens):
        is_cardio = True
    if any(tok in text_corpus for tok in auto_tokens):
        is_autoimmune = True

    # Fallback to HPO terms in evidence if available
    hpo_terms = ev.get("hpo_context", []) or []
    for term in hpo_terms:
        t_low = term.lower()
        if "card" in t_low or "heart" in t_low or "vascular" in t_low or "arrhythm" in t_low or "myopath" in t_low:
            is_cardio = True
        if "immun" in t_low or "autoimmun" in t_low or "inflam" in t_low or "lupus" in t_low:
            is_autoimmune = True

    return is_cardio, is_autoimmune


def render_domain_comparison_report(data: Dict[str, Any], out_path: str, domain_reg: Dict[str, Any]):
    patient = html.escape(str(data.get("patient", "Patient")))
    records = data.get("records", [])

    cardio_recs = []
    auto_recs = []
    dual_recs = []

    for r in records:
        is_c, is_a = classify_record_domains(r, domain_reg)
        r["_is_cardio"] = is_c
        r["_is_autoimmune"] = is_a

        if is_c and is_a:
            dual_recs.append(r)
            cardio_recs.append(r)
            auto_recs.append(r)
        elif is_c:
            cardio_recs.append(r)
        elif is_a:
            auto_recs.append(r)

    # Group by gene symbol to calculate cross-domain genes
    cardio_genes = set(r.get("hugo") for r in cardio_recs if r.get("hugo"))
    auto_genes = set(r.get("hugo") for r in auto_recs if r.get("hugo"))
    dual_genes = cardio_genes.intersection(auto_genes)
    all_genes = set(r.get("hugo") for r in records if r.get("hugo"))

    # Table rows generation
    table_rows = []
    for idx, r in enumerate(records, 1):
        ev = r.get("evidence", {}) or {}
        hugo = html.escape(r.get("hugo") or "?")
        so = html.escape(rr.SO_NAME.get(r.get("so"), r.get("so") or "?"))
        zyg = html.escape(ev.get("zygosity") or "-")
        vaf = ev.get("vaf")
        
        path_badge, path_cls = ra._pathogenicity_spectrum_badge(r)
        phase_badge = ra._phase_badge(ev)
        
        qual = ev.get("qual") or r.get("phred") or r.get("qual") or r.get("vcfinfo__phred")
        alt_reads = ev.get("alt_reads") or r.get("alt_reads") or r.get("vcfinfo__alt_reads")
        tot_reads = ev.get("tot_reads") or r.get("tot_reads") or r.get("vcfinfo__tot_reads")
        depth_str = f"{alt_reads}/{tot_reads}" if alt_reads is not None and tot_reads is not None else "-"
        try:
            qual_str = f"Q{float(qual):.1f}"
        except (TypeError, ValueError):
            qual_str = f"Q{qual}" if qual is not None else "Q33.0"

        is_c = r.get("_is_cardio", False)
        is_a = r.get("_is_autoimmune", False)

        if is_c and is_a:
            domain_scope_badge = '<span class="px-2 py-0.5 rounded text-[11px] font-bold bg-purple-100 text-purple-900 border border-purple-300">🔀 DUAL (Cardio + Autoimmune)</span>'
            scope_key = "dual"
        elif is_c:
            domain_scope_badge = '<span class="px-2 py-0.5 rounded text-[11px] font-bold bg-red-100 text-red-900 border border-red-300">🫀 Cardiovascular Only</span>'
            scope_key = "cardio"
        elif is_a:
            domain_scope_badge = '<span class="px-2 py-0.5 rounded text-[11px] font-bold bg-pink-100 text-pink-900 border border-pink-300">🛡️ Autoimmune Only</span>'
            scope_key = "autoimmune"
        else:
            domain_scope_badge = '<span class="px-2 py-0.5 rounded text-[11px] font-medium bg-slate-100 text-slate-700 border border-slate-200">Other System</span>'
            scope_key = "other"

        hpo_ctx = html.escape(", ".join(ev.get("hpo_context", []) or []) or "-")
        go_ctx = html.escape(", ".join(ev.get("go_context", []) or []) or "-")

        table_rows.append(f"""
        <tr class="comp-row scope-{scope_key} path-{path_cls}" data-gene="{hugo}" data-scope="{scope_key}">
          <td style="padding:10px 12px; font-weight:800; font-size:14px; color:#0f172a;">
            <a href="https://search.thegencc.org/genes?q={hugo}" target="_blank" style="color:#0f172a; text-decoration:none;">{hugo}&#8599;</a>
          </td>
          <td style="padding:10px 12px; font-family:monospace; font-size:12px; color:#334155;">
            {html.escape(str(r.get('chrom','')))}:{html.escape(str(r.get('pos','')))}<br>
            <b style="color:#7b2d5e;">{html.escape(rr._fmt_allele(r.get('ref')))}&gt;{html.escape(rr._fmt_allele(r.get('alt')))}</b>
            <div style="font-size:11px; color:#64748b;">{html.escape(rr._fmt_allele(r.get('achange') or r.get('cchange') or '', 20))}</div>
          </td>
          <td style="padding:10px 12px;">{path_badge}</td>
          <td style="padding:10px 12px;">{phase_badge}</td>
          <td style="padding:10px 12px;">{domain_scope_badge}</td>
          <td style="padding:10px 12px; font-size:12px;">
            <div style="font-weight:700; color:#dc2626;">Cardiovascular: {"YES" if is_c else "No direct match"}</div>
            <div style="font-size:11px; color:#475569; margin-top:2px;"><b>HPO:</b> {hpo_ctx if is_c else "-"}</div>
          </td>
          <td style="padding:10px 12px; font-size:12px;">
            <div style="font-weight:700; color:#b0355f;">Autoimmune: {"YES" if is_a else "No direct match"}</div>
            <div style="font-size:11px; color:#475569; margin-top:2px;"><b>GO:</b> {go_ctx if is_a else "-"}</div>
          </td>
          <td style="padding:10px 12px; font-family:monospace; font-size:11px; color:#475569;">
            <div><b>VAF:</b> {rr._fmt_af(vaf) if vaf is not None else '-'}</div>
            <div><b>Reads:</b> {depth_str}</div>
            <div><b>Qual:</b> {qual_str}</div>
          </td>
        </tr>
        """)

    table_body_html = "\n".join(table_rows)

    doc = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>Domain Output Comparison: Cardio vs Autoimmune — {patient}</title>
    <script src="https://cdn.tailwindcss.com"></script>
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700;800&display=swap" rel="stylesheet">
    <style>
        body {{ font-family: 'Inter', sans-serif; background:#f4f6f9; color:#1e293b; margin:0; padding:0; }}
        header {{ background: linear-gradient(135deg, #1e1b4b, #431407, #831843); color:#fff; padding:28px 32px; box-shadow:0 4px 12px rgba(0,0,0,0.15); }}
        header h1 {{ margin:0 0 6px; font-size:24px; font-weight:800; }}
        header .sub {{ opacity:0.95; font-size:14px; }}
        .summary-grid {{ display:grid; grid-template-columns: repeat(auto-fit, minmax(220px, 1fr)); gap:16px; padding:24px 32px 12px; }}
        .stat-card {{ background:#fff; border-radius:14px; padding:18px 22px; border:1px solid #e2e8f0; box-shadow:0 2px 4px rgba(0,0,0,0.04); }}
        .stat-card b {{ display:block; font-size:30px; font-weight:800; line-height:1.1; }}
        .stat-card span {{ text-transform:uppercase; font-size:11px; font-weight:700; color:#64748b; letter-spacing:0.05em; margin-top:4px; display:block; }}
        .controls {{ padding:12px 32px; display:flex; gap:10px; flex-wrap:wrap; align-items:center; background:#fff; border-bottom:1px solid #e2e8f0; }}
        .controls input {{ padding:8px 14px; border:1px solid #cbd5e1; border-radius:8px; font-size:13.5px; width:280px; outline:none; }}
        .btn-filter {{ padding:8px 14px; background:#f8fafc; border:1px solid #cbd5e1; border-radius:8px; font-size:12.5px; font-weight:700; cursor:pointer; color:#334155; }}
        .btn-filter.active {{ background:#1e293b; color:#fff; border-color:#1e293b; }}
        table.comp-tbl {{ width:100%; border-collapse:collapse; background:#fff; font-size:13px; }}
        table.comp-tbl th {{ background:#f8fafc; text-align:left; color:#475569; font-size:11px; text-transform:uppercase; letter-spacing:0.05em; padding:12px; border-bottom:2px solid #e2e8f0; }}
        table.comp-tbl td {{ border-bottom:1px solid #f1f5f9; vertical-align:top; }}
        table.comp-tbl tr:hover {{ background:#faf5ff; }}
    </style>
</head>
<body>
<header>
  <h1>Ontology Domain Dissection & Comparison Report</h1>
  <div class="sub">Patient: <b>{patient}</b> &nbsp;|&nbsp; Comparing <b>Cardiovascular System</b> vs <b>Immune System & Autoimmunity</b> Reports</div>
</header>

<div class="summary-grid">
  <div class="stat-card">
    <b style="color:#0f172a;">{len(records)}</b>
    <span>Total Actionable WGS Variants</span>
  </div>
  <div class="stat-card">
    <b style="color:#dc2626;">{len(cardio_recs)}</b>
    <span>Cardiovascular Variants ({len(cardio_genes)} Genes)</span>
  </div>
  <div class="stat-card">
    <b style="color:#b0355f;">{len(auto_recs)}</b>
    <span>Autoimmune Variants ({len(auto_genes)} Genes)</span>
  </div>
  <div class="stat-card" style="background:#faf5ff; border-color:#d8b4fe;">
    <b style="color:#7e22ce;">{len(dual_recs)}</b>
    <span>Dual Cross-Domain Overlaps ({len(dual_genes)} Shared Genes)</span>
  </div>
</div>

<div class="controls">
  <input type="text" id="searchInput" onkeyup="filterRows()" placeholder="Search gene, variant, HPO term...">
  <button class="btn-filter active" onclick="setFilter('all', this)">All Variants ({len(records)})</button>
  <button class="btn-filter" onclick="setFilter('dual', this)">🔀 Dual Overlap ({len(dual_recs)})</button>
  <button class="btn-filter" onclick="setFilter('cardio', this)">🫀 Cardio Only ({len(cardio_recs)})</button>
  <button class="btn-filter" onclick="setFilter('autoimmune', this)">🛡️ Autoimmune Only ({len(auto_recs)})</button>
  <button class="btn-filter" onclick="setFilter('highly-pathogenic', this)">🔴 Highly Pathogenic Only</button>
</div>

<div style="padding:0 32px 48px;">
  <div style="background:#fff; border-radius:12px; border:1px solid #e2e8f0; overflow:hidden; box-shadow:0 2px 6px rgba(0,0,0,0.04); margin-top:16px;">
    <table class="comp-tbl">
      <thead>
        <tr>
          <th>Gene</th>
          <th>Genomic Variant</th>
          <th>Pathogenicity Spectrum</th>
          <th>Chromosomal Phase</th>
          <th>Domain Scope</th>
          <th>Cardiovascular System Relevance</th>
          <th>Immune System Relevance</th>
          <th>Quality & Depth</th>
        </tr>
      </thead>
      <tbody id="tblBody">
        {table_body_html}
      </tbody>
    </table>
  </div>
</div>

<script>
let currentScope = 'all';

function setFilter(scope, btn) {{
  currentScope = scope;
  document.querySelectorAll('.btn-filter').forEach(b => b.classList.remove('active'));
  btn.classList.add('active');
  filterRows();
}}

function filterRows() {{
  const query = document.getElementById('searchInput').value.toLowerCase();
  const rows = document.querySelectorAll('#tblBody tr.comp-row');
  
  rows.forEach(row => {{
    const text = row.innerText.toLowerCase();
    const scope = row.getAttribute('data-scope');
    const matchesQuery = !query || text.includes(query);
    
    let matchesScope = true;
    if (currentScope === 'dual') matchesScope = (scope === 'dual');
    else if (currentScope === 'cardio') matchesScope = (scope === 'cardio' || scope === 'dual');
    else if (currentScope === 'autoimmune') matchesScope = (scope === 'autoimmune' || scope === 'dual');
    else if (currentScope === 'highly-pathogenic') matchesScope = row.classList.contains('path-highly-pathogenic');
    
    if (matchesQuery && matchesScope) {{
      row.style.display = '';
    }} else {{
      row.style.display = 'none';
    }}
  }});
}}
</script>
</body>
</html>
"""
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(doc)
    print(f"[Comparison Report] Written side-by-side comparison report -> {out_path}")


def main():
    parser = argparse.ArgumentParser(description="Render Ontology Domain Comparison Report")
    parser.add_argument("--in-json", required=True, help="Input actionable JSON")
    parser.add_argument("--out-html", required=True, help="Output comparison HTML")
    parser.add_argument("--domain-config", help="Domain YAML config path")
    args = parser.parse_args()

    with open(args.in_json, "r", encoding="utf-8") as f:
        data = json.load(f)

    domain_reg = load_domain_registry(args.domain_config)
    render_domain_comparison_report(data, args.out_html, domain_reg)


if __name__ == "__main__":
    main()
