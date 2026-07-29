#!/usr/bin/env python3
"""
render_report.py
Turn the actionable-variant JSON (from ontology_filter.py) into human
deliverables: a styled interactive HTML report, a flat TSV, and a text summary.

This version features a modern glassmorphism aesthetic, sleek UI layout, 
and supports a wider array of clinical/functional annotators.
"""
import argparse
import html
import json
import os

TIER_LABEL = {
    "Tier1": "Tier 1 — Reportable / Pathogenic-grade",
    "Tier2": "Tier 2 — VUS of interest",
    "Tier3": "Tier 3 — Monitor / regulatory",
}
TIER_COLOR = {"Tier1": "#ef4444", "Tier2": "#f59e0b", "Tier3": "#3b82f6"}

SO_NAME = {
    "MIS": "missense", "STG": "stop-gained", "FSD": "frameshift-del",
    "FSI": "frameshift-ins", "SPL": "splice-site", "MLO": "start-lost",
    "STL": "stop-lost", "EXL": "exon-loss", "TAB": "transcript-ablation",
    "IND": "inframe-del", "INI": "inframe-ins", "CSS": "complex-sub",
    "SYN": "synonymous", "INT": "intronic", "UT3": "3'UTR", "UT5": "5'UTR",
    "2KU": "upstream", "2KD": "downstream", "NMD": "NMD-transcript",
}


def _fmt_allele(v, maxlen=14):
    if v is None or v == "":
        return "-"
    s = str(v)
    if len(s) <= maxlen:
        return s
    return f"{s[:maxlen]}…({len(s)}bp)"


def _fmt_af(v):
    if v is None or v == "":
        return "absent"
    try:
        f = float(v)
    except (TypeError, ValueError):
        return str(v)
    if f == 0:
        return "0"
    if f < 1e-4:
        return f"{f:.2e}"
    return f"{f:.4f}"


def _clinvar_link(cid, sig):
    if cid:
        url = f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{cid}/"
        return f'<a href="{url}" target="_blank">{html.escape(sig or "see ClinVar")}</a>'
    return html.escape(sig or "-")


def write_tsv(records, path):
    cols = ["tier", "hugo", "gene_description", "zygosity", "vaf", "rsid",
            "chrom", "pos", "ref", "alt", "so", "achange",
            "cchange", "transcript", "gnomad4_af", "allofus_af", "clinvar_sig",
            "clinvar_id", "revel", "am_path", "cadd_phred", "spliceai_max",
            "pharmgkb__chemicals", "civic__clinical_significance", "interpro__domain",
            "reason_codes"]
    with open(path, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in records:
            ev = r.get("evidence", {})
            row = []
            for c in cols:
                if c == "spliceai_max":
                    v = ev.get("spliceai_max")
                elif c == "zygosity":
                    v = ev.get("zygosity")
                elif c == "vaf":
                    v = ev.get("vaf")
                elif c == "gene_description":
                    v = (r.get("gene_info") or {}).get("description")
                elif c == "reason_codes":
                    v = "|".join(r.get("reason_codes", []))
                else:
                    v = r.get(c)
                row.append("" if v is None else str(v).replace("\t", " "))
            f.write("\t".join(row) + "\n")


def write_text(data, path):
    title = data.get("report_title", "Ontology-Driven Actionable Report")
    lines = []
    lines.append("=" * 70)
    lines.append(title.upper())
    lines.append(f"Patient: {data['patient']}    Domain: {data['domain']}")
    lines.append("=" * 70)
    lines.append(f"Ontology gene panel size : {data['panel_gene_count']}")
    lines.append(f"Panel-gene variants scan : {data['scanned_panel_variants']}")
    lines.append(f"Actionable variants kept : {data['actionable_count']}")
    tc = data["tier_counts"]
    lines.append(f"  Tier 1: {tc['Tier1']}   Tier 2: {tc['Tier2']}   Tier 3: {tc['Tier3']}")
    lines.append("")
    for tier in ("Tier1", "Tier2", "Tier3"):
        recs = [r for r in data["records"] if r["tier"] == tier]
        if not recs:
            continue
        lines.append("-" * 70)
        lines.append(f"{TIER_LABEL[tier]}  ({len(recs)})")
        lines.append("-" * 70)
        for r in recs:
            so = SO_NAME.get(r.get("so"), r.get("so") or "?")
            ev = r.get("evidence", {})
            zyg = ev.get("zygosity")
            lines.append(
                f"  {r['hugo']:10} {r.get('chrom','')}:{r.get('pos','')} "
                f"{_fmt_allele(r.get('ref'))}>{_fmt_allele(r.get('alt'))}  [{so}] "
                f"{_fmt_allele(r.get('achange'), 24) if r.get('achange') else ''}"
                f"{('  ' + zyg) if zyg else ''}"
            )
            gi = r.get("gene_info") or {}
            if gi.get("description"):
                lines.append(f"      NCBI: {gi['description']}"
                             f"{('  [' + gi['map_location'] + ']') if gi.get('map_location') else ''}")
            lines.append(f"      gnomAD4={_fmt_af(r.get('gnomad4_af'))} "
                         f"AoU={_fmt_af(r.get('allofus_af'))} "
                         f"ClinVar={r.get('clinvar_sig') or '-'}")
            lines.append(f"      reasons: {', '.join(r.get('reason_codes', []))}")
        lines.append("")
    with open(path, "w") as f:
        f.write("\n".join(lines))


def _reason_badges(reasons):
    out = []
    for rc in reasons:
        cls = "bg-blue-100 text-blue-800 border-blue-200"
        if rc.startswith("HPO_") or rc.startswith("GO_") or rc.startswith("CLIN") or rc in ("OMIM_DISEASE", "ARRVARS_KNOWN", "PHARMGKB_DRUG", "DENOVO_EVIDENCE"):
            cls = "bg-green-100 text-green-800 border-green-200"
        if rc in ("PVS1_HAPLOINSUFFICIENT", "PP3_CONSENSUS", "SPLICEAI_HIGH", "PM2_RARE"):
            cls = "bg-red-100 text-red-800 border-red-200 font-semibold"
        if rc == "COMMON_AF_FLAG":
            cls = "bg-orange-100 text-orange-800 border-orange-200"
        out.append(f'<span class="px-2 py-0.5 rounded-full text-[10px] border {cls}">{html.escape(rc)}</span>')
    return "".join(out)


def _zyg_badge(zyg):
    if not zyg:
        return ""
    cls = {"Heterozygous": "bg-yellow-100 text-yellow-800 border-yellow-200", 
           "Homozygous": "bg-red-100 text-red-800 border-red-200",
           "Hemizygous": "bg-purple-100 text-purple-800 border-purple-200"}.get(zyg, "bg-gray-100 text-gray-800 border-gray-200")
    return f'<span class="px-2 py-0.5 rounded-md text-xs font-semibold border {cls}">{html.escape(zyg)}</span>'


def _gene_desc_block(r):
    gi = r.get("gene_info")
    if not gi:
        return ""
    gid = gi.get("ncbi_gene_id")
    desc = gi.get("description") or ""
    loc = gi.get("map_location")
    link = (f'<a href="https://www.ncbi.nlm.nih.gov/gene/{html.escape(str(gid))}" '
            f'target="_blank" class="text-indigo-500 hover:text-indigo-600 transition-colors">NCBI&#8599;</a>') if gid else ""
    loc_txt = f' <span class="font-mono text-gray-500">{html.escape(loc)}</span>' if loc else ""
    if not (desc or gid):
        return ""
    return (f'<div class="text-sm text-gray-700 mt-2 mb-1 leading-relaxed"><span class="text-[10px] uppercase tracking-wider text-white bg-indigo-600 rounded px-1.5 py-0.5 mr-1.5">NCBI Gene</span> '
            f'{html.escape(desc)}{loc_txt} {link}</div>')


def _card(r):
    ev = r.get("evidence", {})
    so = SO_NAME.get(r.get("so"), r.get("so") or "?")
    hpo_ctx = ", ".join(ev.get("hpo_context", []) or []) or "-"
    go_ctx = ", ".join(ev.get("go_context", []) or []) or "-"
    gene = html.escape(r.get("hugo") or "?")
    gene_link = f'https://search.thegencc.org/genes?q={gene}'
    hpo_gene_link = f'https://hpo.jax.org/app/browse/search?q={gene}&navFilter=all'
    zyg = ev.get("zygosity")
    vaf = ev.get("vaf")
    rsid = r.get("rsid")
    rsid_html = (f'<a href="https://www.ncbi.nlm.nih.gov/snp/{html.escape(str(rsid))}" '
                 f'target="_blank" class="text-blue-500 hover:underline">{html.escape(str(rsid))}</a>') if rsid and str(rsid).startswith("rs") else (html.escape(str(rsid)) if rsid else "-")
    
    # Extra annotations
    pharm = r.get("pharmgkb__chemicals")
    pharm_html = f'<div><label>PharmGKB</label><span class="truncate block" title="{html.escape(str(pharm))}">{html.escape(str(pharm))}</span></div>' if pharm else ""
    civic = r.get("civic__clinical_significance")
    civic_html = f'<div><label>CIViC</label><span class="truncate block" title="{html.escape(str(civic))}">{html.escape(str(civic))}</span></div>' if civic else ""
    interpro = r.get("interpro__domain")
    interpro_html = f'<div><label>InterPro Domain</label><span class="truncate block text-xs" title="{html.escape(str(interpro))}">{html.escape(str(interpro))}</span></div>' if interpro else ""
    denovo = r.get("denovo__PubmedID")
    denovo_html = f'<div><label>DeNovo PMID</label><span>{html.escape(str(denovo))}</span></div>' if denovo else ""
    
    return f"""
    <div class="variant-card bg-white/80 backdrop-blur-xl border border-white/40 shadow-[0_8px_30px_rgb(0,0,0,0.04)] rounded-2xl p-5 mb-5 transition-all duration-300 hover:shadow-[0_8px_30px_rgb(0,0,0,0.08)] hover:-translate-y-0.5" data-gene="{gene}" data-reasons="{html.escape(' '.join(r.get('reason_codes', [])))}">
      <div class="flex flex-wrap items-baseline gap-3 pb-3 border-b border-slate-100">
        <span class="font-bold text-lg text-slate-800">{gene}</span>
        {_zyg_badge(zyg)}
        <span class="font-mono text-slate-500 text-sm">{html.escape(str(r.get('chrom','')))}:{html.escape(str(r.get('pos','')))}
          {html.escape(_fmt_allele(r.get('ref')))}&gt;{html.escape(_fmt_allele(r.get('alt')))}</span>
        <span class="bg-slate-100 text-slate-600 rounded-lg px-2.5 py-0.5 text-xs font-medium">{html.escape(so)}</span>
        <span class="text-indigo-700 font-mono text-sm">{html.escape(_fmt_allele(r.get('achange') or r.get('cchange') or '', 28))}</span>
      </div>
      {_gene_desc_block(r)}
      <div class="grid grid-cols-2 md:grid-cols-4 lg:grid-cols-5 xl:grid-cols-6 gap-x-4 gap-y-3 mt-4 mb-4 metrics-grid">
        <div><label>Zygosity</label><span>{html.escape(zyg or '-')}</span></div>
        <div><label>Variant allele frac</label><span>{_fmt_af(vaf) if vaf is not None else '-'}</span></div>
        <div><label>dbSNP</label><span>{rsid_html}</span></div>
        <div><label>gnomAD4 AF</label><span>{_fmt_af(r.get('gnomad4_af'))}</span></div>
        <div><label>All of Us AF</label><span>{_fmt_af(r.get('allofus_af'))}</span></div>
        <div><label>ClinVar</label><span>{_clinvar_link(r.get('clinvar_id'), r.get('clinvar_sig'))}</span></div>
        <div><label>REVEL</label><span>{html.escape(str(r.get('revel') or '-'))}</span></div>
        <div><label>AlphaMissense</label><span>{html.escape(str(r.get('am_path') or '-'))}</span></div>
        <div><label>SpliceAI max</label><span>{html.escape(str(ev.get('spliceai_max') if ev.get('spliceai_max') is not None else '-'))}</span></div>
        <div><label>CADD phred</label><span>{html.escape(str(r.get('cadd_phred') or '-'))}</span></div>
        <div><label>Panel support</label><span>{html.escape(str(ev.get('panel_support') or '-'))}/2</span></div>
        {pharm_html}
        {civic_html}
        {interpro_html}
        {denovo_html}
      </div>
      <div class="grid grid-cols-1 md:grid-cols-2 gap-x-4 gap-y-2 mt-3 mb-4 text-sm bg-slate-50/50 p-3 rounded-xl border border-slate-100/50">
        <div><label class="text-[10px] uppercase tracking-wider text-slate-500 font-semibold block mb-0.5">HPO phenotype context</label><span class="text-slate-700">{html.escape(hpo_ctx)}</span>
          &nbsp;<a href="{hpo_gene_link}" target="_blank" class="text-indigo-500 hover:underline text-xs">HPO&#8599;</a></div>
        <div><label class="text-[10px] uppercase tracking-wider text-slate-500 font-semibold block mb-0.5">GO function context</label><span class="text-slate-700">{html.escape(go_ctx)}</span>
          &nbsp;<a href="{gene_link}" target="_blank" class="text-indigo-500 hover:underline text-xs">GenCC&#8599;</a></div>
      </div>
      <div class="flex flex-wrap gap-1.5 mt-2">{_reason_badges(r.get('reason_codes', []))}</div>
    </div>"""


def write_html(data, path):
    tc = data["tier_counts"]
    title = data.get("report_title", "Ontology-Driven Actionable Report")
    sections = []
    for tier in ("Tier1", "Tier2", "Tier3"):
        recs = [r for r in data["records"] if r["tier"] == tier]
        if not recs:
            continue
        cards = "\n".join(_card(r) for r in recs)
        sections.append(f"""
        <section class="mb-10" id="{tier}">
          <h2 class="text-xl font-bold mb-5 flex items-center gap-3" style="color: {TIER_COLOR[tier]}">
            {html.escape(TIER_LABEL[tier])} 
            <span class="bg-white/60 text-slate-700 text-sm px-3 py-0.5 rounded-full border border-slate-200 shadow-sm">{len(recs)}</span>
          </h2>
          {cards}
        </section>""")
    body = "\n".join(sections)
    
    # Modern Glassmorphism CSS + Tailwind via CDN
    doc = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>{html.escape(title)} — {html.escape(data['patient'])}</title>
    <script src="https://cdn.tailwindcss.com"></script>
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&display=swap" rel="stylesheet">
    <style>
        body {{
            font-family: 'Inter', sans-serif;
            background: radial-gradient(circle at 15% 50%, #f1f5f9, #e2e8f0);
            background-attachment: fixed;
            color: #1e293b;
        }}
        /* Pattern overlay for texture */
        body::before {{
            content: '';
            position: fixed;
            top: 0; left: 0; right: 0; bottom: 0;
            background-image: url("data:image/svg+xml,%3Csvg width='60' height='60' viewBox='0 0 60 60' xmlns='http://www.w3.org/2000/svg'%3E%3Cg fill='none' fill-rule='evenodd'%3E%3Cg fill='%2394a3b8' fill-opacity='0.05'%3E%3Cpath d='M36 34v-4h-2v4h-4v2h4v4h2v-4h4v-2h-4zm0-30V0h-2v4h-4v2h4v4h2V6h4V4h-4zM6 34v-4H4v4H0v2h4v4h2v-4h4v-2H6zM6 4V0H4v4H0v2h4v4h2V6h4V4H6z'/%3E%3C/g%3E%3C/g%3E%3C/svg%3E");
            z-index: -1;
        }}
        .metrics-grid label {{
            display: block;
            font-size: 10px;
            text-transform: uppercase;
            letter-spacing: 0.04em;
            color: #64748b;
            font-weight: 600;
            margin-bottom: 2px;
        }}
        .metrics-grid span {{
            font-size: 14px;
            color: #334155;
            font-weight: 400;
        }}
        @media print {{ 
            .noprint {{ display: none; }} 
            .variant-card {{ page-break-inside: avoid; box-shadow: none !important; border: 1px solid #ccc !important; }}
            body {{ background: white !important; }}
            body::before {{ display: none; }}
        }}
    </style>
</head>
<body class="antialiased min-h-screen">
    
    <header class="relative overflow-hidden bg-slate-900 text-white px-8 py-10 shadow-xl">
        <!-- Decorative gradients -->
        <div class="absolute top-0 right-0 -mr-20 -mt-20 w-96 h-96 rounded-full bg-indigo-500 blur-3xl opacity-20 pointer-events-none"></div>
        <div class="absolute bottom-0 left-0 -ml-20 -mb-20 w-80 h-80 rounded-full bg-blue-500 blur-3xl opacity-20 pointer-events-none"></div>
        
        <div class="relative z-10 max-w-7xl mx-auto">
            <h1 class="text-3xl font-bold tracking-tight mb-2">{html.escape(title)}</h1>
            <div class="text-slate-300 font-medium tracking-wide text-sm flex items-center flex-wrap gap-x-6 gap-y-2">
                <span>Patient: <b class="text-white">{html.escape(data['patient'])}</b></span>
                <span class="text-slate-600">|</span>
                <span>Domain: <span class="text-indigo-300">{html.escape(str(data['domain']))}</span></span>
                <span class="text-slate-600">|</span>
                <span>Ontology Panel: <b class="text-white">{data['panel_gene_count']}</b> genes</span>
            </div>
        </div>
    </header>

    <div class="max-w-7xl mx-auto px-6 py-8">
        
        <!-- Summary Stats -->
        <div class="grid grid-cols-2 md:grid-cols-5 gap-4 mb-8">
            <div class="bg-white/60 backdrop-blur-lg border border-white/40 shadow-sm rounded-2xl p-5 flex flex-col justify-center">
                <span class="text-3xl font-bold text-slate-800">{data['actionable_count']}</span>
                <span class="text-[11px] uppercase tracking-wider font-semibold text-slate-500 mt-1">Actionable Variants</span>
            </div>
            <div class="bg-white/60 backdrop-blur-lg border border-white/40 shadow-sm rounded-2xl p-5 flex flex-col justify-center">
                <span class="text-3xl font-bold" style="color: {TIER_COLOR['Tier1']}">{tc['Tier1']}</span>
                <span class="text-[11px] uppercase tracking-wider font-semibold text-slate-500 mt-1">Tier 1</span>
            </div>
            <div class="bg-white/60 backdrop-blur-lg border border-white/40 shadow-sm rounded-2xl p-5 flex flex-col justify-center">
                <span class="text-3xl font-bold" style="color: {TIER_COLOR['Tier2']}">{tc['Tier2']}</span>
                <span class="text-[11px] uppercase tracking-wider font-semibold text-slate-500 mt-1">Tier 2</span>
            </div>
            <div class="bg-white/60 backdrop-blur-lg border border-white/40 shadow-sm rounded-2xl p-5 flex flex-col justify-center">
                <span class="text-3xl font-bold" style="color: {TIER_COLOR['Tier3']}">{tc['Tier3']}</span>
                <span class="text-[11px] uppercase tracking-wider font-semibold text-slate-500 mt-1">Tier 3</span>
            </div>
            <div class="bg-white/60 backdrop-blur-lg border border-white/40 shadow-sm rounded-2xl p-5 flex flex-col justify-center">
                <span class="text-3xl font-bold text-slate-700">{data['scanned_panel_variants']}</span>
                <span class="text-[11px] uppercase tracking-wider font-semibold text-slate-500 mt-1">Panel Variants Scanned</span>
            </div>
        </div>

        <div class="flex items-center gap-4 mb-8 noprint">
            <div class="relative flex-1 max-w-md">
                <svg class="absolute left-3 top-1/2 -translate-y-1/2 w-4 h-4 text-slate-400" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M21 21l-6-6m2-5a7 7 0 11-14 0 7 7 0 0114 0z"></path></svg>
                <input id="flt" type="text" placeholder="Filter by gene or reason code..."
                    class="w-full pl-10 pr-4 py-2.5 bg-white/70 border border-slate-200 rounded-xl shadow-sm focus:outline-none focus:ring-2 focus:ring-indigo-500 focus:border-transparent transition-all"
                    oninput="filterCards(this.value)">
            </div>
            <button onclick="window.print()" class="px-5 py-2.5 bg-white/80 border border-slate-200 text-slate-700 rounded-xl shadow-sm hover:bg-white hover:shadow transition-all font-medium text-sm flex items-center gap-2">
                <svg class="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M17 17h2a2 2 0 002-2v-4a2 2 0 00-2-2H5a2 2 0 00-2 2v4a2 2 0 002 2h2m2 4h6a2 2 0 002-2v-4a2 2 0 00-2-2H9a2 2 0 00-2 2v4a2 2 0 002 2zm8-12V5a2 2 0 00-2-2H9a2 2 0 00-2 2v4h10z"></path></svg>
                Print / Export PDF
            </button>
        </div>

        {body}
        
    </div>

    <footer class="bg-slate-50 border-t border-slate-200 py-10 px-8 text-center text-slate-500 text-xs leading-relaxed mt-10">
        <div class="max-w-4xl mx-auto">
            Generated by <code class="bg-slate-100 text-slate-600 px-1 py-0.5 rounded">ontology_report</code> (domain: {html.escape(str(data['domain']))}).
            Gene selection is derived from OpenCRAVAT <code class="bg-slate-100 px-1 py-0.5 rounded">hpo</code> and <code class="bg-slate-100 px-1 py-0.5 rounded">go</code> annotators.
            Clinical evidence from <code class="bg-slate-100 px-1 py-0.5 rounded">clinvar</code>, <code class="bg-slate-100 px-1 py-0.5 rounded">clingen</code>, <code class="bg-slate-100 px-1 py-0.5 rounded">omim</code>, <code class="bg-slate-100 px-1 py-0.5 rounded">pharmgkb</code>.
            Deleteriousness from REVEL, AlphaMissense, BayesDel, MetaRNN, ESM1b, VARITY, SpliceAI, CADD plus configured domain-specific predictors.
            <br><span class="font-semibold text-slate-600 mt-2 block">Research/screening use — not a substitute for clinical diagnostic interpretation.</span>
        </div>
    </footer>

    <script>
    function filterCards(q){{
        q = q.trim().toLowerCase();
        document.querySelectorAll('.variant-card').forEach(function(c){{
            var hay = (c.dataset.gene + ' ' + c.dataset.reasons).toLowerCase();
            c.style.display = (!q || hay.indexOf(q) >= 0) ? '' : 'none';
        }});
    }}
    </script>
</body>
</html>"""
    with open(path, "w") as f:
        f.write(doc)


def main():
    ap = argparse.ArgumentParser(description="Render actionable report")
    ap.add_argument("--in-json", required=True)
    ap.add_argument("--out-html", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-text", required=True)
    args = ap.parse_args()
    data = json.load(open(args.in_json))
    write_html(data, args.out_html)
    write_tsv(data["records"], args.out_tsv)
    write_text(data, args.out_text)
    print(f"[render] HTML  -> {args.out_html}")
    print(f"[render] TSV   -> {args.out_tsv}")
    print(f"[render] text  -> {args.out_text}")


if __name__ == "__main__":
    main()
