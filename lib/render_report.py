#!/usr/bin/env python3
"""
render_report.py
Turn the actionable-variant JSON (from ontology_filter.py) into human
deliverables: a styled interactive HTML report, a flat TSV, and a text summary.
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
TIER_COLOR = {"Tier1": "#c0392b", "Tier2": "#d68910", "Tier3": "#2471a3"}

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
        cls = "b-geno"
        if rc.startswith("HPO_") or rc.startswith("GO_") or rc.startswith("CLIN") or rc in ("OMIM_DISEASE", "ARRVARS_KNOWN"):
            cls = "b-pheno"
        if rc in ("PVS1_HAPLOINSUFFICIENT", "PP3_CONSENSUS", "SPLICEAI_HIGH", "PM2_RARE"):
            cls = "b-strong"
        if rc == "COMMON_AF_FLAG":
            cls = "b-warn"
        out.append(f'<span class="badge {cls}">{html.escape(rc)}</span>')
    return " ".join(out)


def _zyg_badge(zyg):
    if not zyg:
        return ""
    cls = {"Heterozygous": "z-het", "Homozygous": "z-hom",
           "Hemizygous": "z-hemi"}.get(zyg, "z-other")
    return f'<span class="zyg {cls}">{html.escape(zyg)}</span>'


def _gene_desc_block(r):
    """NCBI Gene description subtitle (present only after enrichment)."""
    gi = r.get("gene_info")
    if not gi:
        return ""
    gid = gi.get("ncbi_gene_id")
    desc = gi.get("description") or ""
    loc = gi.get("map_location")
    link = (f'<a href="https://www.ncbi.nlm.nih.gov/gene/{html.escape(str(gid))}" '
            f'target="_blank">NCBI&#8599;</a>') if gid else ""
    loc_txt = f' <span class="cyto">{html.escape(loc)}</span>' if loc else ""
    if not (desc or gid):
        return ""
    return (f'<div class="gene-desc"><span class="gd-label">NCBI Gene</span> '
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
                 f'target="_blank">{html.escape(str(rsid))}</a>') if rsid and str(rsid).startswith("rs") else (html.escape(str(rsid)) if rsid else "-")
    return f"""
    <div class="card" data-gene="{gene}" data-reasons="{html.escape(' '.join(r.get('reason_codes', [])))}">
      <div class="card-head">
        <span class="gene">{gene}</span>
        {_zyg_badge(zyg)}
        <span class="loc">{html.escape(str(r.get('chrom','')))}:{html.escape(str(r.get('pos','')))}
          {html.escape(_fmt_allele(r.get('ref')))}&gt;{html.escape(_fmt_allele(r.get('alt')))}</span>
        <span class="so">{html.escape(so)}</span>
        <span class="ach">{html.escape(_fmt_allele(r.get('achange') or r.get('cchange') or '', 28))}</span>
      </div>
      {_gene_desc_block(r)}
      <div class="grid">
        <div><label>Zygosity</label>{html.escape(zyg or '-')}</div>
        <div><label>Variant allele frac</label>{_fmt_af(vaf) if vaf is not None else '-'}</div>
        <div><label>dbSNP</label>{rsid_html}</div>
        <div><label>gnomAD4 AF</label>{_fmt_af(r.get('gnomad4_af'))}</div>
        <div><label>All of Us AF</label>{_fmt_af(r.get('allofus_af'))}</div>
        <div><label>ClinVar</label>{_clinvar_link(r.get('clinvar_id'), r.get('clinvar_sig'))}</div>
        <div><label>REVEL</label>{html.escape(str(r.get('revel') or '-'))}</div>
        <div><label>AlphaMissense</label>{html.escape(str(r.get('am_path') or '-'))}</div>
        <div><label>SpliceAI max</label>{html.escape(str(ev.get('spliceai_max') if ev.get('spliceai_max') is not None else '-'))}</div>
        <div><label>CADD phred</label>{html.escape(str(r.get('cadd_phred') or '-'))}</div>
        <div><label>Panel support</label>{html.escape(str(ev.get('panel_support') or '-'))}/2</div>
      </div>
      <div class="onto">
        <div><label>HPO phenotype context</label>{html.escape(hpo_ctx)}
          &nbsp;<a href="{hpo_gene_link}" target="_blank">HPO&#8599;</a></div>
        <div><label>GO function context</label>{html.escape(go_ctx)}
          &nbsp;<a href="{gene_link}" target="_blank">GenCC&#8599;</a></div>
      </div>
      <div class="reasons">{_reason_badges(r.get('reason_codes', []))}</div>
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
        <section class="tier" id="{tier}">
          <h2 style="border-left:6px solid {TIER_COLOR[tier]}">
            {html.escape(TIER_LABEL[tier])} <span class="count">{len(recs)}</span></h2>
          {cards}
        </section>""")
    body = "\n".join(sections)
    doc = f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html.escape(title)} — {html.escape(data['patient'])}</title>
<style>
:root {{ --bg:#f4f6f9; --card:#fff; --ink:#22303f; --muted:#6b7c8f; }}
* {{ box-sizing:border-box; }}
body {{ font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;
  margin:0; background:var(--bg); color:var(--ink); }}
header {{ background:linear-gradient(135deg,#1b2a3a,#2c5364); color:#fff; padding:28px 32px; }}
header h1 {{ margin:0 0 6px; font-size:22px; }}
header .sub {{ opacity:.85; font-size:14px; }}
.summary {{ display:flex; gap:16px; flex-wrap:wrap; padding:18px 32px; }}
.stat {{ background:var(--card); border-radius:10px; padding:14px 18px; min-width:120px;
  box-shadow:0 1px 3px rgba(0,0,0,.08); }}
.stat b {{ display:block; font-size:26px; }}
.stat span {{ color:var(--muted); font-size:12px; text-transform:uppercase; letter-spacing:.04em; }}
.controls {{ padding:0 32px 8px; }}
.controls input {{ padding:8px 12px; border:1px solid #cdd6e0; border-radius:8px; width:280px; }}
section.tier {{ padding:8px 32px 24px; }}
section.tier h2 {{ font-size:17px; padding-left:12px; }}
section.tier h2 .count {{ background:#e8edf3; border-radius:12px; padding:1px 10px; font-size:13px; margin-left:8px; }}
.card {{ background:var(--card); border-radius:10px; padding:14px 16px; margin:10px 0;
  box-shadow:0 1px 3px rgba(0,0,0,.07); }}
.card-head {{ display:flex; gap:14px; align-items:baseline; flex-wrap:wrap; border-bottom:1px solid #eef2f6; padding-bottom:8px; }}
.card-head .gene {{ font-weight:700; font-size:16px; }}
.card-head .loc {{ font-family:ui-monospace,Menlo,monospace; color:var(--muted); font-size:13px; }}
.card-head .so {{ background:#eef2f6; border-radius:6px; padding:1px 8px; font-size:12px; }}
.card-head .ach {{ color:#2c5364; font-size:13px; font-family:ui-monospace,monospace; }}
.zyg {{ font-size:11px; font-weight:600; border-radius:6px; padding:1px 8px; }}
.z-het {{ background:#fef5e7; color:#b9770e; }}
.z-hom {{ background:#fdecea; color:#c0392b; }}
.z-hemi {{ background:#f4ecf7; color:#7d3c98; }}
.z-other {{ background:#eef2f6; color:#566573; }}
.gene-desc {{ font-size:13px; color:#34495e; margin:6px 0 2px; line-height:1.4; }}
.gene-desc .gd-label {{ font-size:10px; text-transform:uppercase; letter-spacing:.04em;
  color:#fff; background:#2c5364; border-radius:4px; padding:1px 6px; margin-right:6px; }}
.gene-desc .cyto {{ color:var(--muted); font-family:ui-monospace,monospace; }}
.grid {{ display:grid; grid-template-columns:repeat(4,1fr); gap:6px 18px; margin:10px 0; }}
.grid label, .onto label {{ display:block; font-size:11px; color:var(--muted); text-transform:uppercase; letter-spacing:.03em; }}
.grid div, .onto div {{ font-size:14px; }}
.onto {{ display:grid; grid-template-columns:1fr 1fr; gap:6px 18px; margin:6px 0 10px; }}
.reasons {{ margin-top:6px; }}
.badge {{ display:inline-block; font-size:11px; padding:2px 8px; border-radius:10px; margin:2px 3px 0 0; }}
.b-pheno {{ background:#eafaf1; color:#1e8449; }}
.b-geno {{ background:#eaf2fb; color:#2471a3; }}
.b-strong {{ background:#fdecea; color:#c0392b; font-weight:600; }}
.b-warn {{ background:#fef5e7; color:#b9770e; }}
footer {{ padding:18px 32px 40px; color:var(--muted); font-size:12px; }}
.noprint {{}} @media print {{ .noprint {{ display:none; }} .card {{ page-break-inside:avoid; }} }}
</style></head>
<body>
<header>
  <h1>{html.escape(title)}</h1>
  <div class="sub">Patient: <b>{html.escape(data['patient'])}</b> &nbsp;|&nbsp; Domain: {html.escape(str(data['domain']))}
   &nbsp;|&nbsp; Gene panel derived from HPO + GO ontologies ({data['panel_gene_count']} genes)</div>
</header>
<div class="summary">
  <div class="stat"><b>{data['actionable_count']}</b><span>Actionable</span></div>
  <div class="stat"><b style="color:{TIER_COLOR['Tier1']}">{tc['Tier1']}</b><span>Tier 1</span></div>
  <div class="stat"><b style="color:{TIER_COLOR['Tier2']}">{tc['Tier2']}</b><span>Tier 2</span></div>
  <div class="stat"><b style="color:{TIER_COLOR['Tier3']}">{tc['Tier3']}</b><span>Tier 3</span></div>
  <div class="stat"><b>{data['scanned_panel_variants']}</b><span>Panel variants scanned</span></div>
</div>
<div class="controls noprint">
  <input id="flt" type="text" placeholder="Filter by gene or reason code…"
    oninput="filterCards(this.value)">
  <button onclick="window.print()">Print / Export PDF</button>
</div>
{body}
<footer>
  Generated by <code>ontology_report</code> (domain: {html.escape(str(data['domain']))}).
  Gene selection is derived from the OpenCRAVAT <code>hpo</code> and <code>go</code>
  annotators; clinical evidence from <code>clinvar</code>, <code>clingen</code>,
  <code>omim</code>; deleteriousness from REVEL, AlphaMissense, BayesDel, MetaRNN,
  ESM1b, VARITY, SpliceAI, CADD plus any configured domain-specific predictors.
  <br>Research/screening use — not a substitute for clinical diagnostic interpretation.
</footer>
<script>
function filterCards(q){{
  q=q.trim().toLowerCase();
  document.querySelectorAll('.card').forEach(function(c){{
    var hay=(c.dataset.gene+' '+c.dataset.reasons).toLowerCase();
    c.style.display = (!q || hay.indexOf(q)>=0) ? '' : 'none';
  }});
}}
</script>
</body></html>"""
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
