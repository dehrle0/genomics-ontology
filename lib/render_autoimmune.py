#!/usr/bin/env python3
"""
render_autoimmune.py
A domain-specific renderer for the autoimmunity report. It reuses the generic
renderer's helpers (formatting, TSV, text) but adds two things the autoimmune
story needs and the generic tier view does not:

  1. A TRAIT-BURDEN VISUALIZATION — an inline SVG chart (no external libraries,
     works offline) summarising, across all kept variants, how many catalogued
     GWAS risk-allele associations point at each autoimmune trait, coloured by
     the strongest (smallest) p-value seen for that trait.

  2. LIVE STUDY-EVIDENCE CARDS — per variant, the current GWAS Catalog
     associations pulled by lib/enrich_report.py: trait, p-value, odds
     ratio / beta, risk allele, and a PubMed link to the underlying study.

Gene cards also surface the NCBI Gene description (from the same enrichment).
Everything degrades gracefully when enrichment was skipped or offline: the
report still renders, just without the live layers.
"""
import argparse
import html
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import render_report as rr  # noqa: E402  (shared helpers, TSV/text writers)

TIER_LABEL = {
    "Tier1": "Tier 1 — Monogenic / high-impact",
    "Tier2": "Tier 2 — Supported risk / VUS",
    "Tier3": "Tier 3 — Catalogued risk allele / monitor",
}
TIER_COLOR = rr.TIER_COLOR


# --------------------------------------------------------------------------- #
# Aggregation for the visualization
# --------------------------------------------------------------------------- #
def collect_traits(records):
    """Aggregate GWAS associations across all records into per-trait rollups."""
    traits = {}
    for r in records:
        se = r.get("study_evidence") or {}
        snp = se.get("snp") or {}
        gene = r.get("hugo")
        for study in snp.get("top", []) or []:
            pv = _pv(study.get("pvalue"))
            for t in study.get("traits", []) or []:
                d = traits.setdefault(t, {"count": 0, "min_p": None,
                                          "genes": set(), "rsids": set()})
                d["count"] += 1
                if gene:
                    d["genes"].add(gene)
                if snp.get("rsid"):
                    d["rsids"].add(snp["rsid"])
                if pv is not None and (d["min_p"] is None or pv < d["min_p"]):
                    d["min_p"] = pv
    out = []
    for name, d in traits.items():
        out.append({
            "trait": name, "count": d["count"], "min_p": d["min_p"],
            "genes": sorted(d["genes"]), "rsids": sorted(d["rsids"]),
        })
    out.sort(key=lambda x: (-x["count"], x["min_p"] if x["min_p"] is not None else 1.0))
    return out


def _pv(x):
    try:
        f = float(x)
        return f if f > 0 else None
    except (TypeError, ValueError):
        return None


def _p_color(min_p):
    """Colour a bar by association strength: deeper red = stronger (smaller p)."""
    if not min_p or min_p <= 0:
        return "#95a5a6"
    nlp = -math.log10(min_p)            # e.g. p=5e-8 -> 7.3
    nlp = max(0.0, min(nlp, 50.0))
    # interpolate light-blue (weak) -> crimson (strong) over 0..30
    t = min(nlp / 30.0, 1.0)
    r = int(41 + t * (192 - 41))
    g = int(128 + t * (57 - 128))
    b = int(185 + t * (43 - 185))
    return f"rgb({r},{g},{b})"


def _fmt_p(x):
    p = _pv(x)
    if p is None:
        return "n/a"
    if p < 1e-3:
        return f"{p:.1e}"
    return f"{p:.3g}"


def trait_chart_svg(trait_rows, top_n=14):
    """Inline SVG horizontal bar chart of the top autoimmune traits by number of
    catalogued risk-allele associations across this genome."""
    rows = trait_rows[:top_n]
    if not rows:
        return ('<p class="empty">No catalogued GWAS trait associations were '
                'available (run with study enrichment online to populate this '
                'chart).</p>')
    max_count = max(r["count"] for r in rows)
    row_h, gap, left, right, top = 26, 8, 260, 90, 14
    width = 900
    bar_area = width - left - right
    height = top + len(rows) * (row_h + gap) + 10
    parts = [f'<svg viewBox="0 0 {width} {height}" width="100%" '
             f'preserveAspectRatio="xMinYMin meet" role="img" '
             f'aria-label="Autoimmune trait burden chart" class="trait-chart">']
    y = top
    for r in rows:
        w = max(2, int(bar_area * (r["count"] / max_count)))
        color = _p_color(r["min_p"])
        label = html.escape(_truncate(r["trait"], 38))
        genes = html.escape(", ".join(r["genes"][:4]) + ("…" if len(r["genes"]) > 4 else ""))
        title = html.escape(f"{r['trait']} — {r['count']} associations, "
                            f"best p={_fmt_p(r['min_p'])}, genes: {', '.join(r['genes'])}")
        parts.append(
            f'<g><title>{title}</title>'
            f'<text x="{left - 8}" y="{y + row_h * 0.68}" text-anchor="end" '
            f'class="t-lbl">{label}</text>'
            f'<rect x="{left}" y="{y}" width="{w}" height="{row_h}" rx="4" '
            f'fill="{color}"><title>{title}</title></rect>'
            f'<text x="{left + w + 6}" y="{y + row_h * 0.68}" class="t-cnt">'
            f'{r["count"]}</text>'
            f'<text x="{left + 4}" y="{y + row_h * 0.68}" class="t-gene">{genes}</text>'
            f'</g>')
        y += row_h + gap
    parts.append('</svg>')
    return "\n".join(parts)


def _truncate(s, n):
    s = str(s)
    return s if len(s) <= n else s[: n - 1] + "…"


# --------------------------------------------------------------------------- #
# Per-variant study evidence + gene description
# --------------------------------------------------------------------------- #
def _pubmed_link(pub):
    if not pub:
        return ""
    pid = pub.get("pubmed_id")
    cite = " ".join(x for x in [pub.get("author"), pub.get("journal"),
                                (pub.get("date") or "")[:4]] if x)
    label = html.escape(cite or (f"PMID {pid}" if pid else "study"))
    if pid:
        return (f'<a href="https://pubmed.ncbi.nlm.nih.gov/{html.escape(str(pid))}/" '
                f'target="_blank" title="{html.escape(pub.get("title") or "")}">'
                f'{label} &#8599;</a>')
    return label


def _study_rows(rec):
    se = rec.get("study_evidence") or {}
    snp = se.get("snp") or {}
    top = snp.get("top", []) or []
    if not top:
        return ""
    rows = []
    for s in top:
        traits = html.escape(", ".join(s.get("traits", []) or []) or "-")
        pcell = _fmt_p(s.get("pvalue"))
        orv = s.get("or_beta")
        orcell = html.escape(f"{orv}") if orv not in (None, "") else "-"
        risk = html.escape(str(s.get("risk_allele") or "-"))
        rows.append(
            f"<tr><td class='trait'>{traits}</td><td class='pv'>{pcell}</td>"
            f"<td>{orcell}</td><td class='ra'>{risk}</td>"
            f"<td class='pub'>{_pubmed_link(s.get('pubmed'))}</td></tr>")
    n = snp.get("n_associations")
    more = (f'<div class="study-more">Showing {len(top)} of {n} catalogued '
            f'associations for {html.escape(str(snp.get("rsid")))}.</div>'
            if n and n > len(top) else "")
    return (f'<div class="studies"><div class="studies-h">Current GWAS Catalog '
            f'evidence</div><table class="study-tbl"><thead><tr>'
            f'<th>Trait</th><th>p-value</th><th>OR/&beta;</th><th>Risk allele</th>'
            f'<th>Study</th></tr></thead><tbody>'
            + "".join(rows) + f'</tbody></table>{more}</div>')


def _gene_summary_block(rec):
    gi = rec.get("gene_info") or {}
    desc = gi.get("description")
    summ = gi.get("summary")
    gid = gi.get("ncbi_gene_id")
    loc = gi.get("map_location")
    if not (desc or summ or gid):
        return ""
    link = (f' <a href="https://www.ncbi.nlm.nih.gov/gene/{html.escape(str(gid))}" '
            f'target="_blank">NCBI&#8599;</a>') if gid else ""
    loc_txt = f' <span class="cyto">{html.escape(loc)}</span>' if loc else ""
    head = (f'<div class="gene-desc"><span class="gd-label">NCBI Gene</span>'
            f'{html.escape(desc or "")}{loc_txt}{link}</div>') if (desc or gid) else ""
    body = (f'<div class="gene-summary">{html.escape(_truncate(summ, 420))}</div>'
            if summ else "")
    return head + body


def _card(r):
    ev = r.get("evidence", {})
    so = rr.SO_NAME.get(r.get("so"), r.get("so") or "?")
    gene = html.escape(r.get("hugo") or "?")
    zyg = ev.get("zygosity")
    rsid = r.get("rsid")
    rsid_html = (f'<a href="https://www.ncbi.nlm.nih.gov/snp/{html.escape(str(rsid))}" '
                 f'target="_blank">{html.escape(str(rsid))}</a>'
                 ) if rsid and str(rsid).startswith("rs") else "-"
    hpo_ctx = ", ".join(ev.get("hpo_context", []) or []) or "-"
    go_ctx = ", ".join(ev.get("go_context", []) or []) or "-"
    return f"""
    <div class="card" data-gene="{gene}" data-reasons="{html.escape(' '.join(r.get('reason_codes', [])))}">
      <div class="card-head">
        <span class="gene">{gene}</span>
        {rr._zyg_badge(zyg)}
        <span class="loc">{html.escape(str(r.get('chrom','')))}:{html.escape(str(r.get('pos','')))}
          {html.escape(rr._fmt_allele(r.get('ref')))}&gt;{html.escape(rr._fmt_allele(r.get('alt')))}</span>
        <span class="so">{html.escape(so)}</span>
        <span class="ach">{html.escape(rr._fmt_allele(r.get('achange') or r.get('cchange') or '', 28))}</span>
      </div>
      {_gene_summary_block(r)}
      <div class="grid">
        <div><label>Zygosity</label>{html.escape(zyg or '-')}</div>
        <div><label>dbSNP</label>{rsid_html}</div>
        <div><label>gnomAD4 AF</label>{rr._fmt_af(r.get('gnomad4_af'))}</div>
        <div><label>All of Us AF</label>{rr._fmt_af(r.get('allofus_af'))}</div>
        <div><label>ClinVar</label>{rr._clinvar_link(r.get('clinvar_id'), r.get('clinvar_sig'))}</div>
        <div><label>REVEL</label>{html.escape(str(r.get('revel') or '-'))}</div>
        <div><label>AlphaMissense</label>{html.escape(str(r.get('am_path') or '-'))}</div>
        <div><label>Panel support</label>{html.escape(str(ev.get('panel_support') or '-'))}/2</div>
      </div>
      <div class="onto">
        <div><label>HPO phenotype context</label>{html.escape(hpo_ctx)}</div>
        <div><label>GO function context</label>{html.escape(go_ctx)}</div>
      </div>
      {_study_rows(r)}
      <div class="reasons">{rr._reason_badges(r.get('reason_codes', []))}</div>
    </div>"""


# --------------------------------------------------------------------------- #
# Document
# --------------------------------------------------------------------------- #
def write_html(data, path):
    tc = data["tier_counts"]
    title = data.get("report_title", "Autoimmune Risk & Evidence Report")
    records = data["records"]
    trait_rows = collect_traits(records)
    chart = trait_chart_svg(trait_rows)
    enr = data.get("enrichment", {})

    n_traits = len(trait_rows)
    n_studies = sum(len((r.get("study_evidence", {}).get("snp") or {}).get("top", []) or [])
                    for r in records)

    sections = []
    for tier in ("Tier1", "Tier2", "Tier3"):
        recs = [r for r in records if r["tier"] == tier]
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

    enr_note = ""
    if enr:
        if enr.get("offline"):
            enr_note = ("Live enrichment ran in <b>offline</b> mode — gene "
                        "descriptions and study evidence shown are from cache only.")
        else:
            enr_note = (f"Live enrichment: {enr.get('records_with_gene_info', 0)} gene "
                        f"descriptions and {enr.get('records_with_study_evidence', 0)} "
                        f"variants with study evidence "
                        f"({enr.get('remote_calls', 0)} API calls, "
                        f"{enr.get('remote_errors', 0)} errors).")

    doc = f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html.escape(title)} — {html.escape(data['patient'])}</title>
<style>
:root {{ --bg:#f4f6f9; --card:#fff; --ink:#22303f; --muted:#6b7c8f; }}
* {{ box-sizing:border-box; }}
body {{ font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;
  margin:0; background:var(--bg); color:var(--ink); }}
header {{ background:linear-gradient(135deg,#4a1c40,#7b2d5e,#b0355f); color:#fff; padding:28px 32px; }}
header h1 {{ margin:0 0 6px; font-size:22px; }}
header .sub {{ opacity:.9; font-size:14px; }}
.summary {{ display:flex; gap:16px; flex-wrap:wrap; padding:18px 32px 6px; }}
.stat {{ background:var(--card); border-radius:10px; padding:14px 18px; min-width:120px;
  box-shadow:0 1px 3px rgba(0,0,0,.08); }}
.stat b {{ display:block; font-size:26px; }}
.stat span {{ color:var(--muted); font-size:12px; text-transform:uppercase; letter-spacing:.04em; }}
.viz {{ margin:8px 32px 4px; background:var(--card); border-radius:12px; padding:18px 22px;
  box-shadow:0 1px 3px rgba(0,0,0,.08); }}
.viz h2 {{ margin:0 0 4px; font-size:16px; }}
.viz .cap {{ color:var(--muted); font-size:12.5px; margin:0 0 12px; }}
.trait-chart .t-lbl {{ font-size:12.5px; fill:var(--ink); }}
.trait-chart .t-cnt {{ font-size:12px; fill:var(--muted); }}
.trait-chart .t-gene {{ font-size:10.5px; fill:#fff; opacity:.85; }}
.legend {{ font-size:11.5px; color:var(--muted); margin-top:8px; display:flex; gap:8px; align-items:center; }}
.legend .bar {{ height:10px; width:160px; border-radius:5px;
  background:linear-gradient(90deg, rgb(41,128,185), rgb(192,57,43)); display:inline-block; }}
.enr-note {{ margin:0 32px 6px; font-size:12.5px; color:var(--muted); }}
.controls {{ padding:6px 32px 8px; }}
.controls input {{ padding:8px 12px; border:1px solid #cdd6e0; border-radius:8px; width:280px; }}
section.tier {{ padding:8px 32px 24px; }}
section.tier h2 {{ font-size:17px; padding-left:12px; }}
section.tier h2 .count {{ background:#e8edf3; border-radius:12px; padding:1px 10px; font-size:13px; margin-left:8px; }}
.card {{ background:var(--card); border-radius:10px; padding:14px 16px; margin:10px 0;
  box-shadow:0 1px 3px rgba(0,0,0,.07); }}
.card-head {{ display:flex; gap:12px; align-items:baseline; flex-wrap:wrap; border-bottom:1px solid #eef2f6; padding-bottom:8px; }}
.card-head .gene {{ font-weight:700; font-size:16px; }}
.card-head .loc {{ font-family:ui-monospace,Menlo,monospace; color:var(--muted); font-size:13px; }}
.card-head .so {{ background:#eef2f6; border-radius:6px; padding:1px 8px; font-size:12px; }}
.card-head .ach {{ color:#7b2d5e; font-size:13px; font-family:ui-monospace,monospace; }}
.zyg {{ font-size:11px; font-weight:600; border-radius:6px; padding:1px 8px; }}
.z-het {{ background:#fef5e7; color:#b9770e; }}
.z-hom {{ background:#fdecea; color:#c0392b; }}
.z-hemi {{ background:#f4ecf7; color:#7d3c98; }}
.z-other {{ background:#eef2f6; color:#566573; }}
.gene-desc {{ font-size:13px; color:#34495e; margin:8px 0 2px; line-height:1.4; }}
.gene-desc .gd-label {{ font-size:10px; text-transform:uppercase; letter-spacing:.04em;
  color:#fff; background:#7b2d5e; border-radius:4px; padding:1px 6px; margin-right:6px; }}
.gene-desc .cyto {{ color:var(--muted); font-family:ui-monospace,monospace; }}
.gene-summary {{ font-size:12.5px; color:var(--muted); margin:2px 0 6px; line-height:1.45; }}
.grid {{ display:grid; grid-template-columns:repeat(4,1fr); gap:6px 18px; margin:10px 0; }}
.grid label, .onto label {{ display:block; font-size:11px; color:var(--muted); text-transform:uppercase; letter-spacing:.03em; }}
.grid div, .onto div {{ font-size:14px; }}
.onto {{ display:grid; grid-template-columns:1fr 1fr; gap:6px 18px; margin:6px 0 10px; }}
.studies {{ margin:8px 0 6px; }}
.studies-h {{ font-size:11px; text-transform:uppercase; letter-spacing:.04em; color:#7b2d5e; font-weight:600; margin-bottom:4px; }}
.study-tbl {{ width:100%; border-collapse:collapse; font-size:12.5px; }}
.study-tbl th {{ text-align:left; color:var(--muted); font-weight:600; border-bottom:1px solid #eef2f6; padding:3px 6px; font-size:11px; }}
.study-tbl td {{ padding:3px 6px; border-bottom:1px solid #f4f6f9; vertical-align:top; }}
.study-tbl td.trait {{ max-width:320px; }}
.study-tbl td.pv {{ font-family:ui-monospace,monospace; color:#c0392b; }}
.study-tbl td.ra {{ font-family:ui-monospace,monospace; }}
.study-more {{ font-size:11px; color:var(--muted); margin-top:4px; }}
.reasons {{ margin-top:6px; }}
.badge {{ display:inline-block; font-size:11px; padding:2px 8px; border-radius:10px; margin:2px 3px 0 0; }}
.b-pheno {{ background:#eafaf1; color:#1e8449; }}
.b-geno {{ background:#eaf2fb; color:#2471a3; }}
.b-strong {{ background:#fdecea; color:#c0392b; font-weight:600; }}
.b-warn {{ background:#fef5e7; color:#b9770e; }}
.empty {{ color:var(--muted); font-style:italic; }}
footer {{ padding:18px 32px 40px; color:var(--muted); font-size:12px; }}
.noprint {{}} @media print {{ .noprint {{ display:none; }} .card {{ page-break-inside:avoid; }} }}
</style></head>
<body>
<header>
  <h1>{html.escape(title)}</h1>
  <div class="sub">Patient: <b>{html.escape(data['patient'])}</b> &nbsp;|&nbsp; Domain: {html.escape(str(data['domain']))}
   &nbsp;|&nbsp; Panel derived from HPO autoimmune phenotypes + GO immune functions ({data['panel_gene_count']} genes)</div>
</header>
<div class="summary">
  <div class="stat"><b>{data['actionable_count']}</b><span>Risk variants</span></div>
  <div class="stat"><b style="color:{TIER_COLOR['Tier1']}">{tc['Tier1']}</b><span>Tier 1</span></div>
  <div class="stat"><b style="color:{TIER_COLOR['Tier2']}">{tc['Tier2']}</b><span>Tier 2</span></div>
  <div class="stat"><b style="color:{TIER_COLOR['Tier3']}">{tc['Tier3']}</b><span>Tier 3</span></div>
  <div class="stat"><b>{n_traits}</b><span>Traits implicated</span></div>
  <div class="stat"><b>{n_studies}</b><span>GWAS associations</span></div>
</div>
{f'<div class="enr-note noprint">{enr_note}</div>' if enr_note else ''}
<div class="viz">
  <h2>Autoimmune trait burden across your risk variants</h2>
  <p class="cap">Bars count catalogued GWAS risk-allele associations pointing at each
  trait; colour encodes the strongest reported association (deeper red = smaller p-value).
  Hover a bar for the contributing genes and best p-value.</p>
  {chart}
  <div class="legend">weak <span class="bar"></span> strong association &nbsp;·&nbsp;
    p-value scale (&minus;log<sub>10</sub> p, 0&ndash;30)</div>
</div>
<div class="controls noprint">
  <input id="flt" type="text" placeholder="Filter by gene or reason code…"
    oninput="filterCards(this.value)">
  <button onclick="window.print()">Print / Export PDF</button>
</div>
{body}
<footer>
  Generated by <code>ontology_report</code> (domain: {html.escape(str(data['domain']))},
  renderer: autoimmune). Gene panel derived from the OpenCRAVAT <code>hpo</code> and
  <code>go</code> annotators; catalogued risk alleles from <code>gwas_catalog</code>;
  live study evidence from the EBI GWAS Catalog REST API; gene descriptions from
  NCBI Gene (E-utilities). Autoimmune risk is polygenic and context-dependent —
  this report summarises evidence, it is not a diagnosis.
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
    ap = argparse.ArgumentParser(description="Render autoimmune report (viz + live studies)")
    ap.add_argument("--in-json", required=True)
    ap.add_argument("--out-html", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-text", required=True)
    args = ap.parse_args()
    data = json.load(open(args.in_json))
    write_html(data, args.out_html)
    rr.write_tsv(data["records"], args.out_tsv)
    rr.write_text(data, args.out_text)
    print(f"[render-autoimmune] HTML  -> {args.out_html}")
    print(f"[render-autoimmune] TSV   -> {args.out_tsv}")
    print(f"[render-autoimmune] text  -> {args.out_text}")


if __name__ == "__main__":
    main()
