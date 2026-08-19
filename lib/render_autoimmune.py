#!/usr/bin/env python3
# render_autoimmune.py
import argparse
import html
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import render_report as rr

TIER_LABEL = {
    "Tier1": "Tier 1 — Monogenic / high-impact",
    "Tier2": "Tier 2 — Supported risk / VUS",
    "Tier3": "Tier 3 — Catalogued risk allele / monitor",
}
TIER_COLOR = rr.TIER_COLOR


def collect_traits(records):
    traits = {}
    for r in records:
        ev = r.get("evidence", {}) or {}
        gene = r.get("hugo")
        rsid = r.get("rsid")
        
        # 1. OpenCRAVAT gwas_catalog table in SQLite
        gwas_dis = ev.get("gwas_disease") or r.get("gwas_disease") or r.get("gwas_catalog__disease")
        gwas_pv = ev.get("gwas_pval") or r.get("gwas_pval") or r.get("gwas_catalog__pval")
        if gwas_dis and str(gwas_dis).strip() not in ("", "-", "None"):
            for t in str(gwas_dis).split(";"):
                t_clean = t.strip()
                if t_clean:
                    d = traits.setdefault(t_clean, {"count": 0, "min_p": None, "genes": set(), "rsids": set()})
                    d["count"] += 1
                    if gene:
                        d["genes"].add(gene)
                    if rsid:
                        d["rsids"].add(rsid)
                    pv = _pv(gwas_pv)
                    if pv is not None and (d["min_p"] is None or pv < d["min_p"]):
                        d["min_p"] = pv

        # 2. Live online GWAS enrichment fallback
        se = r.get("study_evidence") or {}
        snp = se.get("snp") or {}
        for study in snp.get("top", []) or []:
            pv = _pv(study.get("pvalue"))
            for t in study.get("traits", []) or []:
                d = traits.setdefault(t, {"count": 0, "min_p": None, "genes": set(), "rsids": set()})
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
    if not min_p or min_p <= 0:
        return "#95a5a6"
    nlp = -math.log10(min_p)
    nlp = max(0.0, min(nlp, 50.0))
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
    rows = trait_rows[:top_n]
    if not rows:
        return ('<p class="empty">No catalogued GWAS trait associations were '
                'available (run with study enrichment online to populate this chart).</p>')
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


def _pubmed_link(pub):
    if not pub:
        return ""
    pid = pub.get("pubmed_id")
    cite = " ".join(x for x in [pub.get("author"), pub.get("journal"), (pub.get("date") or "")[:4]] if x)
    label = html.escape(cite or (f"PMID {pid}" if pid else "study"))
    if pid:
        return (f'<a href="https://pubmed.ncbi.nlm.nih.gov/{html.escape(str(pid))}/" '
                f'target="_blank" title="{html.escape(pub.get("title") or "")}">'
                f'{label} &#8599;</a>')
    return label


def _study_rows(rec):
    ev = rec.get("evidence", {}) or {}
    gwas_dis = ev.get("gwas_disease") or rec.get("gwas_disease") or rec.get("gwas_catalog__disease")
    gwas_pv = ev.get("gwas_pval") or rec.get("gwas_pval") or rec.get("gwas_catalog__pval")
    gwas_or = ev.get("gwas_or_beta") or rec.get("gwas_or_beta") or rec.get("gwas_catalog__or_beta")
    gwas_pmid = ev.get("gwas_pmid") or rec.get("gwas_pmid") or rec.get("gwas_catalog__pmid")
    gwas_risk = ev.get("gwas_risk_allele") or rec.get("gwas_risk_allele") or rec.get("gwas_catalog__risk_allele")
    rsid = rec.get("rsid")

    rows = []
    if gwas_dis and str(gwas_dis).strip() not in ("", "-", "None"):
        traits = html.escape(str(gwas_dis).replace(";", ", "))
        pcell = _fmt_p(gwas_pv)
        orcell = html.escape(str(gwas_or)) if gwas_or not in (None, "") else "-"
        risk = html.escape(str(gwas_risk or "-"))
        rows.append(
            f"<tr><td class='trait'>{traits}</td><td class='pv'>{pcell}</td>"
            f"<td>{orcell}</td><td class='ra'>{risk}</td>"
            f"<td class='pub'>{rr._pmid_link(gwas_pmid)}</td></tr>")

    se = rec.get("study_evidence") or {}
    snp = se.get("snp") or {}
    top = snp.get("top", []) or []
    for s in top:
        traits = html.escape(", ".join(s.get("traits", []) or []) or "-")
        pcell = _fmt_p(s.get("pvalue"))
        orv = s.get("or_beta")
        orcell = html.escape(f"{orv}") if orv not in (None, "") else "-"
        risk = html.escape(str(s.get("risk_allele") or "-"))
        rows.append(
            f"<tr><td class='trait'>{traits}</td><td class='pv'>{pcell}</td>"
            f"<td>{orcell}</td><td class='ra'>{risk}</td>"
            f"<td class='pub'>{rr._pmid_link(s.get('pubmed'))}</td></tr>")

    if not rows:
        return ""
    return (f'<div class="studies" style="margin-top:12px; background:#f8fafc; padding:12px; border-radius:8px; border:1px solid #e2e8f0;"><div class="studies-h" style="font-weight:700; font-size:12px; color:#334155; margin-bottom:6px;">Current GWAS Catalog Evidence</div><table class="study-tbl" style="width:100%; font-size:12px; border-collapse:collapse;"><thead><tr style="text-align:left; color:#64748b; font-size:11px;">'
            f'<th>Trait</th><th>p-value</th><th>OR/&beta;</th><th>Risk allele</th>'
            f'<th>Study</th></tr></thead><tbody>'
            + "".join(rows) + f'</tbody></table></div>')


def _gene_summary_block(rec):
    gi = rec.get("gene_info") or {}
    desc = gi.get("description")
    summ = gi.get("summary")
    gid = gi.get("ncbi_gene_id")
    loc = gi.get("map_location")
    omim_id = rec.get("omim_id")
    clin_dis = rec.get("clinvar_disease")
    
    parts = []
    if desc or gid:
        link = (f' <a href="https://www.ncbi.nlm.nih.gov/gene/{html.escape(str(gid))}" '
                f'target="_blank" style="color:#2563eb; font-weight:600; text-decoration:none;">NCBI&#8599;</a>') if gid else ""
        loc_txt = f' <span style="font-family:monospace; color:#64748b; font-weight:normal;">[{html.escape(loc)}]</span>' if loc else ""
        parts.append(f'<div class="gene-desc-bold"><span class="gd-label">NCBI Gene</span>'
                     f'<strong>{html.escape(desc or "")}</strong>{loc_txt}{link}</div>')
    
    if omim_id or clin_dis:
        omim_parts = []
        if omim_id:
            omim_link = f'<a href="https://omim.org/entry/{html.escape(str(omim_id))}" target="_blank" style="color:#7c3aed; font-weight:700; text-decoration:none;">OMIM #{html.escape(str(omim_id))}&#8599;</a>'
            omim_parts.append(omim_link)
        if clin_dis:
            omim_parts.append(f'<span style="color:#334155; font-weight:500;">{html.escape(str(clin_dis))}</span>')
        parts.append(f'<div class="omim-block"><span class="omim-label">OMIM Clinical Synopsis</span>'
                     f'{" &nbsp;|&nbsp; ".join(omim_parts)}</div>')
                     
    if summ and summ != desc:
        parts.append(f'<div class="gene-summary">{html.escape(_truncate(summ, 450))}</div>')
        
    return "".join(parts)


def _phase_badge(ev):
    ph = ev.get("phasing") or ""
    origin = ev.get("phase_origin")
    if not ph or ph == "Unphased (Short-Read WGS)":
        return '<span class="px-2 py-0.5 rounded text-[11px] font-medium bg-slate-100 text-slate-600 border border-slate-200">Unphased (40x WGS)</span>'
    cls = "bg-purple-100 text-purple-800 border-purple-200"
    if origin == "Maternal":
        cls = "bg-pink-100 text-pink-800 border-pink-200 font-semibold"
    elif origin == "Paternal":
        cls = "bg-blue-100 text-blue-800 border-blue-200 font-semibold"
    label = html.escape(ph)
    return f'<span class="px-2 py-0.5 rounded text-[11px] border {cls}" title="Parental / chromosomal phase assignment">{label}</span>'


def _card(r):
    ev = r.get("evidence", {})
    so = rr.SO_NAME.get(r.get("so"), r.get("so") or "?")
    gene = html.escape(r.get("hugo") or "?")
    gene_link = f'https://search.thegencc.org/genes?q={gene}'
    hpo_gene_link = f'https://hpo.jax.org/app/browse/search?q={gene}&navFilter=all'
    zyg = ev.get("zygosity")
    vaf = ev.get("vaf")
    
    qual = ev.get("qual") or r.get("phred") or r.get("qual") or r.get("vcfinfo__phred")
    alt_reads = ev.get("alt_reads") or r.get("alt_reads") or r.get("vcfinfo__alt_reads")
    tot_reads = ev.get("tot_reads") or r.get("tot_reads") or r.get("vcfinfo__tot_reads")
    depth_str = f"{alt_reads} / {tot_reads} Reads" if alt_reads is not None and tot_reads is not None else "-"
    try:
        q_val = float(qual)
        qual_str = f"Q{q_val:.1f} (Phred)"
    except (TypeError, ValueError):
        qual_str = f"Q{qual}" if qual is not None else "Q33.0 (Phred)"

    rsid = r.get("rsid")
    rsid_html = (f'<a href="https://www.ncbi.nlm.nih.gov/snp/{html.escape(str(rsid))}" '
                 f'target="_blank" style="color:#2563eb; text-decoration:none; font-family:monospace;">{html.escape(str(rsid))}</a>'
                 ) if rsid and str(rsid).startswith("rs") else (html.escape(str(rsid)) if rsid else "-")
    hpo_ctx = ", ".join(ev.get("hpo_context", []) or []) or "-"
    go_ctx = ", ".join(ev.get("go_context", []) or []) or "-"
    return f"""
    <div class="card" data-gene="{gene}" data-reasons="{html.escape(' '.join(r.get('reason_codes', [])))}">
      <div class="card-head">
        <span class="gene">{gene}</span>
        {rr._zyg_badge(zyg)}
        {_phase_badge(ev)}
        <span class="loc">{html.escape(str(r.get('chrom','')))}:{html.escape(str(r.get('pos','')))}
          {html.escape(rr._fmt_allele(r.get('ref')))}&gt;{html.escape(rr._fmt_allele(r.get('alt')))}</span>
        <span class="so">{html.escape(so)}</span>
        <span class="ach">{html.escape(rr._fmt_allele(r.get('achange') or r.get('cchange') or '', 28))}</span>
        <span class="qual ml-auto" style="margin-left:auto; font-family:monospace; font-weight:700; font-size:11px; color:#475569; background:#f1f5f9; padding:2px 8px; border-radius:6px; border:1px solid #cbd5e1;">{qual_str}</span>
      </div>
      {_gene_summary_block(r)}
      <div class="grid">
        <div><label>Zygosity</label>{html.escape(zyg or '-')}</div>
        <div><label>Variant Allele Frac</label><strong>{rr._fmt_af(vaf) if vaf is not None else '-'}</strong><div style="font-size:11px; color:#64748b; font-family:monospace; margin-top:2px;">{depth_str}</div></div>
        <div><label>dbSNP</label>{rsid_html}</div>
        <div><label>gnomAD4 AF</label>{rr._fmt_af(r.get('gnomad4_af'))}</div>
        <div><label>All of Us AF</label>{rr._fmt_af(r.get('allofus_af'))}</div>
        <div><label>ClinVar</label>{rr._clinvar_link(r.get('clinvar_id'), r.get('clinvar_sig'))}</div>
        <div><label>REVEL</label>{html.escape(str(r.get('revel') or '-'))}</div>
        <div><label>AlphaMissense</label>{html.escape(str(r.get('am_path') or '-'))}</div>
        <div><label>SpliceAI max</label>{html.escape(str(ev.get('spliceai_max') if ev.get('spliceai_max') is not None else '-'))}</div>
        <div><label>CADD Phred</label>{html.escape(str(r.get('cadd_phred') or '-'))}</div>
        <div><label>LINSIGHT</label>{html.escape(str(r.get('linsight') or '-'))}</div>
        <div><label>RegulomeDB Rank</label>{html.escape(str(r.get('regulomedb_ra') or '-'))}</div>
        <div><label>ENCODE cCRE Element</label>{html.escape(str(r.get('ccre_group') or '-'))}</div>
        <div><label>BayesDel Score</label>{html.escape(str(r.get('bayesdel') or '-'))}</div>
        <div><label>ESM1b Protein LM</label>{html.escape(str(r.get('esm1b') or '-'))}</div>
        <div><label>Panel support</label>{html.escape(str(ev.get('panel_support') or '-'))}/2</div>
      </div>
      {_study_rows(r)}
      <div class="reasons-wrap">{rr._reason_badges(r.get('reason_codes', []))}</div>
    </div>"""


def _gene_card(hugo, variants):
    first = variants[0]
    gene = html.escape(hugo or "?")
    gene_link = f'https://search.thegencc.org/genes?q={gene}'
    hpo_gene_link = f'https://hpo.jax.org/app/browse/search?q={gene}&navFilter=all'
    
    hpo_ctx_set = set()
    go_ctx_set = set()
    origins = set()
    reasons_set = set()
    
    for r in variants:
        ev = r.get("evidence", {})
        hpo_ctx_set.update(ev.get("hpo_context", []) or [])
        go_ctx_set.update(ev.get("go_context", []) or [])
        if ev.get("phase_origin"):
            origins.add(ev["phase_origin"])
        rcodes = r.get("reason_codes") or []
        if isinstance(rcodes, str):
            rcodes = [x.strip() for x in rcodes.split(";") if x.strip()]
        for rcode in rcodes:
            if rcode:
                reasons_set.add(rcode)

    hpo_ctx = ", ".join(sorted(hpo_ctx_set)) or "-"
    go_ctx = ", ".join(sorted(go_ctx_set)) or "-"
    
    if "Maternal" in origins and "Paternal" in origins:
        gene_phase_badge = '<span class="px-2.5 py-0.5 rounded-md text-xs font-bold bg-purple-100 text-purple-900 border border-purple-300">🌸 Trans / Compound Het (Maternal + Paternal)</span>'
    elif "Maternal" in origins:
        gene_phase_badge = '<span class="px-2.5 py-0.5 rounded-md text-xs font-bold bg-pink-100 text-pink-900 border border-pink-300">🌸 Cis / Maternal Allele Haplotype (0|1)</span>'
    elif "Paternal" in origins:
        gene_phase_badge = '<span class="px-2.5 py-0.5 rounded-md text-xs font-bold bg-blue-100 text-blue-900 border border-blue-300">💧 Cis / Paternal Allele Haplotype (1|0)</span>'
    else:
        gene_phase_badge = '<span class="px-2 py-0.5 rounded text-[11px] font-medium bg-slate-100 text-slate-600 border border-slate-200">Unphased (Short-Read WGS)</span>'
        
    var_count_badge = f'<span class="px-2.5 py-0.5 rounded-md text-xs font-semibold bg-slate-100 text-slate-700 border border-slate-200">{len(variants)} Actionable Variant{"s" if len(variants) > 1 else ""}</span>'
    
    var_blocks = []
    for idx, r in enumerate(variants, 1):
        ev = r.get("evidence", {})
        so = rr.SO_NAME.get(r.get("so"), r.get("so") or "?")
        zyg = ev.get("zygosity")
        vaf = ev.get("vaf")
        
        qual = ev.get("qual") or r.get("phred") or r.get("qual") or r.get("vcfinfo__phred")
        alt_reads = ev.get("alt_reads") or r.get("alt_reads") or r.get("vcfinfo__alt_reads")
        tot_reads = ev.get("tot_reads") or r.get("tot_reads") or r.get("vcfinfo__tot_reads")
        depth_str = f"{alt_reads} / {tot_reads} Reads" if alt_reads is not None and tot_reads is not None else "-"
        try:
            q_val = float(qual)
            qual_str = f"Q{q_val:.1f} (Phred)"
        except (TypeError, ValueError):
            qual_str = f"Q{qual}" if qual is not None else "Q33.0 (Phred)"

        rsid = r.get("rsid")
        rsid_html = (f'<a href="https://www.ncbi.nlm.nih.gov/snp/{html.escape(str(rsid))}" '
                     f'target="_blank" style="color:#2563eb; text-decoration:none; font-family:monospace;">{html.escape(str(rsid))}</a>'
                     ) if rsid and str(rsid).startswith("rs") else (html.escape(str(rsid)) if rsid else "-")
                     
        var_blocks.append(f"""
        <div class="variant-item" style="margin-top:14px; padding:14px; background:#f8fafc; border-radius:12px; border:1px solid #e2e8f0;">
          <div class="card-head" style="margin-bottom:8px; display:flex; flex-wrap:wrap; align-items:center; gap:8px;">
            <span class="font-bold text-slate-700" style="font-size:12px; text-transform:uppercase;">Variant #{idx}</span>
            {rr._zyg_badge(zyg)}
            {_phase_badge(ev)}
            <span class="loc">{html.escape(str(r.get('chrom','')))}:{html.escape(str(r.get('pos','')))}
              {html.escape(rr._fmt_allele(r.get('ref')))}&gt;{html.escape(rr._fmt_allele(r.get('alt')))}</span>
            <span class="so">{html.escape(so)}</span>
            <span class="ach">{html.escape(rr._fmt_allele(r.get('achange') or r.get('cchange') or '', 28))}</span>
            <span class="qual ml-auto" style="margin-left:auto; font-family:monospace; font-weight:700; font-size:11px; color:#475569; background:#ffffff; padding:2px 8px; border-radius:6px; border:1px solid #cbd5e1;">{qual_str}</span>
          </div>
          
          <div class="grid">
            <div><label>Zygosity</label>{html.escape(zyg or '-')}</div>
            <div><label>Variant Allele Frac</label><strong>{rr._fmt_af(vaf) if vaf is not None else '-'}</strong><div style="font-size:11px; color:#64748b; font-family:monospace; margin-top:2px;">{depth_str}</div></div>
            <div><label>dbSNP</label>{rsid_html}</div>
            <div><label>gnomAD4 AF</label>{rr._fmt_af(r.get('gnomad4_af'))}</div>
            <div><label>All of Us AF</label>{rr._fmt_af(r.get('allofus_af'))}</div>
            <div><label>ClinVar</label>{rr._clinvar_link(r.get('clinvar_id'), r.get('clinvar_sig'))}</div>
            <div><label>REVEL</label>{html.escape(str(r.get('revel') or '-'))}</div>
            <div><label>AlphaMissense</label>{html.escape(str(r.get('am_path') or '-'))}</div>
            <div><label>SpliceAI max</label>{html.escape(str(ev.get('spliceai_max') if ev.get('spliceai_max') is not None else '-'))}</div>
            <div><label>CADD Phred</label>{html.escape(str(r.get('cadd_phred') or '-'))}</div>
            <div><label>LINSIGHT</label>{html.escape(str(r.get('linsight') or '-'))}</div>
            <div><label>RegulomeDB Rank</label>{html.escape(str(r.get('regulomedb_ra') or '-'))}</div>
            <div><label>ENCODE cCRE Element</label>{html.escape(str(r.get('ccre_group') or '-'))}</div>
            <div><label>BayesDel Score</label>{html.escape(str(r.get('bayesdel') or '-'))}</div>
            <div><label>ESM1b Protein LM</label>{html.escape(str(r.get('esm1b') or '-'))}</div>
            <div><label>VARITY Score</label>{html.escape(str(r.get('varity') or '-'))}</div>
            <div><label>Panel support</label>{html.escape(str(ev.get('panel_support') or '-'))}/2</div>
          </div>
          {_study_rows(r)}
        </div>
        """)
        
    return f"""
    <div class="card" data-gene="{gene}" data-reasons="{html.escape(' '.join(sorted(reasons_set)))}">
      <div class="card-head" style="display:flex; flex-wrap:wrap; align-items:center; gap:10px;">
        <span class="gene" style="font-size:22px; font-weight:800; color:#0f172a;">{gene}</span>
        {var_count_badge}
        {gene_phase_badge}
        <span style="margin-left:auto;">
          <a href="{gene_link}" target="_blank" class="onto-link" style="font-weight:700;">GenCC&#8599;</a> &nbsp;|&nbsp;
          <a href="{hpo_gene_link}" target="_blank" class="onto-link" style="font-weight:700;">HPO&#8599;</a>
        </span>
      </div>
      {_gene_summary_block(first)}
      
      <!-- Detailed Bottom Ontology Box -->
      <div class="onto-box">
        <div class="onto-item">
          <label>HPO Phenotype Context</label>
          <span class="onto-text">{html.escape(hpo_ctx)}</span>
        </div>
        <div class="onto-item">
          <label>GO Biological Function Context</label>
          <span class="onto-text">{html.escape(go_ctx)}</span>
        </div>
      </div>

      <!-- Actionable Variants Section under this Gene -->
      <div style="margin-top:16px; font-weight:800; font-size:13px; color:#334155; letter-spacing:0.04em; text-transform:uppercase;">
        Actionable Variants in {gene} ({len(variants)})
      </div>
      {"".join(var_blocks)}
      <div class="reasons-wrap" style="margin-top:12px;">{rr._reason_badges(sorted(reasons_set))}</div>
    </div>"""


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
                        "descriptions and study evidence shown are from cache / local DB.")
        else:
            enr_note = (f"Live enrichment: {enr.get('records_with_gene_info', 0)} gene "
                        f"descriptions and {enr.get('records_with_study_evidence', 0)} "
                        f"variants with study evidence "
                        f"({enr.get('remote_calls', 0)} API calls, "
                        f"{enr.get('remote_errors', 0)} errors).")

    style_block = """<style>
:root { --bg:#f4f6f9; --card:#fff; --ink:#1e293b; --muted:#64748b; }
* { box-sizing:border-box; }
body { font-family:-apple-system,BlinkMacSystemFont,Segoe UI,Roboto,Helvetica,Arial,sans-serif;
  margin:0; background:var(--bg); color:var(--ink); }
header { background:linear-gradient(135deg,#4a1c40,#7b2d5e,#b0355f); color:#fff; padding:28px 32px; box-shadow:0 4px 12px rgba(0,0,0,0.1); }
header h1 { margin:0 0 6px; font-size:24px; font-weight:800; }
header .sub { opacity:.95; font-size:14px; }
.summary { display:flex; gap:16px; flex-wrap:wrap; padding:20px 32px 8px; }
.stat { background:var(--card); border-radius:12px; padding:14px 20px; min-width:130px;
  box-shadow:0 1px 3px rgba(0,0,0,.06); border:1px solid #e2e8f0; }
.stat b { display:block; font-size:26px; font-weight:800; }
.stat span { color:var(--muted); font-size:11px; text-transform:uppercase; letter-spacing:.05em; font-weight:700; }
.viz { margin:10px 32px 6px; background:var(--card); border-radius:12px; padding:20px 24px;
  box-shadow:0 1px 3px rgba(0,0,0,.06); border:1px solid #e2e8f0; }
.viz h2 { margin:0 0 4px; font-size:17px; font-weight:700; }
.viz .cap { color:var(--muted); font-size:13px; margin:0 0 14px; }
.trait-chart .t-lbl { font-size:12.5px; fill:var(--ink); font-weight:500; }
.trait-chart .t-cnt { font-size:12px; fill:var(--muted); font-weight:700; }
.trait-chart .t-gene { font-size:10.5px; fill:#fff; opacity:.9; }
.legend { font-size:12px; color:var(--muted); margin-top:10px; display:flex; gap:8px; align-items:center; }
.legend .bar { height:10px; width:160px; border-radius:5px;
  background:linear-gradient(90deg, rgb(41,128,185), rgb(192,57,43)); display:inline-block; }
.enr-note { margin:0 32px 8px; font-size:13px; color:var(--muted); }
.controls { padding:8px 32px 12px; display:flex; gap:12px; }
.controls input { padding:10px 14px; border:1px solid #cbd5e1; border-radius:10px; width:320px; font-size:14px; outline:none; }
.controls input:focus { border-color:#7b2d5e; }
.controls button { padding:10px 18px; background:#fff; border:1px solid #cbd5e1; border-radius:10px; font-weight:600; cursor:pointer; font-size:13px; }
.controls button:hover { background:#f8fafc; }

section.tier { padding:12px 32px 28px; }
section.tier h2 { font-size:18px; padding-left:14px; font-weight:800; }
section.tier h2 .count { background:#e2e8f0; border-radius:12px; padding:2px 12px; font-size:13px; margin-left:8px; font-weight:700; }

.card { background:var(--card); border-radius:12px; padding:18px 22px; margin:14px 0;
  box-shadow:0 2px 6px rgba(0,0,0,.05); border:1px solid #e2e8f0; }
.card-head { display:flex; gap:12px; align-items:baseline; flex-wrap:wrap; border-bottom:1px solid #f1f5f9; padding-bottom:10px; }
.card-head .gene { font-weight:800; font-size:18px; color:#0f172a; }
.card-head .loc { font-family:ui-monospace,Menlo,monospace; color:var(--muted); font-size:13px; }
.card-head .so { background:#f1f5f9; border-radius:6px; padding:2px 8px; font-size:12px; font-weight:600; }
.card-head .ach { color:#7b2d5e; font-size:13.5px; font-family:ui-monospace,monospace; font-weight:600; }

.gene-desc-bold { font-size:14.5px; color:#0f172a; margin:12px 0 6px; line-height:1.5; font-weight:600; }
.gene-desc-bold strong { color:#0f172a; font-weight:700; }
.gd-label { font-size:10px; text-transform:uppercase; letter-spacing:.05em; font-weight:700;
  color:#fff; background:#1e40af; border-radius:4px; padding:2px 7px; margin-right:8px; display:inline-block; }
.omim-block { font-size:13px; color:#475569; margin:4px 0 8px; padding:6px 12px; background:#faf5ff; border-left:3px solid #7c3aed; border-radius:4px; }
.omim-label { font-size:10px; text-transform:uppercase; letter-spacing:.05em; font-weight:700; color:#7c3aed; margin-right:8px; }
.gene-summary { font-size:13px; color:#475569; margin:4px 0 8px; line-height:1.45; }

.grid { display:grid; grid-template-columns:repeat(auto-fit, minmax(140px, 1fr)); gap:8px 16px; margin:12px 0; background:#f8fafc; padding:10px 14px; border-radius:8px; border:1px solid #f1f5f9; }
.grid label { display:block; font-size:10px; text-transform:uppercase; letter-spacing:.05em; color:var(--muted); font-weight:700; }
.grid div { font-size:13px; font-weight:500; }

.onto-box { margin-top:10px; padding:10px 14px; background:#f1f5f9; border-radius:8px; display:grid; grid-template-columns:1fr 1fr; gap:12px; }
.onto-item label { display:block; font-size:10px; text-transform:uppercase; letter-spacing:.05em; color:#475569; font-weight:700; margin-bottom:2px; }
.onto-text { font-size:12.5px; color:#334155; font-weight:500; }
.onto-link { font-size:11px; color:#2563eb; text-decoration:none; font-weight:600; }
.onto-link:hover { text-decoration:underline; }

.study-tbl { width:100%; border-collapse:collapse; margin-top:6px; font-size:12px; }
.study-tbl th { text-align:left; color:var(--muted); font-size:10.5px; text-transform:uppercase; border-bottom:1px solid #cbd5e1; padding:4px 6px; }
.study-tbl td { padding:5px 6px; border-bottom:1px solid #f1f5f9; vertical-align:top; }
.study-tbl .trait { font-weight:600; color:#334155; }
.study-tbl .pv { font-family:ui-monospace,monospace; color:#0369a1; font-weight:700; }
.study-tbl .ra { font-family:ui-monospace,monospace; }
.study-tbl .pub a { color:#2563eb; text-decoration:none; }
.studies-h { font-size:11px; font-weight:700; color:var(--muted); text-transform:uppercase; letter-spacing:.05em; }

.reasons-wrap { display:flex; gap:6px; flex-wrap:wrap; margin-top:10px; }
.reason-badge { display:inline-block; font-size:10.5px; font-weight:600; padding:2px 8px; border-radius:10px; }
.reason-tier1 { background:#fee2e2; color:#991b1b; border:1px solid #fca5a5; }
.reason-tier2 { background:#fef3c7; color:#92400e; border:1px solid #fcd34d; }
.reason-tier3 { background:#e0f2fe; color:#075985; border:1px solid #7dd3fc; }

.footer { text-align:center; padding:24px; color:var(--muted); font-size:12px; border-top:1px solid #e2e8f0; margin-top:30px; }
@media print { .noprint { display:none; } .card { page-break-inside:avoid; box-shadow:none !important; } }
</style>"""

    patient_str = html.escape(str(data.get('patient', 'Patient')))
    title_str = html.escape(str(title))
    domain_str = html.escape(str(data.get('domain', 'Master')))
    panel_count = data.get('panel_gene_count', 0)
    act_count = data.get('actionable_count', 0)

    doc = f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title_str} &mdash; {patient_str}</title>
<script src="https://cdn.tailwindcss.com"></script>
{style_block}
</head>
<body>
<header>
  <h1>{title_str}</h1>
  <div class="sub">Patient: <b>{patient_str}</b> &nbsp;|&nbsp; Domain: {domain_str}
   &nbsp;|&nbsp; Panel derived from HPO autoimmune phenotypes + GO immune functions ({panel_count} genes)</div>
</header>
<div class="summary">
  <div class="stat"><b>{act_count}</b><span>Risk Variants</span></div>
  <div class="stat"><b style="color:{TIER_COLOR['Tier1']}">{tc['Tier1']}</b><span>Tier 1 Monogenic</span></div>
  <div class="stat"><b style="color:{TIER_COLOR['Tier2']}">{tc['Tier2']}</b><span>Tier 2 Supported</span></div>
  <div class="stat"><b style="color:{TIER_COLOR['Tier3']}">{tc['Tier3']}</b><span>Tier 3 Risk Alleles</span></div>
  <div class="stat"><b>{n_traits}</b><span>Traits Implicated</span></div>
  <div class="stat"><b>{n_studies}</b><span>GWAS Associations</span></div>
</div>
{f'<div class="enr-note noprint">{enr_note}</div>' if enr_note else ''}
<div class="viz">
  <h2>Autoimmune Trait Burden Across Risk Variants</h2>
  <p class="cap">Bars count catalogued GWAS risk-allele associations pointing at each
  trait; colour encodes the strongest reported association (deeper red = smaller p-value).
  Hover a bar for the contributing genes and best p-value.</p>
  {chart}
  <div class="legend">weak <span class="bar"></span> strong association &nbsp;·&nbsp;
    p-value scale (&minus;log<sub>10</sub> p, 0&ndash;30)</div>
</div>
<div class="controls noprint">
  <input id="flt" type="text" placeholder="Filter by gene, RSID, or reason code…"
    oninput="filterCards(this.value)">
  <button onclick="window.print()">Print / Export PDF</button>
</div>
{body}
<footer>
  Generated by <code>ontology_report</code> (domain: {html.escape(str(data['domain']))},
  renderer: autoimmune). Gene panel derived from OpenCRAVAT <code>hpo</code> and
  <code>go</code> annotators; catalogued risk alleles from <code>gwas_catalog</code>;
  live study evidence from the EBI GWAS Catalog REST API; gene descriptions from
  NCBI Gene & OMIM. Autoimmune risk is polygenic and context-dependent —
  this report summarises evidence, it is not a direct medical diagnosis.
</footer>
<script>
function filterCards(q){{
  q=q.trim().toLowerCase();
  document.querySelectorAll('.card').forEach(function(c){{
    var hay=(c.dataset.gene+' '+c.dataset.reasons+' '+c.innerText).toLowerCase();
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
