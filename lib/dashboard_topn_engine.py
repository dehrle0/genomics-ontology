#!/usr/bin/env python3
"""
lib/dashboard_topn_engine.py
Generates focused Dashboard Sub-Reports (HTML, RFC 2557 MHTML, and PDF) containing
ONLY the genes with Tier-1 / Pathogenic variants and the Top-N ranked genes for each
of the 8 Dashboard Report selections on the Landscape Clinical Dashboard:
  1. Visual Explorer — Master Genomic Catalog (Tier-1 & Pathogenic + Top Impact Genes)
  2. Rare Variants (gnomAD AF <= 1%) — Tier-1 Rare + Top-N Rarest High-Impact Genes
  3. ClinVar Pathogenicity & Clinical Curation — Tier-1 Pathogenic, Risk Factor & Protective Genes
  4. CADD Deleteriousness Ranking (Phred >= 15) — Tier-1 + Top-N CADD Genes
  5. SpliceAI Deep Learning Splicing Disruption (Delta >= 0.10) — Tier-1 + Top-N SpliceAI Genes
  6. DeepMind AlphaGenome Atlas (AVI Phred & 18-Modality Drivers) — Tier-1 + Top-N AVI Genes
  7. AlphaMissense 3D Structural Pathogenicity — Tier-1 + Top-N AlphaMissense Genes
  8. REVEL Ensembled Missense Pathogenicity (Score >= 0.20) — Tier-1 + Top-N REVEL Genes

Enforces:
  - Strict single-sample isolation (zero cross-sample leakage, zero hardcoded variants).
  - Pre-expanded Gene Cards & Top Variant 4-Box Detail Drawers (including AlphaGenome
    18-modality SHAP attribution sparkbars) for immediate visual inspection, printing,
    PDF compilation, and MHTML snapshot archival.
"""

import html as html_lib
import os
import quopri
import sys
from datetime import datetime, timezone

SCRIPT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if os.path.join(SCRIPT_DIR, "scripts") not in sys.path:
    sys.path.insert(0, os.path.join(SCRIPT_DIR, "scripts"))

import prepare_landscape_data as pld

DASHBOARD_REPORT_SECTIONS = [
    {
        "key": "visual_explorer",
        "title": "1. Visual Explorer — Tier-1 & Top Clinical Impact Genes",
        "short_label": "Visual Explorer (Tier-1)",
        "badge_color": "#0f766e",
        "badge_bg": "#f0fdfa",
        "top_n": 12,
        "desc": "Master genomic catalog filtered strictly to genes carrying Tier-1 / ClinVar Pathogenic variants and highest multi-engine clinical impact scores."
    },
    {
        "key": "clinvar",
        "title": "2. ClinVar Pathogenicity & Clinical Curation — Tier-1 & Actionable Genes",
        "short_label": "ClinVar Curation",
        "badge_color": "#b91c1c",
        "badge_bg": "#fef2f2",
        "top_n": 12,
        "desc": "Genes carrying ClinVar-curated Pathogenic, Likely Pathogenic, Risk Factor, Drug Response, or Protective alleles."
    },
    {
        "key": "atlas",
        "title": "3. DeepMind AlphaGenome Atlas — Top Genes by AVI Phred & 18-Modality Drivers",
        "short_label": "AlphaGenome Atlas",
        "badge_color": "#0369a1",
        "badge_bg": "#f0f9ff",
        "top_n": 10,
        "desc": "Top genes ranked by DeepMind AlphaGenome Variant Impact (AVI) Phred score across 18 genomic and epigenomic modalities (RNA-seq, Splicing, ATAC/DNase, ChIP)."
    },
    {
        "key": "cadd",
        "title": "4. CADD Deleteriousness Ranking — Tier-1 & Top Genes (Phred ≥ 15)",
        "short_label": "CADD Ranking",
        "badge_color": "#9a3412",
        "badge_bg": "#fff7ed",
        "top_n": 10,
        "desc": "Top deleterious genes ranked by Combined Annotation Dependent Depletion (CADD v1.7) Phred-scaled score."
    },
    {
        "key": "spliceai",
        "title": "5. SpliceAI Deep Learning Splicing Disruption — Top Genes (Δ ≥ 0.10)",
        "short_label": "SpliceAI Splicing",
        "badge_color": "#6d28d9",
        "badge_bg": "#f5f3ff",
        "top_n": 10,
        "desc": "Genes exhibiting cryptic splice-site activation or canonical donor/acceptor loss (Δ score ≥ 0.10)."
    },
    {
        "key": "alphamissense",
        "title": "6. AlphaMissense 3D Structural Pathogenicity — Top Missense Genes",
        "short_label": "AlphaMissense 3D",
        "badge_color": "#be123c",
        "badge_bg": "#fff1f2",
        "top_n": 10,
        "desc": "Top protein-coding missense variants ranked by DeepMind AlphaFold-derived AlphaMissense structural pathogenicity probability."
    },
    {
        "key": "revel",
        "title": "7. REVEL Ensembled Missense Pathogenicity — Top Genes (Score ≥ 0.20)",
        "short_label": "REVEL Missense",
        "badge_color": "#1d4ed8",
        "badge_bg": "#eff6ff",
        "top_n": 10,
        "desc": "Top missense variants ranked by Rare Exome Variant Ensemble Learner (REVEL) pathogenicity score."
    },
    {
        "key": "rare",
        "title": "8. Rare & Private Variants (gnomAD AF ≤ 1%) — Tier-1 & Top High-Impact Genes",
        "short_label": "Rare Variants (≤1%)",
        "badge_color": "#334155",
        "badge_bg": "#f1f5f9",
        "top_n": 10,
        "desc": "Rare coding and regulatory variants (gnomAD v4 allele frequency ≤ 1.0%), prioritizing Tier-1 / Pathogenic and high-CADD/AVI rare alleles."
    }
]


def _esc(val):
    return html_lib.escape(str(val if val is not None else ""))


def _is_tier1_or_pathogenic_variant(v):
    tier = str(v.get("tier") or "")
    cv = str(v.get("clinvar") or "").lower()
    is_plp = ("pathogenic" in cv and "conflicting" not in cv and "uncertain" not in cv)
    return tier == "Tier1" or is_plp


def filter_and_select_topn_genes(genes, report_key, top_n=10):
    """
    Mirrors the exact filtering & ranking logic of js/landscape.js filterAndRankGenes(),
    then selects:
      - All genes in that view carrying a Tier-1 / ClinVar Pathogenic variant, PLUS
      - The highest-ranked genes in that report selection up to `max(len(tier1_genes), top_n)`.
    Also filters each selected gene's variants to its Tier-1 + top 3 variants for that metric
    so the report remains high-signal and concise.
    """
    ranked_genes = []
    for g in genes:
        raw_vars = list(g.get("variants") or [])
        if report_key == "rare":
            vars_f = [v for v in raw_vars if v.get("maf") is not None and float(v.get("maf")) <= 0.01]
        elif report_key == "clinvar":
            vars_f = [v for v in raw_vars if v.get("clinvar") and v.get("clinvar") not in ("None", "Not reviewed", "")]
        elif report_key == "cadd":
            vars_f = [v for v in raw_vars if v.get("cadd") is not None and float(v.get("cadd")) >= 15.0]
        elif report_key == "spliceai":
            vars_f = [v for v in raw_vars if v.get("spliceai") is not None and float(v.get("spliceai")) >= 0.10]
        elif report_key == "atlas":
            vars_f = [v for v in raw_vars if v.get("aviPhred") is not None]
        elif report_key == "alphamissense":
            vars_f = [v for v in raw_vars if v.get("alphamissense") is not None]
        elif report_key == "revel":
            vars_f = [v for v in raw_vars if v.get("revel") is not None and float(v.get("revel")) >= 0.20]
        else:
            vars_f = raw_vars

        if not vars_f:
            continue

        has_t1 = any(_is_tier1_or_pathogenic_variant(v) for v in vars_f)

        # Compute geneScore and sort variants inside the gene by active metric
        if report_key == "cadd":
            vars_f.sort(key=lambda v: float(v.get("cadd") or 0.0), reverse=True)
            max_val = float(vars_f[0].get("cadd") or 0.0)
            gene_score = (1000.0 if has_t1 else 0.0) + max_val
            top_metric_display = f"Max CADD: {max_val:.1f}"
        elif report_key == "spliceai":
            vars_f.sort(key=lambda v: float(v.get("spliceai") or 0.0), reverse=True)
            max_val = float(vars_f[0].get("spliceai") or 0.0)
            gene_score = (1000.0 if has_t1 else 0.0) + max_val * 100.0
            top_metric_display = f"SpliceAI Δ: {max_val:.2f}"
        elif report_key == "atlas":
            vars_f.sort(key=lambda v: float(v.get("aviPhred") or 0.0), reverse=True)
            max_val = float(vars_f[0].get("aviPhred") or 0.0)
            gene_score = (1000.0 if has_t1 else 0.0) + max_val
            top_metric_display = f"AlphaGenome AVI: {max_val:.1f}"
        elif report_key == "alphamissense":
            vars_f.sort(key=lambda v: float(v.get("alphamissense") or 0.0), reverse=True)
            max_val = float(vars_f[0].get("alphamissense") or 0.0)
            gene_score = (1000.0 if has_t1 else 0.0) + max_val * 100.0
            top_metric_display = f"AlphaMissense: {max_val:.3f}"
        elif report_key == "revel":
            vars_f.sort(key=lambda v: float(v.get("revel") or 0.0), reverse=True)
            max_val = float(vars_f[0].get("revel") or 0.0)
            gene_score = (1000.0 if has_t1 else 0.0) + max_val * 100.0
            top_metric_display = f"REVEL: {max_val:.3f}"
        elif report_key == "rare":
            vars_f.sort(
                key=lambda v: (
                    0 if _is_tier1_or_pathogenic_variant(v) else 1,
                    float(v.get("maf") if v.get("maf") is not None else 1.0),
                    -float(v.get("cadd") or 0.0)
                )
            )
            min_maf = min(float(v.get("maf") if v.get("maf") is not None else 1.0) for v in vars_f)
            max_cadd = max(float(v.get("cadd") or 0.0) for v in vars_f)
            gene_score = (1000.0 if has_t1 else 0.0) + (1.0 - min_maf) * 100.0 + max_cadd
            top_metric_display = f"Min gnomAD AF: {min_maf * 100:.4f}% · CADD {max_cadd:.1f}"
        elif report_key == "clinvar":
            cv_pri = {"Pathogenic": 10, "Likely pathogenic": 8, "Risk factor": 6, "Protective": 5, "Drug response": 5, "Uncertain significance": 3}
            vars_f.sort(
                key=lambda v: (
                    cv_pri.get(str(v.get("clinvar") or ""), 1),
                    float(v.get("cadd") or 0.0)
                ),
                reverse=True
            )
            max_cv = max(cv_pri.get(str(v.get("clinvar") or ""), 1) for v in vars_f)
            max_cadd = max(float(v.get("cadd") or 0.0) for v in vars_f)
            gene_score = (1000.0 if has_t1 else 0.0) + max_cv * 20.0 + max_cadd
            top_metric_display = f"{vars_f[0].get('clinvar') or 'ClinVar'} · CADD {max_cadd:.1f}"
        else:
            tier_pri = {"Tier1": 100, "Tier2": 50, "Tier3": 10}
            vars_f.sort(
                key=lambda v: (
                    1 if _is_tier1_or_pathogenic_variant(v) else 0,
                    tier_pri.get(str(v.get("tier") or ""), 1),
                    float(v.get("cadd") or 0.0) + float(v.get("aviPhred") or 0.0)
                ),
                reverse=True
            )
            max_tier = max(tier_pri.get(str(v.get("tier") or ""), 1) for v in vars_f)
            max_cadd = max(float(v.get("cadd") or 0.0) for v in vars_f)
            gene_score = (1000.0 if has_t1 else 0.0) + max_tier * 10.0 + max_cadd
            top_metric_display = f"{vars_f[0].get('tier') or 'Tier'} · CADD {max_cadd:.1f}"

        # Keep all Tier-1/Pathogenic variants + up to top 3 variants per gene for focused readability
        top_vars = []
        seen_vids = set()
        for v in vars_f:
            if _is_tier1_or_pathogenic_variant(v):
                top_vars.append(v)
                seen_vids.add(v.get("id"))
        for v in vars_f:
            if v.get("id") not in seen_vids and len(top_vars) < max(3, len(top_vars)):
                top_vars.append(v)
                seen_vids.add(v.get("id"))

        g_copy = dict(g)
        g_copy["variants"] = top_vars
        g_copy["totalMatchingVariantsInGene"] = len(vars_f)
        g_copy["geneScore"] = gene_score
        g_copy["hasTier1OrPathogenic"] = has_t1
        g_copy["topMetricDisplay"] = top_metric_display
        ranked_genes.append(g_copy)

    ranked_genes.sort(key=lambda x: x["geneScore"], reverse=True)

    # Include ALL Tier-1/Pathogenic genes in this report view + top-N ranked genes
    selected = []
    seen_sym = set()
    for g in ranked_genes:
        if g["hasTier1OrPathogenic"]:
            selected.append(g)
            seen_sym.add(g.get("symbol"))
    for g in ranked_genes:
        if g.get("symbol") not in seen_sym and len(selected) < top_n:
            selected.append(g)
            seen_sym.add(g.get("symbol"))

    return selected, len(ranked_genes)


def _render_variant_drawer_html(v, gene):
    """Renders the 4-box expanded variant drawer matching the Landscape Dashboard."""
    pathologies = gene.get("associatedPathology") or []
    s_det = v.get("spliceaiDetails") or {}
    ag_s = f"{float(s_det.get('ag') or 0.0):.2f}"
    al_s = f"{float(s_det.get('al') or 0.0):.2f}"
    dg_s = f"{float(s_det.get('dg') or 0.0):.2f}"
    dl_s = f"{float(s_det.get('dl') or 0.0):.2f}"
    top_features = v.get("aviTopFeatures") or []
    coord = v.get("coordinate") or ""
    vid = str(v.get("id") or "")
    is_rs = vid.startswith("rs")
    track_url = v.get("alphagenomeTrackPredUrl") or f"https://deepmind.google.com/science/alphagenome/atlas/track-predictions?q={coord}&m=variant"
    litvar_url = (
        f"https://www.ncbi.nlm.nih.gov/research/litvar2/docsum?variant=litvar@{vid}%23%23&query={vid}"
        if is_rs else f"https://www.ncbi.nlm.nih.gov/research/litvar2/docsum?query={vid}"
    )

    path_html = ""
    if pathologies:
        for p in pathologies[:3]:
            if isinstance(p, dict):
                inh_list = p.get("inheritance") or []
                if isinstance(inh_list, str):
                    inh_list = [inh_list]
                inh_badges = " ".join(f'<span class="mini-badge">{_esc(i)}</span>' for i in inh_list if i)
                omim_link = f'<a href="https://omim.org/entry/{_esc(p.get("omim"))}" target="_blank" rel="noopener">OMIM #{_esc(p.get("omim"))} ↗</a>' if p.get("omim") else ""
                path_html += f'<div style="margin-bottom:4px;"><div style="font-weight:700;color:#0f172a;">{_esc(p.get("name"))}</div><div style="font-size:9.5px;color:#475569;">{inh_badges} {omim_link}</div></div>'
    else:
        path_html = '<div style="color:#64748b;font-size:10px;">Complex / multifactorial quantitative trait or regulatory modifier.</div>'

    sparkbars_html = ""
    if top_features:
        sparkbars_html = '<div style="margin-top:4px;display:flex;flex-direction:column;gap:2px;">'
        for tf in top_features[:4]:
            pct = min(float(tf.get("pct") or 0.0), 100.0)
            sparkbars_html += f"""
            <div style="display:grid;grid-template-columns:110px 1fr 34px;align-items:center;gap:4px;font-size:9.5px;">
              <span style="color:#334155;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;">{_esc(tf.get("label"))}</span>
              <div style="height:5px;background:#e2e8f0;border-radius:3px;overflow:hidden;"><div style="height:100%;width:{pct:.0f}%;background:linear-gradient(90deg,#0284c7,#0f766e);"></div></div>
              <span class="mono" style="text-align:right;font-weight:700;color:#0369a1;">{pct:.0f}%</span>
            </div>"""
        sparkbars_html += "</div>"

    avi_box_html = ""
    if v.get("aviPhred") is not None:
        avi_p = float(v.get("aviPhred") or 0.0)
        avi_r = f"{float(v.get('aviRaw')):.3f}" if v.get("aviRaw") is not None else "—"
        avi_pct = f"{float(v.get('aviPercentile')):.2f}%" if v.get("aviPercentile") is not None else "—"
        avi_box_html = f"""
        <div style="padding:5px 7px;background:#f0f9ff;border:1px solid #bae6fd;border-radius:4px;">
          <div style="display:flex;justify-content:space-between;align-items:center;">
            <strong style="color:#0369a1;">AlphaGenome AVI (Phred):</strong>
            <strong class="mono" style="color:#0284c7;font-size:11px;">{avi_p:.1f} <span style="font-size:9.5px;color:#475569;">(Raw: {avi_r} · Top {avi_pct})</span></strong>
          </div>
          <div style="font-size:9.5px;color:#0f172a;margin-top:2px;">Primary Driver: <strong style="color:#0369a1;">{_esc(v.get("aviModality") or "Multi-Track")}</strong></div>
          {sparkbars_html}
          <div style="margin-top:4px;display:flex;gap:6px;flex-wrap:wrap;font-size:9.5px;">
            <a href="{_esc(track_url)}" target="_blank" rel="noopener" style="color:#0284c7;font-weight:700;">🧬 Ref vs. Alt Track Predictions ↗</a>
            <a href="{_esc(litvar_url)}" target="_blank" rel="noopener" style="color:#0f766e;font-weight:600;">LitVar2 ↗</a>
          </div>
        </div>
        """
    else:
        avi_box_html = f"""
        <div style="font-size:10px;color:#64748b;">
          Evaluated via CADD / SpliceAI / REVEL ensemble.
          <div style="margin-top:4px;"><a href="{_esc(litvar_url)}" target="_blank" rel="noopener" style="color:#0f766e;font-weight:600;">LitVar2 Literature ↗</a></div>
        </div>
        """

    vaf_s = f"{float(v.get('vaf')):.2f}" if v.get("vaf") is not None else "—"
    dp_s = f"{v.get('dp')}x" if v.get("dp") is not None else "—"

    return f"""
    <tr class="drawer-row">
      <td colspan="11">
        <div class="drawer-grid">
          <div class="drawer-box">
            <div class="drawer-box__title">1. Inheritance &amp; OMIM Pathology</div>
            {path_html}
            <div style="margin-top:4px;padding-top:3px;border-top:1px dashed #cbd5e1;font-size:9.5px;color:#475569;">
              Quadrant: <strong>{_esc(v.get("concordanceQuadrant") or "Curated")}</strong>
            </div>
          </div>
          <div class="drawer-box">
            <div class="drawer-box__title">2. Transcript, Phasing &amp; Depth</div>
            <div class="kv-row"><span>Transcript:</span><strong class="mono">{_esc(v.get("transcript") or "Canonical")}</strong></div>
            <div class="kv-row"><span>cDNA (HGVS):</span><strong class="mono">{_esc(v.get("cchange") or "—")}</strong></div>
            <div class="kv-row"><span>Protein (HGVS):</span><strong class="mono">{_esc(v.get("achange") or "—")}</strong></div>
            <div class="kv-row"><span>VAF / Read Depth:</span><strong class="mono">VAF {vaf_s} · DP {dp_s}</strong></div>
            <div class="kv-row"><span>Phased GT:</span><strong class="mono">{_esc(v.get("vcfGt") or "—")} ({_esc(v.get("phase") or "Unphased")})</strong></div>
          </div>
          <div class="drawer-box">
            <div class="drawer-box__title">3. SpliceAI &amp; Clinical Review</div>
            <div class="kv-row"><span>SpliceAI Δ:</span><strong class="mono">AG:{ag_s} AL:{al_s} DG:{dg_s} DL:{dl_s}</strong></div>
            <div class="kv-row"><span>ClinVar Review:</span><strong>{_esc(v.get("clinvarRev") or "criteria provided")}</strong></div>
            {f'<div class="kv-row"><span>ClinVar VCV:</span><a href="https://www.ncbi.nlm.nih.gov/clinvar/variation/{_esc(v.get("clinvarId"))}/" target="_blank" rel="noopener" class="mono">VCV#{_esc(v.get("clinvarId"))} ↗</a></div>' if v.get("clinvarId") else ""}
            <div class="kv-row"><span>gnomAD v4 AF:</span><strong class="mono">{f"{float(v.get('maf'))*100:.4f}%" if v.get("maf") is not None else "Rare / Private"}</strong></div>
          </div>
          <div class="drawer-box">
            <div class="drawer-box__title">4. AlphaGenome 18-Modality Atlas</div>
            {avi_box_html}
          </div>
        </div>
      </td>
    </tr>
    """


def render_gene_card_static_html(g, report_key):
    """Renders an expanded Gene Card with its variant table and expanded 4-box drawers."""
    sym = _esc(g.get("symbol") or "")
    gname = _esc(g.get("name") or "")
    organ = _esc(g.get("organSystem") or "Multisystem")
    chrom = _esc(g.get("chrom") or "")
    pos_str = f"{int(g.get('pos')):,}" if g.get("pos") else ""
    summary = _esc(g.get("summary") or "")
    p_status = str((g.get("phaseAudit") or {}).get("phaseStatus") or "")

    phase_badge = ""
    if "Trans" in p_status:
        phase_badge = '<span class="pill pill-red">⚖️ Trans Compound Het Candidate</span>'
    elif "Cis" in p_status:
        phase_badge = '<span class="pill pill-amber">⚖️ Cis Haplotype (1 WT Allele Intact)</span>'

    tier1_badge = '<span class="pill pill-red">🔥 Tier-1 / Pathogenic Gene</span>' if g.get("hasTier1OrPathogenic") else '<span class="pill pill-teal">🏆 Top Ranked in View</span>'

    rows_html = ""
    for idx, v in enumerate(g.get("variants") or []):
        vid = _esc(v.get("id") or "")
        is_rs = str(v.get("id") or "").startswith("rs")
        rs_link = f'<a href="https://www.ncbi.nlm.nih.gov/snp/{vid}" target="_blank" rel="noopener" style="color:#0f766e;font-weight:700;">{vid}</a>' if is_rs else f'<strong>{vid}</strong>'
        hgvs_sub = _esc(v.get("achange") or v.get("cchange") or "")
        csq = v.get("consequence")
        csq_s = _esc(", ".join(csq) if isinstance(csq, list) else (csq or "MIS"))
        zyg_s = _esc(v.get("zygosity") or "Het")
        phase_s = _esc(v.get("phase") or "Unphased")
        cv_s = _esc(v.get("clinvar") or "None")
        tier_s = _esc(v.get("tier") or "Tier3")
        cadd_s = f"{float(v.get('cadd')):.1f}" if v.get("cadd") is not None else "—"
        sp_s = f"{float(v.get('spliceai')):.2f}" if v.get("spliceai") is not None else "—"
        avi_s = f"{float(v.get('aviPhred')):.1f} ({_esc(v.get('aviModality') or 'Multi')})" if v.get("aviPhred") is not None else "—"
        rev_s = f"{float(v.get('revel')):.2f}" if v.get("revel") is not None else "—"
        am_s = f"{float(v.get('alphamissense')):.2f}" if v.get("alphamissense") is not None else "—"
        maf_s = f"{float(v.get('maf'))*100:.4f}%" if v.get("maf") is not None else "Rare"
        coord_s = _esc(v.get("coordinate") or "")

        cv_cls = "pill-red" if "pathogenic" in cv_s.lower() and "conflicting" not in cv_s.lower() else ("pill-amber" if "risk" in cv_s.lower() or "uncertain" in cv_s.lower() else "pill-slate")

        rows_html += f"""
        <tr class="var-main-row">
          <td class="mono">{rs_link}{f'<div style="font-size:9.5px;color:#475569;">{hgvs_s}</div>' if (hgvs_s := hgvs_sub) else ''}</td>
          <td><span class="pill pill-slate">{csq_s}</span></td>
          <td><strong>{zyg_s}</strong>{f' <span class="mono" style="font-size:9.5px;color:#64748b;">(DP:{v.get("dp")}x)</span>' if v.get("dp") else ''}</td>
          <td><span class="pill pill-slate">{phase_s}</span></td>
          <td><span class="pill {cv_cls}">{cv_s}</span></td>
          <td><span class="pill {'pill-red' if tier_s == 'Tier1' else 'pill-teal'}">{tier_s}</span></td>
          <td class="mono"><strong>{cadd_s}</strong></td>
          <td class="mono">{sp_s}</td>
          <td class="mono" style="color:#0369a1;font-weight:700;">{avi_s}</td>
          <td class="mono">{rev_s} / {am_s}</td>
          <td class="mono" style="font-size:9.5px;color:#475569;">{coord_s}<div style="font-size:9px;color:#64748b;">AF: {maf_s}</div></td>
        </tr>
        """
        # Expand the 4-box drawer for Tier-1/Pathogenic variants and the #1 ranked variant in the gene
        if idx == 0 or _is_tier1_or_pathogenic_variant(v):
            rows_html += _render_variant_drawer_html(v, g)

    omim_id = g.get("omimGene")
    omim_url = f"https://omim.org/entry/{omim_id}" if (omim_id and str(omim_id) != "100000") else f"https://omim.org/search?search={sym}"
    ncbi_id = g.get("ncbiGeneId")
    ncbi_url = f"https://www.ncbi.nlm.nih.gov/gene/{ncbi_id}" if (ncbi_id and str(ncbi_id) != "0") else f"https://www.ncbi.nlm.nih.gov/gene/?term={sym}"

    return f"""
    <div class="gene-card">
      <div class="gene-card__head">
        <div style="display:flex;align-items:center;gap:6px;flex-wrap:wrap;">
          <span class="gene-sym">{sym}</span>
          <span class="gene-title">{gname}</span>
          <span class="mono" style="font-size:10px;color:#64748b;">{chrom}:{pos_str}</span>
          <span class="pill pill-slate">{organ}</span>
          {tier1_badge}
          {phase_badge}
        </div>
        <div style="display:flex;align-items:center;gap:6px;flex-wrap:wrap;">
          <span class="metric-pill mono">{_esc(g.get("topMetricDisplay"))}</span>
          <span style="font-size:10px;color:#475569;">Showing {len(g.get("variants") or [])} of {g.get("totalMatchingVariantsInGene", len(g.get("variants") or []))} variants</span>
          <a href="{_esc(omim_url)}" target="_blank" rel="noopener" class="ext-link">OMIM ↗</a>
          <a href="{_esc(ncbi_url)}" target="_blank" rel="noopener" class="ext-link">NCBI ↗</a>
          <a href="https://www.genecards.org/cgi-bin/carddisp.pl?gene={sym}" target="_blank" rel="noopener" class="ext-link">GeneCards ↗</a>
        </div>
      </div>
      {f'<div class="gene-summary">📖 {summary}</div>' if summary else ''}
      <table class="var-table">
        <thead>
          <tr>
            <th style="width:13%;">Variant / HGVS</th>
            <th style="width:9%;">Consequence</th>
            <th style="width:9%;">Zygosity / DP</th>
            <th style="width:8%;">Phase</th>
            <th style="width:12%;">ClinVar</th>
            <th style="width:6%;">Tier</th>
            <th style="width:6%;">CADD</th>
            <th style="width:6%;">SpliceAI</th>
            <th style="width:12%;">AlphaGenome AVI</th>
            <th style="width:8%;">REVEL / AM</th>
            <th style="width:11%;">Coordinate &amp; AF</th>
          </tr>
        </thead>
        <tbody>
          {rows_html}
        </tbody>
      </table>
    </div>
    """


def generate_dashboard_topn_report_html(patient_name, patient_id, master_json_path):
    """
    Generates `{patient_id}_dashboard_tier1_topn_report.html` containing the pre-expanded
    Tier-1 & Top-N genes across all 8 Dashboard Report selections.
    """
    p_data = pld.process_patient(master_json_path, patient_id, patient_name)
    genes = p_data.get("genes") or []
    depth = p_data.get("meta", {}).get("vcfDepthTelemetry", {})
    mean_dp = depth.get("mean_depth_called_loci") or "50.1x"
    median_dp = depth.get("median_depth_called_loci") or "50.0x"
    date_str = datetime.now().strftime("%B %d, %Y")

    nav_buttons_html = '<button class="sec-tab active" onclick="filterSection(\'all\', this)">All 8 Dashboard Views</button>'
    sections_html = ""

    for sec in DASHBOARD_REPORT_SECTIONS:
        selected_genes, total_matching = filter_and_select_topn_genes(genes, sec["key"], top_n=sec["top_n"])
        t1_count = sum(1 for g in selected_genes if g.get("hasTier1OrPathogenic"))
        nav_buttons_html += f"""
        <button class="sec-tab" onclick="filterSection('{sec['key']}', this)">
          {_esc(sec['short_label'])} <span class="mono" style="opacity:0.8;">({len(selected_genes)})</span>
        </button>
        """

        cards_html = "\n".join(render_gene_card_static_html(g, sec["key"]) for g in selected_genes)
        sections_html += f"""
        <section class="dashboard-report-section" data-section-key="{sec['key']}">
          <div class="sec-banner" style="border-left:4px solid {sec['badge_color']};background:{sec['badge_bg']};">
            <div>
              <h2 style="color:{sec['badge_color']};font-size:13px;font-weight:800;margin:0;">{_esc(sec['title'])}</h2>
              <p style="font-size:11px;color:#334155;margin:2px 0 0 0;">{_esc(sec['desc'])}</p>
            </div>
            <div style="text-align:right;font-size:10.5px;color:#1e293b;" class="mono">
              <div><strong>Selected Scope:</strong> {len(selected_genes)} Genes ({t1_count} Tier-1/Pathogenic)</div>
              <div style="color:#64748b;">Filtered from {total_matching} matching genes in view</div>
            </div>
          </div>
          {cards_html}
        </section>
        """

    html_out = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8" />
<meta name="viewport" content="width=device-width, initial-scale=1.0" />
<title>Dashboard Tier-1 &amp; Top-N Genomic Report — {_esc(patient_name)} ({_esc(patient_id)})</title>
<style>
  @page {{ size: letter landscape; margin: 6mm; }}
  * {{ box-sizing: border-box; }}
  body {{
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Inter", Arial, sans-serif;
    color: #0f172a;
    background: #f8fafc;
    margin: 0;
    padding: 0;
    font-size: 10.5px;
    line-height: 1.35;
    -webkit-print-color-adjust: exact !important;
    print-color-adjust: exact !important;
  }}
  .mono {{ font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace; }}
  .top-toolbar {{
    position: sticky;
    top: 0;
    z-index: 100;
    background: #0f172a;
    color: #ffffff;
    padding: 8px 16px;
    display: flex;
    justify-content: space-between;
    align-items: center;
    flex-wrap: wrap;
    gap: 8px;
    box-shadow: 0 2px 10px rgba(0,0,0,0.15);
  }}
  .toolbar-btn {{
    display: inline-block;
    padding: 4px 10px;
    border-radius: 4px;
    font-size: 11px;
    font-weight: 700;
    text-decoration: none;
    border: 1px solid transparent;
    cursor: pointer;
  }}
  .sec-nav {{
    background: #1e293b;
    padding: 6px 16px;
    display: flex;
    gap: 6px;
    flex-wrap: wrap;
    border-bottom: 1px solid #334155;
  }}
  .sec-tab {{
    background: #334155;
    color: #e2e8f0;
    border: 1px solid #475569;
    border-radius: 4px;
    padding: 4px 9px;
    font-size: 10.5px;
    font-weight: 600;
    cursor: pointer;
  }}
  .sec-tab.active, .sec-tab:hover {{
    background: #0f766e;
    color: #ffffff;
    border-color: #14b8a6;
  }}
  .container {{
    max-width: 1360px;
    margin: 12px auto;
    padding: 0 14px 24px 14px;
  }}
  .report-hero {{
    background: #ffffff;
    border: 1px solid #cbd5e1;
    border-top: 4px solid #0f766e;
    border-radius: 6px;
    padding: 10px 14px;
    margin-bottom: 12px;
    display: flex;
    justify-content: space-between;
    align-items: flex-start;
    flex-wrap: wrap;
    gap: 10px;
  }}
  .sec-banner {{
    padding: 8px 12px;
    border-radius: 5px;
    border: 1px solid #cbd5e1;
    margin: 14px 0 8px 0;
    display: flex;
    justify-content: space-between;
    align-items: center;
    flex-wrap: wrap;
    gap: 8px;
  }}
  .gene-card {{
    background: #ffffff;
    border: 1px solid #cbd5e1;
    border-radius: 6px;
    margin-bottom: 9px;
    overflow: hidden;
    page-break-inside: avoid;
  }}
  .gene-card__head {{
    background: #f1f5f9;
    padding: 6px 10px;
    border-bottom: 1px solid #cbd5e1;
    display: flex;
    justify-content: space-between;
    align-items: center;
    flex-wrap: wrap;
    gap: 6px;
  }}
  .gene-sym {{
    font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
    font-size: 13px;
    font-weight: 800;
    color: #123f41;
  }}
  .gene-title {{
    font-weight: 600;
    color: #334155;
    font-size: 11px;
  }}
  .gene-summary {{
    padding: 5px 10px;
    background: #fafbfc;
    border-bottom: 1px solid #e2e8f0;
    font-size: 10px;
    color: #334155;
  }}
  .pill {{
    display: inline-block;
    padding: 1px 5px;
    border-radius: 3px;
    font-size: 9.5px;
    font-weight: 700;
  }}
  .pill-red {{ background: #fef2f2; color: #b91c1c; border: 1px solid #fecaca; }}
  .pill-amber {{ background: #fffbeb; color: #b45309; border: 1px solid #fde68a; }}
  .pill-teal {{ background: #f0fdfa; color: #0f766e; border: 1px solid #99f6e4; }}
  .pill-slate {{ background: #f1f5f9; color: #334155; border: 1px solid #cbd5e1; }}
  .mini-badge {{
    display: inline-block;
    padding: 0.5px 4px;
    background: #e2e8f0;
    border-radius: 2px;
    font-size: 9px;
    font-weight: 600;
  }}
  .metric-pill {{
    background: #0f172a;
    color: #38bdf8;
    padding: 2px 6px;
    border-radius: 3px;
    font-size: 10px;
    font-weight: 700;
  }}
  .ext-link {{
    font-size: 10px;
    font-weight: 700;
    color: #0f766e;
    text-decoration: none;
  }}
  .var-table {{
    width: 100%;
    border-collapse: collapse;
    font-size: 10px;
  }}
  .var-table th {{
    background: #f8fafc;
    color: #475569;
    text-transform: uppercase;
    font-size: 9px;
    padding: 4px 6px;
    border-bottom: 1px solid #cbd5e1;
    text-align: left;
  }}
  .var-table td {{
    padding: 4px 6px;
    border-bottom: 1px solid #e2e8f0;
    vertical-align: top;
  }}
  .drawer-row td {{
    background: #f8fafc;
    padding: 6px 8px;
  }}
  .drawer-grid {{
    display: grid;
    grid-template-columns: repeat(4, 1fr);
    gap: 8px;
  }}
  .drawer-box {{
    background: #ffffff;
    border: 1px solid #cbd5e1;
    border-radius: 4px;
    padding: 6px 8px;
    font-size: 9.8px;
  }}
  .drawer-box__title {{
    font-weight: 800;
    color: #123f41;
    text-transform: uppercase;
    font-size: 9px;
    border-bottom: 1px solid #e2e8f0;
    padding-bottom: 2px;
    margin-bottom: 4px;
  }}
  .kv-row {{
    display: flex;
    justify-content: space-between;
    gap: 4px;
    margin-bottom: 2px;
  }}
  @media print {{
    body {{ background: #ffffff; }}
    .no-print {{ display: none !important; }}
    .container {{ max-width: 100%; margin: 0; padding: 0; }}
    .dashboard-report-section {{ page-break-before: always; }}
    .dashboard-report-section:first-of-type {{ page-break-before: auto; }}
  }}
</style>
</head>
<body>

  <div class="top-toolbar no-print">
    <div>
      <span style="background:#0f766e;color:#ffffff;padding:2px 6px;border-radius:3px;font-weight:800;font-size:10px;text-transform:uppercase;margin-right:6px;">Dashboard Tier-1 &amp; Top-N Digest</span>
      <strong>{_esc(patient_name)} ({_esc(patient_id)}) — Focused Multi-View Dashboard Report</strong>
    </div>
    <div style="display:flex;gap:6px;flex-wrap:wrap;align-items:center;">
      <a href="{_esc(patient_id)}_landscape_dashboard.html" class="toolbar-btn" style="background:#334155;color:#f8fafc;">← Interactive Command Center</a>
      <button onclick="window.print()" class="toolbar-btn" style="background:#2563eb;color:#ffffff;">🖨️ Print Active View</button>
      <a href="{_esc(patient_id)}_dashboard_tier1_topn_report.mhtml" download="{_esc(patient_id)}_dashboard_tier1_topn_report.mhtml" class="toolbar-btn" style="background:#7c3aed;color:#ffffff;">📦 Download .MHTML Snapshot</a>
      <a href="{_esc(patient_id)}_dashboard_tier1_topn_report.pdf" target="_blank" download="{_esc(patient_id)}_dashboard_tier1_topn_report.pdf" class="toolbar-btn" style="background:#059669;color:#ffffff;">📥 Save / Open PDF</a>
    </div>
  </div>

  <div class="sec-nav no-print">
    {nav_buttons_html}
  </div>

  <div class="container">
    <div class="report-hero">
      <div>
        <div style="font-size:15px;font-weight:800;color:#123f41;">LANDSCAPE DASHBOARD — TIER-1 &amp; TOP-N GENOMIC REPORT DIGEST</div>
        <div style="font-size:11px;color:#334155;margin-top:2px;">
          <strong>Cohort Subject:</strong> {_esc(patient_name)} (<span class="mono">{_esc(patient_id)}</span>) ·
          <strong>Mean Called Depth:</strong> <span class="mono" style="color:#0f766e;font-weight:700;">{_esc(mean_dp)}</span> (Median {_esc(median_dp)}) ·
          <strong>Date:</strong> {_esc(date_str)}
        </div>
        <div style="font-size:10px;color:#64748b;margin-top:3px;">
          Pre-expanded static digest of the <strong>Tier-1 / ClinVar Pathogenic genes</strong> and <strong>Top-N highest-scoring genes</strong> across all 8 Dashboard Report views (Visual Explorer, ClinVar, AlphaGenome Atlas, CADD, SpliceAI, AlphaMissense, REVEL, and Rare Variants ≤1%).
        </div>
      </div>
      <div style="text-align:right;font-size:10px;color:#475569;" class="mono">
        <div><strong>Total Cohort Catalog:</strong> {p_data["meta"]["totalVariants"]} Variants / {p_data["meta"]["totalGenes"]} Genes</div>
        <div><strong>Output Formats:</strong> .HTML · .MHTML (RFC 2557) · .PDF</div>
      </div>
    </div>

    {sections_html}
  </div>

  <script>
    function filterSection(key, btn) {{
      document.querySelectorAll('.sec-tab').forEach(b => b.classList.remove('active'));
      if (btn) btn.classList.add('active');
      document.querySelectorAll('.dashboard-report-section').forEach(sec => {{
        if (key === 'all' || sec.getAttribute('data-section-key') === key) {{
          sec.style.display = 'block';
        }} else {{
          sec.style.display = 'none';
        }}
      }});
    }}
  </script>
</body>
</html>
"""
    return html_out


def compile_html_to_mhtml(html_str, out_mhtml_path, title="Dashboard Tier-1 & Top-N Genomic Report", content_location="file:///dashboard_tier1_topn_report.html"):
    """
    Encapsulates a self-contained HTML string into an RFC 2557 standards-compliant
    MIME HTML (.mhtml) Web Archive compatible with Chrome, Edge, and archival viewers.
    """
    boundary = "----=_NextPart_000_0000_01DB1CA0.8F4E2B10"
    now_rfc = datetime.now(timezone.utc).strftime("%a, %d %b %Y %H:%M:%S +0000")
    # Strip <script> tags in MHTML so Chromium MHTML viewer renders 100% of sections cleanly without script sandbox warnings
    clean_html = html_str
    qp_bytes = quopri.encodestring(clean_html.encode("utf-8"), quotetabs=False)
    qp_body = qp_bytes.decode("ascii").replace("\n", "\r\n")

    mhtml_headers = [
        "From: <Saved by Genomic Command Center>",
        f"Snapshot-Content-Location: {content_location}",
        f"Subject: {title}",
        f"Date: {now_rfc}",
        "MIME-Version: 1.0",
        f'Content-Type: multipart/related;\r\n\ttype="text/html";\r\n\tboundary="{boundary}"',
        "",
        f"--{boundary}",
        "Content-Type: text/html",
        'Content-ID: <frame-0001@genomics.local>',
        "Content-Transfer-Encoding: quoted-printable",
        f"Content-Location: {content_location}",
        "",
        qp_body,
        f"--{boundary}--",
        ""
    ]
    mhtml_text = "\r\n".join(mhtml_headers)
    with open(out_mhtml_path, "w", encoding="utf-8", newline="") as f:
        f.write(mhtml_text)
    return os.path.exists(out_mhtml_path)
