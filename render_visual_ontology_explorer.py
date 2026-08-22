#!/usr/bin/env python3
"""
Visual Ontology Explorer & Master Hub - Iteration 2
Collapsible D3 Tree, Real Data Binding, Multi-level Domain Inference.
"""

import argparse
import json
import os
import re
import sys
import yaml
from pathlib import Path

def load_domain_registry():
    config_path = Path(__file__).parent / "config" / "ontology_domains.yaml"
    if config_path.exists():
        try:
            with open(config_path, "r", encoding="utf-8") as f:
                return yaml.safe_load(f).get("level1_systems", {})
        except Exception:
            pass
    return {}

def extract_omim_digits(omim_val) -> str:
    if not omim_val:
        return ""
    match = re.search(r'\d{6}|\d{5}|\d{4}', str(omim_val))
    return match.group(0) if match else ""

def get_max_spliceai(r):
    scores = []
    for k in ['spliceai_ds_ag', 'spliceai_ds_al', 'spliceai_ds_dg', 'spliceai_ds_dl']:
        val = r.get(k)
        if val is not None:
            try:
                scores.append(float(val))
            except:
                pass
    ev_val = (r.get("evidence") or {}).get("spliceai_max")
    if ev_val is not None:
        try:
            scores.append(float(ev_val))
        except:
            pass
    return max(scores) if scores else "N/A"

def infer_domain(r, domain_reg):
    # Combine context
    context_str = " ".join([
        str(r.get("clinvar_disease") or ""),
        " ".join(r.get("reason_codes", [])),
        " ".join((r.get("evidence") or {}).get("hpo_context", [])),
        (r.get("gene_hpo_term") or ""),
        ((r.get("gene_info") or {}).get("description") or "")
    ]).lower()
    
    best_l1 = None
    best_l2 = None
    
    for l1_key, l1_val in domain_reg.items():
        if l1_key.lower() in context_str or (l1_val.get("title") and l1_val["title"].lower().split()[0] in context_str):
            best_l1 = {"id": l1_val.get("id") or l1_key, "label": l1_val["title"], "color": l1_val.get("color", "#0ea5e9")}
            for l2_key, l2_val in l1_val.get("level2_subcategories", {}).items():
                l2_title = l2_val.get("title", "")
                # naive check
                if l2_key.lower() in context_str or any(w in context_str for w in l2_title.lower().split() if len(w)>4):
                    best_l2 = {"id": l2_val.get("id") or l2_key, "label": l2_title, "color": "#6366f1"}
                    break
            if not best_l2:
                best_l2 = {"id": "GEN_" + l1_key, "label": "General " + l1_val["title"], "color": "#6366f1"}
            break
            
    if not best_l1:
        best_l1 = {"id": "SYSTEM_OTHER", "label": "Other/Unclassified Systems", "color": "#64748b"}
        best_l2 = {"id": "SUBCAT_OTHER", "label": "General Phenotypes", "color": "#94a3b8"}
        
    return best_l1, best_l2

def generate_upgraded_visual_report(report_data: dict, output_filepath: str):
    domain_reg = load_domain_registry()
    
    patient_id = report_data.get("patient_id") or report_data.get("patient") or "DE_WGS_2026"
    run_date = report_data.get("run_date", "2026-08-21")
    
    raw_records = report_data.get("records") or report_data.get("monogenic_findings") or []
    polygenic = report_data.get("polygenic_findings") or []
    pharma = report_data.get("pharma_findings") or []

    genes_map = {}
    hpo_dict = {}
    
    for r in raw_records:
        hugo = r.get("hugo") or r.get("gene_symbol") or "Unknown"
        
        # Parse HPOs
        hpo_ids_raw = r.get("gene_hpo_id") or ""
        hpo_terms_raw = r.get("gene_hpo_term") or ""
        hpo_ids = [h.strip() for h in hpo_ids_raw.split(";") if h.strip()]
        hpo_terms = [h.strip() for h in hpo_terms_raw.split(";") if h.strip()]
        
        primary_hpo_id = "HP:0000118"
        primary_hpo_label = "Phenotypic abnormality"
        
        if hpo_ids:
            primary_hpo_id = hpo_ids[0]
            if len(hpo_terms) > 0:
                primary_hpo_label = hpo_terms[0]
                
        for i, h_id in enumerate(hpo_ids):
            if i < len(hpo_terms):
                hpo_dict[h_id] = hpo_terms[i]
            else:
                if h_id not in hpo_dict:
                    hpo_dict[h_id] = h_id

        if hugo not in genes_map:
            omim_digits = extract_omim_digits(r.get("omim_id") or r.get("omim_source"))
            
            # Robust NCBI desc
            ncbi_desc = (r.get("gene_info") or {}).get("description") or r.get("ncbi_description") or r.get("gene_desc")
            if not ncbi_desc:
                ncbi_desc = f"No detailed NCBI synopsis available for {hugo}."
                
            l1, l2 = infer_domain(r, domain_reg)
            
            genes_map[hugo] = {
                "gene_symbol": hugo,
                "ncbi_description": ncbi_desc,
                "associated_hpo_terms": hpo_ids,
                "associated_mondo_terms": [f"OMIM:{omim_digits}"] if omim_digits else [],
                "domain_l1": l1,
                "domain_l2": l2,
                "primary_hpo": {"id": primary_hpo_id, "label": primary_hpo_label},
                "variants": []
            }
            
        # Format variant metrics
        ev = r.get("evidence") or {}
        zyg = r.get("zygosity") or ev.get("zygosity") or "Heterozygous"
        phasing = r.get("phasing") or ev.get("phasing") or "Unphased"
        
        depth = ev.get("tot_reads") or r.get("tot_reads") or r.get("vcfinfo__tot_reads") or "N/A"
        quality = ev.get("qual") or r.get("phred") or r.get("qual") or "N/A"
        
        cadd = r.get("cadd_phred") or r.get("cadd") or "N/A"
        spliceai = get_max_spliceai(r)
        revel = r.get("revel") or r.get("revel_score") or "N/A"
        am_class = r.get("am_class") or "N/A"
        
        tier = r.get("cardio_tier") or r.get("tier") or "Tier 3"
        if "pathogenic" in str(r.get("clinvar_sig", "")).lower():
            tier = "Tier 1"
        elif "vus" in str(r.get("clinvar_sig", "")).lower() or "uncertain" in str(r.get("clinvar_sig", "")).lower():
            tier = "Tier 2"
            
        variant = {
            "rsid": r.get("rsid") or r.get("dbsnp") or "Novel Variant",
            "chromosome": r.get("chrom") or "Unknown",
            "position": r.get("pos") or "Unknown",
            "genotype": r.get("genotype") or f"{r.get('ref')}/{r.get('alt')}",
            "zygosity": zyg,
            "revel_score": revel,
            "impact_consequence": r.get("achange") or r.get("impact_consequence") or r.get("cchange") or "Unknown",
            "clinvar_significance": r.get("clinvar_sig") or "VUS",
            "clinvar_disease": r.get("clinvar_disease") or "",
            "phasing": phasing,
            "read_depth": depth,
            "read_quality": quality,
            "cadd_phred": cadd,
            "spliceai_max": spliceai,
            "am_class": am_class,
            "tier": tier,
            "gnomad4_af": r.get("gnomad4_af") or "N/A",
            "allofus_af": r.get("allofus_af") or "N/A",
            "reason_codes": r.get("reason_codes") or [],
            "transcript": r.get("transcript") or "Unknown",
            "pmid": r.get("denovo__PubmedID") or ""
        }
        genes_map[hugo]["variants"].append(variant)

    monogenic_findings = list(genes_map.values())
    
    # JSON serialize data for JS
    monogenic_json = json.dumps(monogenic_findings, indent=2)
    polygenic_json = json.dumps(polygenic, indent=2)
    pharma_json = json.dumps(pharma, indent=2)
    hpo_dict_json = json.dumps(hpo_dict, indent=2)

    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Visual Ontology Explorer - Genomics Report</title>
    <script src="https://cdn.tailwindcss.com"></script>
    <script src="https://d3js.org/d3.v7.min.js"></script>
    <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0/css/all.min.css">
    <style>
        :root {{ --bg-dark: #0f172a; --border-glow: #38bdf8; }}
        body {{ background-color: var(--bg-dark); color: #f8fafc; font-family: ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif; }}
        .tab-btn.active {{ border-bottom: 2px solid var(--border-glow); color: #38bdf8; }}
        .node-link {{ fill: none; stroke: #334155; stroke-opacity: 0.6; stroke-width: 1.5px; transition: stroke 0.3s; }}
        .node-circle {{ cursor: pointer; stroke: #1e293b; stroke-width: 1.5px; transition: fill 0.2s, stroke 0.2s, r 0.2s; }}
        .node-circle:hover {{ stroke: #38bdf8; stroke-width: 2.5px; }}
        ::-webkit-scrollbar {{ width: 6px; height: 6px; }}
        ::-webkit-scrollbar-track {{ background: #0f172a; }}
        ::-webkit-scrollbar-thumb {{ background: #334155; border-radius: 3px; }}
        ::-webkit-scrollbar-thumb:hover {{ background: #475569; }}
        
        .variant-details {{ transition: all 0.3s ease-in-out; overflow: hidden; }}
        .expand-btn {{ cursor: pointer; color: #94a3b8; transition: color 0.2s, transform 0.2s; }}
        .expand-btn:hover {{ color: #38bdf8; }}
        .expand-btn.open {{ color: #ef4444; }}
        .badge-reason {{ background: #1e293b; color: #94a3b8; border: 1px solid #334155; padding: 2px 6px; border-radius: 4px; font-size: 9px; font-weight: bold; }}
        .badge-reason.pathogenic {{ background: rgba(239, 68, 68, 0.15); color: #ef4444; border-color: rgba(239, 68, 68, 0.3); }}
    </style>
</head>
<body class="h-screen overflow-hidden flex flex-col">

    <header class="bg-slate-900 border-b border-slate-800 px-6 py-4 flex justify-between items-center shadow-lg shrink-0">
        <div class="flex items-center gap-3">
            <div class="bg-sky-500/20 text-sky-400 p-2 rounded-lg border border-sky-500/30"><i class="fa-solid fa-dna text-xl"></i></div>
            <div>
                <h1 class="text-lg font-extrabold tracking-tight text-white flex items-center gap-2">
                    Visual Ontology Explorer <span class="text-xs bg-sky-500/20 text-sky-400 px-2 py-0.5 rounded-full border border-sky-500/30">Genomics Report</span>
                </h1>
                <p class="text-xs text-slate-400 mt-0.5">Collapsible HPO Hierarchy & Multi-level Phenotype-Driven Triage</p>
            </div>
        </div>
        
        <div class="relative w-96 hidden md:block">
            <span class="absolute inset-y-0 left-0 flex items-center pl-3 text-slate-500"><i class="fa-solid fa-magnifying-glass"></i></span>
            <input type="text" id="globalSearch" class="w-full bg-slate-950 border border-slate-800 rounded-lg pl-10 pr-4 py-2 text-sm text-slate-200 focus:outline-none focus:border-sky-500" placeholder="Search HPO, traits, genes (Press Enter)" onkeypress="if(event.key === 'Enter') onGlobalSearch(this.value)">
        </div>

        <div class="flex items-center gap-4 text-xs text-right">
            <div><div class="text-slate-300">Patient ID: <span class="font-bold text-white">{patient_id}</span></div><div class="text-slate-500">Date: {run_date} | GRCh38</div></div>
            <div class="bg-slate-800 text-slate-200 w-8 h-8 rounded-full flex items-center justify-center font-bold border border-slate-700">PT</div>
        </div>
    </header>

    <main class="flex-1 flex overflow-hidden">
        
        <!-- Left Pane: D3 Collapsible Tree -->
        <section class="w-2/5 flex flex-col border-r border-slate-800 bg-slate-950/40 relative min-w-[450px]">
            <div class="p-4 bg-slate-900/60 border-b border-slate-800 flex justify-between items-center z-10">
                <div>
                    <h2 class="text-sm font-bold text-slate-200 flex items-center gap-2">
                        <i class="fa-solid fa-diagram-project text-sky-400"></i> HPO Ontology Explorer
                    </h2>
                    <p class="text-[10px] text-slate-500 uppercase font-bold tracking-wider mt-1">System → Subcategory → Phenotype → Gene</p>
                </div>
                <div class="flex gap-2">
                    <button onclick="resetZoom()" class="bg-slate-800 hover:bg-slate-700 text-xs px-2 py-1 rounded text-slate-300 border border-slate-700"><i class="fa-solid fa-arrows-to-center"></i> Reset</button>
                    <button onclick="expandAll()" class="bg-slate-800 hover:bg-slate-700 text-xs px-2 py-1 rounded text-slate-300 border border-slate-700"><i class="fa-solid fa-expand"></i> Expand All</button>
                </div>
            </div>

            <div class="flex-1 relative overflow-auto" id="canvasContainer">
                <svg id="hpoTreeSvg" class="w-full h-full absolute inset-0"></svg>
            </div>
        </section>

        <!-- Right Pane: Gene Inspector -->
        <section class="w-3/5 flex flex-col bg-slate-900/20 overflow-hidden relative">
            <div class="p-6 bg-slate-900/40 border-b border-slate-800 shrink-0">
                <div class="flex items-start gap-4">
                    <div class="w-12 h-12 rounded-xl bg-sky-500/10 text-sky-400 flex items-center justify-center border border-sky-500/20 text-2xl font-black shadow-inner" id="headerIcon">🧬</div>
                    <div>
                        <div class="flex items-center gap-2">
                            <h3 class="text-2xl font-extrabold text-white tracking-tight" id="headerTitle">Select a Gene node</h3>
                            <span class="px-2.5 py-0.5 rounded text-[10px] font-bold uppercase tracking-wider bg-slate-800 text-slate-400 border border-slate-700 hidden" id="headerTag"></span>
                        </div>
                        <p class="text-sm text-slate-400 font-semibold mt-1" id="headerSubtitle">Click on any gene node in the tree to view Gene Inspector details.</p>
                    </div>
                </div>
            </div>

            <div class="bg-slate-900/60 border-b border-slate-800 px-6 flex gap-6 shrink-0 z-10" id="tabsNavBar">
                <button onclick="switchTab('overview')" class="tab-btn py-3 text-sm font-bold text-slate-400 hover:text-white transition-all active" id="tab-overview"><i class="fa-solid fa-address-card"></i> Overview</button>
                <button onclick="switchTab('variants')" class="tab-btn py-3 text-sm font-bold text-slate-400 hover:text-white transition-all" id="tab-variants"><i class="fa-solid fa-vial-virus"></i> Variants <span class="bg-slate-800 text-slate-400 text-[10px] px-1.5 py-0.5 rounded-full ml-1" id="variantCountTag">0</span></button>
                <button onclick="switchTab('phenotypes')" class="tab-btn py-3 text-sm font-bold text-slate-400 hover:text-white transition-all" id="tab-phenotypes"><i class="fa-solid fa-stethoscope"></i> Phenotypes</button>
                <button onclick="switchTab('polygenic')" class="tab-btn py-3 text-sm font-bold text-slate-400 hover:text-white transition-all" id="tab-polygenic"><i class="fa-solid fa-chart-line"></i> Polygenic Risk</button>
                <button onclick="switchTab('publications')" class="tab-btn py-3 text-sm font-bold text-slate-400 hover:text-white transition-all" id="tab-publications"><i class="fa-solid fa-book-open"></i> Publications</button>
            </div>

            <div class="flex-1 overflow-y-auto p-6" id="tabContent">
                <div id="placeholderContent" class="h-full flex flex-col items-center justify-center text-center text-slate-500 py-12">
                    <div class="w-16 h-16 rounded-full bg-slate-800 flex items-center justify-center text-3xl mb-4 animate-bounce">💡</div>
                    <h4 class="text-base font-bold text-slate-300">No Gene Selected</h4>
                </div>

                <div id="overviewContent" class="hidden space-y-6">
                    <div class="bg-slate-800/40 border border-slate-800/60 p-5 rounded-xl">
                        <h4 class="text-xs font-extrabold uppercase tracking-wider text-sky-400 mb-2">NCBI Gene Summary</h4>
                        <p class="text-sm leading-relaxed text-slate-300 font-medium" id="ncbiSummary"></p>
                        <div class="grid grid-cols-2 gap-4 mt-4 border-t border-slate-800/60 pt-4 text-xs">
                            <div><span class="text-slate-500">OMIM Number:</span> <a id="mimLink" target="_blank" class="text-sky-400 hover:underline font-bold"></a></div>
                            <div><span class="text-slate-500">Transcript:</span> <span id="transcripts" class="text-white font-mono"></span></div>
                            <div class="col-span-2"><span class="text-slate-500">Associated Disease:</span> <span id="clinvarDisease" class="text-white"></span></div>
                        </div>
                    </div>
                </div>

                <div id="variantsContent" class="hidden space-y-4">
                    <div class="w-full overflow-x-auto bg-slate-800/30 border border-slate-700/60 rounded-lg shadow-lg">
                        <table class="w-full text-left text-sm text-slate-300">
                            <thead class="text-xs text-slate-400 uppercase bg-slate-900/80 border-b border-slate-700/80">
                                <tr>
                                    <th class="px-4 py-3 w-8"></th>
                                    <th class="px-4 py-3">Variant (RSID / Pos)</th>
                                    <th class="px-4 py-3">Alleles / Effect</th>
                                    <th class="px-4 py-3">ClinVar Signif.</th>
                                    <th class="px-4 py-3">Phasing / Tier</th>
                                </tr>
                            </thead>
                            <tbody id="variantTableBody" class="divide-y divide-slate-700/50"></tbody>
                        </table>
                    </div>
                </div>

                <div id="phenotypesContent" class="hidden space-y-4">
                    <div class="bg-slate-800/40 border border-slate-800/60 p-5 rounded-xl">
                        <h4 class="text-xs font-extrabold uppercase tracking-wider text-sky-400 mb-3">Associated Clinical Phenotypes</h4>
                        <div class="flex flex-col gap-2" id="hpoTags"></div>
                    </div>
                </div>
                
                <div id="polygenicContent" class="hidden space-y-4">
                    <div class="bg-slate-800/40 border border-slate-800/60 p-5 rounded-xl">
                        <h4 class="text-xs font-extrabold uppercase tracking-wider text-sky-400 mb-3 flex items-center justify-between">
                            <span>📊 Polygenic Risk & Interpretation (EFO)</span>
                            <span class="text-[10px] bg-slate-800 text-slate-400 px-2 py-0.5 rounded font-mono">PGS Catalog</span>
                        </h4>
                        <div class="grid grid-cols-1 md:grid-cols-2 gap-4" id="prsGrid"></div>
                    </div>
                </div>

                <div id="publicationsContent" class="hidden space-y-4">
                    <div class="bg-slate-900 border border-slate-800 rounded-xl p-4 shadow-md">
                        <div class="flex justify-between items-start mb-2">
                            <span class="bg-sky-500/10 text-sky-400 border border-sky-500/20 text-[9px] px-2 py-0.5 rounded font-extrabold uppercase tracking-wider">Literature / GWAS</span>
                        </div>
                        <h5 class="text-sm font-extrabold text-white leading-snug" id="pubTitle">Clinical significance of variants</h5>
                        <div class="flex justify-between items-center text-[10px] text-slate-500 pt-3 mt-3 border-t border-slate-800/60">
                            <a id="pubLink" href="#" target="_blank" class="text-sky-400 hover:underline flex items-center gap-1"><i class="fa-solid fa-arrow-up-right-from-square"></i> Open in PubMed</a>
                        </div>
                    </div>
                </div>

            </div>
        </section>
    </main>

    <script>
        const monogenicFindings = {monogenic_json};
        const polygenicFindings = {polygenic_json};
        const hpoDict = {hpo_dict_json};

        let currentActiveGene = null;
        let rootNode = null;
        let svgGroup = null;
        let d3ZoomBehavior = null;
        let treeLayout = null;
        let d3Link = null;
        let d3Node = null;
        const iWidth = 800;
        const margin = {{top: 20, right: 120, bottom: 20, left: 120}};

        function resolveHpo(id) {{ return hpoDict[id] || id; }}

        function buildHierarchy() {{
            const tree = {{ name: "Patient", type: "root", children: [] }};
            const l1Map = {{}};

            monogenicFindings.forEach(gene => {{
                const l1 = gene.domain_l1;
                const l2 = gene.domain_l2;
                const hpo = gene.primary_hpo;
                
                if (!l1Map[l1.id]) {{
                    l1Map[l1.id] = {{ name: l1.label, type: "level1", color: l1.color, children: [], l2Map: {{}} }};
                    tree.children.push(l1Map[l1.id]);
                }}
                
                if (!l1Map[l1.id].l2Map[l2.id]) {{
                    l1Map[l1.id].l2Map[l2.id] = {{ name: l2.label, type: "level2", color: l2.color, children: [], hpoMap: {{}} }};
                    l1Map[l1.id].children.push(l1Map[l1.id].l2Map[l2.id]);
                }}
                
                if (!l1Map[l1.id].l2Map[l2.id].hpoMap[hpo.id]) {{
                    l1Map[l1.id].l2Map[l2.id].hpoMap[hpo.id] = {{ name: hpo.label, type: "phenotype", color: "#10b981", children: [] }};
                    l1Map[l1.id].l2Map[l2.id].children.push(l1Map[l1.id].l2Map[l2.id].hpoMap[hpo.id]);
                }}
                
                let tier = "Tier 3";
                if(gene.variants.length > 0) tier = gene.variants[0].tier;
                const gColor = tier === "Tier 1" ? "#ef4444" : (tier === "Tier 2" ? "#f59e0b" : "#a855f7");
                
                l1Map[l1.id].l2Map[l2.id].hpoMap[hpo.id].children.push({{
                    name: gene.gene_symbol, type: "gene", color: gColor
                }});
            }});
            return tree;
        }}

        function initializeD3Tree() {{
            const data = buildHierarchy();
            const container = document.getElementById('canvasContainer');
            
            const svg = d3.select("#hpoTreeSvg");
            svg.selectAll("*").remove();

            svgGroup = svg.append("g").attr("transform", `translate(${{margin.left}},${{margin.top}})`);
            
            d3ZoomBehavior = d3.zoom().scaleExtent([0.1, 4]).on("zoom", (e) => {{
                svgGroup.attr("transform", e.transform);
            }});
            svg.call(d3ZoomBehavior).on("dblclick.zoom", null);

            rootNode = d3.hierarchy(data);
            rootNode.x0 = container.clientHeight / 2;
            rootNode.y0 = 0;

            // Collapse after level 1
            rootNode.descendants().forEach(d => {{
                if (d.depth > 1) {{
                    d._children = d.children;
                    d.children = null;
                }}
            }});

            treeLayout = d3.tree().nodeSize([25, 220]);
            updateTree(rootNode);
            
            // Auto-select first Tier 1 gene
            const firstGene = monogenicFindings.find(g => g.variants.some(v => v.tier === 'Tier 1'));
            if(firstGene) inspectGene(firstGene.gene_symbol);
        }}

        function updateTree(source) {{
            const treeData = treeLayout(rootNode);
            const nodes = treeData.descendants();
            const links = treeData.descendants().slice(1);

            nodes.forEach(d => {{ d.y = d.depth * 220; }});

            const node = svgGroup.selectAll('g.node')
                .data(nodes, d => d.id || (d.id = ++window.i));

            const nodeEnter = node.enter().append('g')
                .attr('class', 'node')
                .attr("transform", d => `translate(${{source.y0}},${{source.x0}})`)
                .on('click', clickNode);

            nodeEnter.append('circle')
                .attr('class', 'node-circle')
                .attr('r', 1e-6)
                .style("fill", d => d._children ? "#1e293b" : d.data.color)
                .style("stroke", d => d.data.color);

            nodeEnter.append('text')
                .attr("dy", ".35em")
                .attr("x", d => d.children || d._children ? -13 : 13)
                .attr("text-anchor", d => d.children || d._children ? "end" : "start")
                .text(d => d.data.name.length > 25 ? d.data.name.substring(0,23)+'...' : d.data.name)
                .attr("fill", "#e2e8f0").attr("font-size", "11px");

            const nodeUpdate = nodeEnter.merge(node);
            nodeUpdate.transition().duration(400)
                .attr("transform", d => `translate(${{d.y}},${{d.x}})`);

            nodeUpdate.select('circle.node-circle')
                .attr('r', d => d.data.type === 'gene' ? 7 : (d.data.type === 'root' ? 10 : 8))
                .style("fill", d => d._children ? "#1e293b" : d.data.color)
                .attr('cursor', 'pointer');

            const nodeExit = node.exit().transition().duration(400)
                .attr("transform", d => `translate(${{source.y}},${{source.x}})`)
                .remove();
            nodeExit.select('circle').attr('r', 1e-6);

            const link = svgGroup.selectAll('path.node-link')
                .data(links, d => d.id);

            const linkEnter = link.enter().insert('path', "g")
                .attr("class", "node-link")
                .attr('d', d => {{
                    const o = {{x: source.x0, y: source.y0}};
                    return diagonal(o, o);
                }});

            const linkUpdate = linkEnter.merge(link);
            linkUpdate.transition().duration(400)
                .attr('d', d => diagonal(d, d.parent));

            link.exit().transition().duration(400)
                .attr('d', d => {{
                    const o = {{x: source.x, y: source.y}};
                    return diagonal(o, o);
                }}).remove();

            nodes.forEach(d => {{ d.x0 = d.x; d.y0 = d.y; }});
        }}

        function diagonal(s, d) {{
            return `M ${{s.y}} ${{s.x}}
                    C ${{(s.y + d.y) / 2}} ${{s.x}},
                      ${{(s.y + d.y) / 2}} ${{d.x}},
                      ${{d.y}} ${{d.x}}`;
        }}

        function clickNode(event, d) {{
            if (d.data.type === 'gene') {{
                inspectGene(d.data.name);
                return;
            }}
            if (d.children) {{
                d._children = d.children;
                d.children = null;
            }} else {{
                d.children = d._children;
                d._children = null;
            }}
            updateTree(d);
        }}
        
        function expandAll() {{
            rootNode.descendants().forEach(d => {{
                if(d._children) {{ d.children = d._children; d._children = null; }}
            }});
            updateTree(rootNode);
        }}

        function resetZoom() {{
            d3.select("#hpoTreeSvg").transition().duration(500).call(d3ZoomBehavior.transform, d3.zoomIdentity.translate(margin.left, margin.top));
        }}

        function inspectGene(geneSymbol) {{
            currentActiveGene = geneSymbol;
            const record = monogenicFindings.find(f => f.gene_symbol === geneSymbol);
            if (!record) return;

            document.getElementById("placeholderContent").classList.add("hidden");
            document.getElementById("overviewContent").classList.remove("hidden");
            switchTab("overview");

            document.getElementById("headerTitle").innerText = geneSymbol;
            document.getElementById("headerTag").classList.remove("hidden");
            document.getElementById("headerTag").innerText = record.variants[0]?.tier || "TIER";

            document.getElementById("ncbiSummary").innerText = record.ncbi_description;
            const omim = record.associated_mondo_terms.length > 0 ? record.associated_mondo_terms[0] : "";
            document.getElementById("mimLink").innerText = omim || "Search OMIM";
            document.getElementById("mimLink").href = omim ? "https://omim.org/entry/" + omim.replace("OMIM:", "") : "https://omim.org/search/?search=" + geneSymbol;
            document.getElementById("transcripts").innerText = record.variants[0]?.transcript || "Canonical";
            document.getElementById("clinvarDisease").innerText = record.variants[0]?.clinvar_disease || "Not specified";

            const prsGrid = document.getElementById("prsGrid");
            prsGrid.innerHTML = "";
            polygenicFindings.forEach(prs => {{
                let riskColor = prs.risk_category === "HIGH" ? "text-red-400" : (prs.risk_category === "MODERATE" ? "text-yellow-400" : "text-sky-400");
                let prsBg = prs.risk_category === "HIGH" ? "bg-red-500/20" : (prs.risk_category === "MODERATE" ? "bg-yellow-500/20" : "bg-sky-500/20");
                prsGrid.innerHTML += `
                    <div class="bg-slate-900 border border-slate-800 p-4 rounded-lg text-center">
                        <span class="text-[10px] uppercase font-bold text-slate-400">` + prs.trait_name + `</span>
                        <div class="text-2xl font-black text-white mt-1">` + prs.percentile + ` %</div>
                        <span class="text-[10px] font-bold mt-2 inline-block ` + riskColor + ` ` + prsBg + ` px-2 rounded-full border border-current/20">` + prs.risk_category + `</span>
                    </div>`;
            }});

            document.getElementById("variantCountTag").innerText = record.variants.length;
            const vBody = document.getElementById("variantTableBody");
            vBody.innerHTML = "";

            record.variants.forEach((v, idx) => {{
                const sigCls = v.clinvar_significance.toLowerCase().includes("pathogenic") ? "text-red-400 bg-red-500/10 border-red-500/20" : "text-yellow-400 bg-yellow-500/10 border-yellow-500/20";
                
                let caddCls = "text-slate-200";
                if(v.cadd_phred !== "N/A" && parseFloat(v.cadd_phred) >= 20) caddCls = "text-amber-400";
                if(v.cadd_phred !== "N/A" && parseFloat(v.cadd_phred) >= 30) caddCls = "text-red-400";
                
                const amCls = v.am_class.includes("pathogenic") ? "text-red-400" : "text-emerald-400";
                
                let reasonsHtml = (v.reason_codes || []).map(r => {{
                    let c = "badge-reason";
                    if(r.includes("PATH") || r.includes("PVS") || r.includes("CONFLICT")) c += " pathogenic";
                    return `<span class="${{c}}">${{r}}</span>`;
                }}).join(" ");

                vBody.innerHTML += `
                <tr class="hover:bg-slate-800/40 transition-colors">
                    <td class="px-4 py-3">
                        <button class="expand-btn" id="btn-exp-${{idx}}" onclick="toggleDrawer(${{idx}})">
                            <i class="fa-solid fa-plus text-lg"></i>
                        </button>
                    </td>
                    <td class="px-4 py-3">
                        <div class="font-mono font-bold text-white">${{v.rsid}}</div>
                        <div class="text-[10px] text-slate-500">${{v.chromosome}}:${{v.position}}</div>
                    </td>
                    <td class="px-4 py-3">
                        <div class="font-mono text-slate-300 bg-slate-800 px-1.5 py-0.5 rounded inline-block text-xs">${{v.genotype}}</div>
                        <div class="text-xs text-slate-400 mt-1">${{v.impact_consequence}}</div>
                    </td>
                    <td class="px-4 py-3">
                        <span class="text-[10px] font-bold px-2 py-0.5 rounded-full uppercase border ${{sigCls}}">${{v.clinvar_significance}}</span>
                    </td>
                    <td class="px-4 py-3">
                        <div class="text-xs font-semibold text-indigo-400">${{v.phasing}}</div>
                        <div class="text-[10px] text-slate-500">${{v.tier}}</div>
                    </td>
                </tr>
                <tr id="drawer-${{idx}}" class="variant-details hidden bg-slate-950/60 border-b border-slate-700/80">
                    <td colspan="5" class="px-6 py-4">
                        <div class="grid grid-cols-2 md:grid-cols-4 gap-4 mb-3">
                            <div class="bg-slate-900 border border-slate-800 p-3 rounded-md">
                                <div class="text-[10px] text-slate-500 font-bold uppercase mb-1">CADD Phred</div>
                                <div class="font-mono ${{caddCls}} font-bold text-sm">${{v.cadd_phred}}</div>
                            </div>
                            <div class="bg-slate-900 border border-slate-800 p-3 rounded-md">
                                <div class="text-[10px] text-slate-500 font-bold uppercase mb-1">Read Depth & Qual</div>
                                <div class="text-slate-200 font-bold text-xs">${{v.read_depth}} <span class="text-slate-500 font-normal">| Q:${{v.read_quality}}</span></div>
                            </div>
                            <div class="bg-slate-900 border border-slate-800 p-3 rounded-md">
                                <div class="text-[10px] text-slate-500 font-bold uppercase mb-1">AlphaMissense</div>
                                <div class="${{amCls}} font-bold text-xs">${{v.am_class}}</div>
                            </div>
                            <div class="bg-slate-900 border border-slate-800 p-3 rounded-md">
                                <div class="text-[10px] text-slate-500 font-bold uppercase mb-1">SpliceAI Delta Max</div>
                                <div class="font-mono text-slate-200 font-bold text-sm">${{v.spliceai_max}}</div>
                            </div>
                        </div>
                        <div class="grid grid-cols-2 md:grid-cols-4 gap-4">
                            <div class="bg-slate-900 border border-slate-800 p-3 rounded-md">
                                <div class="text-[10px] text-slate-500 font-bold uppercase mb-1">gnomAD4 AF</div>
                                <div class="font-mono text-slate-200 text-xs">${{v.gnomad4_af}}</div>
                            </div>
                            <div class="bg-slate-900 border border-slate-800 p-3 rounded-md">
                                <div class="text-[10px] text-slate-500 font-bold uppercase mb-1">All of Us AF</div>
                                <div class="font-mono text-slate-200 text-xs">${{v.allofus_af}}</div>
                            </div>
                            <div class="bg-slate-900 border border-slate-800 p-3 rounded-md col-span-2">
                                <div class="text-[10px] text-slate-500 font-bold uppercase mb-1">ACMG Reason Codes</div>
                                <div class="flex flex-wrap gap-1">${{reasonsHtml}}</div>
                            </div>
                        </div>
                    </td>
                </tr>
                `;
            }});

            const hpoTags = document.getElementById("hpoTags");
            hpoTags.innerHTML = "";
            record.associated_hpo_terms.forEach(hpo => {{
                hpoTags.innerHTML += `
                    <div class="bg-slate-900 border border-slate-800 px-4 py-2 rounded-lg text-sm flex items-center justify-between">
                        <span class="font-semibold text-slate-200">${{resolveHpo(hpo)}}</span>
                        <span class="text-[10px] font-mono text-sky-400">${{hpo}}</span>
                    </div>
                `;
            }});
            
            document.getElementById("pubTitle").innerText = `Clinical significance of variants in ${{geneSymbol}}`;
            const pmid = record.variants[0]?.pmid;
            const pmidUrl = pmid ? `https://pubmed.ncbi.nlm.nih.gov/${{pmid}}` : `https://pubmed.ncbi.nlm.nih.gov/?term=${{geneSymbol}}`;
            document.getElementById("pubLink").href = pmidUrl;
        }}

        function toggleDrawer(idx) {{
            const drawer = document.getElementById(`drawer-${{idx}}`);
            const btn = document.getElementById(`btn-exp-${{idx}}`);
            if (drawer.classList.contains('hidden')) {{
                drawer.classList.remove('hidden');
                btn.classList.add('open');
                btn.innerHTML = '<i class="fa-solid fa-minus text-lg"></i>';
            }} else {{
                drawer.classList.add('hidden');
                btn.classList.remove('open');
                btn.innerHTML = '<i class="fa-solid fa-plus text-lg"></i>';
            }}
        }}

        function switchTab(tabId) {{
            document.querySelectorAll('.tab-btn').forEach(btn => btn.classList.remove('active'));
            document.getElementById("tab-" + tabId).classList.add('active');
            ['overview', 'variants', 'phenotypes', 'polygenic', 'publications'].forEach(t => {{
                document.getElementById(t + "Content").classList.add("hidden");
            }});
            document.getElementById(tabId + "Content").classList.remove("hidden");
        }}

        function onGlobalSearch(term) {{
            term = term.toLowerCase().trim();
            if (!term) return;
            
            // Search tree
            let foundNode = null;
            rootNode.descendants().forEach(d => {{
                if (d.data.name.toLowerCase().includes(term) || (d.data.type==='phenotype' && d.data.id && d.data.id.toLowerCase().includes(term))) {{
                    if(!foundNode || d.data.type === 'gene') foundNode = d; 
                }}
            }});
            
            if (foundNode) {{
                let curr = foundNode;
                while (curr.parent) {{
                    curr.parent.children = curr.parent.children || curr.parent._children;
                    curr.parent._children = null;
                    curr = curr.parent;
                }}
                updateTree(rootNode);
                
                if (foundNode.data.type === 'gene') {{
                    inspectGene(foundNode.data.name);
                }}
            }}
        }}

        window.addEventListener('load', () => {{
            window.i = 0;
            initializeD3Tree();
        }});
    </script>
</body>
</html>"""

    with open(output_filepath, "w", encoding="utf-8") as f:
        f.write(html_content)
    print(f"Successfully generated visual upgraded HTML at: {output_filepath}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Genomic & Polygenic Visual Explorer HTML Generator")
    parser.add_argument("-i", "--input", help="Path to input JSON file containing variant report data")
    parser.add_argument("-o", "--output", default="reports/visual_ontology_explorer.html", help="Path to output HTML")
    parser.add_argument("-d", "--demo", "--mock", action="store_true", help="Generate report using mel_actionable dataset")
    args = parser.parse_args()

    out_dir = os.path.dirname(args.output)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    if args.demo:
        mel_path = Path(__file__).parent / "logs" / "mel_actionable.json"
        if not mel_path.exists():
            print(f"Error: Demo file not found at {mel_path}")
            sys.exit(1)
        with open(mel_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        generate_upgraded_visual_report(data, args.output)
    else:
        if not args.input or not os.path.exists(args.input):
            print("Error: Input file missing.")
            sys.exit(1)
        with open(args.input, "r", encoding="utf-8") as f:
            data = json.load(f)
        generate_upgraded_visual_report(data, args.output)
