#!/usr/bin/env python3
import argparse, json, os, re, sys, yaml
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

def extract_omim_digits(omim_val):
    if not omim_val: return ""
    match = re.search(r'\d{6}|\d{5}|\d{4}', str(omim_val))
    return match.group(0) if match else ""

def get_max_spliceai(r):
    scores = []
    for k in ['spliceai_ds_ag', 'spliceai_ds_al', 'spliceai_ds_dg', 'spliceai_ds_dl']:
        val = r.get(k)
        if val is not None:
            try: scores.append(float(val))
            except: pass
    ev_val = (r.get("evidence") or {}).get("spliceai_max")
    if ev_val is not None:
        try: scores.append(float(ev_val))
        except: pass
    return max(scores) if scores else "N/A"

def generate_upgraded_visual_report(report_data: dict, output_filepath: str):
    domain_reg = load_domain_registry()
    patient_id = report_data.get("patient_id") or report_data.get("patient") or "DE_WGS_2026"
    run_date = report_data.get("run_date", "2026-08-21")
    
    raw_records = report_data.get("records") or report_data.get("monogenic_findings") or []

    genes_map = {}
    
    for r in raw_records:
        hugo = r.get("hugo") or r.get("gene_symbol") or "Unknown"
        hpo_ids_raw = r.get("gene_hpo_id") or ""
        hpo_terms_raw = r.get("gene_hpo_term") or ""
        hpo_ids = [h.strip() for h in hpo_ids_raw.split(";") if h.strip()]
        hpo_terms = [h.strip() for h in hpo_terms_raw.split(";") if h.strip()]

        if hugo not in genes_map:
            omim = extract_omim_digits(r.get("omim_id") or r.get("omim_source"))
            ncbi_desc = (r.get("gene_info") or {}).get("description") or r.get("ncbi_description") or r.get("gene_desc") or ""
            
            mapped_l1, mapped_l2 = None, None
            for l1_key, l1_val in domain_reg.items():
                for l2_key, l2_val in l1_val.get("level2_subcategories", {}).items():
                    domain_hpos = set(l2_val.get("hpo_terms", []))
                    if any(hid in domain_hpos for hid in hpo_ids):
                        mapped_l1 = {"id": l1_val.get("id", l1_key), "label": l1_val.get("title", l1_key), "color": l1_val.get("color", "#0ea5e9")}
                        mapped_l2 = {"id": l2_val.get("id", l2_key), "label": l2_val.get("title", l2_key), "color": "#6366f1"}
                        break
                if mapped_l1: break
                
            if not mapped_l1:
                mapped_l1 = {"id": "SYSTEM_OTHER", "label": "Other/Unclassified", "color": "#64748b"}
                mapped_l2 = {"id": "SUBCAT_OTHER", "label": "General Findings", "color": "#94a3b8"}
                
            genes_map[hugo] = {
                "gene_symbol": hugo,
                "ncbi_description": ncbi_desc,
                "omim": omim,
                "domain_l1": mapped_l1,
                "domain_l2": mapped_l2,
                "variants": []
            }
            
        ev = r.get("evidence") or {}
        zyg = r.get("zygosity") or ev.get("zygosity") or "Heterozygous"
        depth = ev.get("tot_reads") or r.get("tot_reads") or r.get("vcfinfo__tot_reads") or "N/A"
        quality = ev.get("qual") or r.get("phred") or r.get("qual") or "N/A"
        cadd = r.get("cadd_phred") or r.get("cadd") or "N/A"
        spliceai = get_max_spliceai(r)
        revel = r.get("revel") or r.get("revel_score") or "N/A"
        
        tier = r.get("cardio_tier") or r.get("tier") or "Tier 3"
        sig = str(r.get("clinvar_sig", "")).lower()
        if "pathogenic" in sig and "conflicting" not in sig: tier = "Tier 1"
        elif "vus" in sig or "uncertain" in sig or "conflicting" in sig: tier = "Tier 2"
            
        genes_map[hugo]["variants"].append({
            "rsid": r.get("rsid") or r.get("dbsnp") or "Novel",
            "genotype": r.get("genotype") or f"{r.get('ref')}>{r.get('alt')}",
            "zygosity": zyg,
            "impact_consequence": r.get("achange") or r.get("impact_consequence") or r.get("cchange") or "Unknown",
            "clinvar_significance": r.get("clinvar_sig") or "VUS",
            "read_depth": depth,
            "read_quality": quality,
            "cadd_phred": cadd,
            "spliceai_max": spliceai,
            "revel_score": revel,
            "tier": tier
        })

    monogenic_findings = list(genes_map.values())
    
    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Ontology Master Hub</title>
    <script src="https://cdn.tailwindcss.com"></script>
    <script src="https://d3js.org/d3.v7.min.js"></script>
    <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0/css/all.min.css">
    <style>
        :root {{ --bg-dark: #0f172a; --border-glow: #38bdf8; }}
        body {{ background-color: var(--bg-dark); color: #f8fafc; font-family: ui-sans-serif, system-ui, -apple-system, sans-serif; }}
        .tab-btn.active {{ border-bottom: 2px solid var(--border-glow); color: #38bdf8; }}
        .node-link {{ fill: none; stroke: #334155; stroke-opacity: 0.6; stroke-width: 1.5px; transition: stroke 0.3s; }}
        .node-circle {{ cursor: pointer; stroke: #1e293b; stroke-width: 1.5px; transition: fill 0.2s, stroke 0.2s, r 0.2s; }}
        .node-circle:hover {{ stroke: #38bdf8; stroke-width: 2.5px; }}
        ::-webkit-scrollbar {{ width: 6px; height: 6px; }}
        ::-webkit-scrollbar-track {{ background: #0f172a; }}
        ::-webkit-scrollbar-thumb {{ background: #334155; border-radius: 3px; }}
        ::-webkit-scrollbar-thumb:hover {{ background: #475569; }}
    </style>
</head>
<body class="h-screen overflow-hidden flex flex-col">

    <header class="bg-slate-900 border-b border-slate-800 px-6 py-4 flex justify-between items-center shadow-lg shrink-0">
        <div class="flex items-center gap-3">
            <div class="bg-sky-500/20 text-sky-400 p-2 rounded-lg border border-sky-500/30"><i class="fa-solid fa-dna text-xl"></i></div>
            <div>
                <h1 class="text-lg font-extrabold tracking-tight text-white flex items-center gap-2">
                    Genomic Ontology Hub
                </h1>
                <p class="text-xs text-slate-400 mt-0.5">Explore by System, Phenotype, or Gene</p>
            </div>
        </div>
        
        <div class="relative w-96 hidden md:block">
            <span class="absolute inset-y-0 left-0 flex items-center pl-3 text-slate-500"><i class="fa-solid fa-magnifying-glass"></i></span>
            <input type="text" id="globalSearch" class="w-full bg-slate-950 border border-slate-800 rounded-lg pl-10 pr-4 py-2 text-sm text-slate-200 focus:outline-none focus:border-sky-500" placeholder="Search systems, phenotypes, genes..." oninput="onGlobalSearch(this.value)">
        </div>
    </header>

    <main class="flex-1 flex overflow-hidden">
        
        <!-- Left Pane: D3 Collapsible Tree (3 Levels: L1 -> L2 -> Gene) -->
        <section class="w-1/3 flex flex-col border-r border-slate-800 bg-slate-950/40 relative min-w-[350px]">
            <div class="p-4 bg-slate-900/60 border-b border-slate-800 flex justify-between items-center z-10">
                <div>
                    <h2 class="text-sm font-bold text-slate-200 flex items-center gap-2">
                        <i class="fa-solid fa-diagram-project text-sky-400"></i> Domain Explorer
                    </h2>
                    <p class="text-[10px] text-slate-500 uppercase font-bold tracking-wider mt-1">System → Subcategory → Gene</p>
                </div>
            </div>
            <div class="flex-1 relative overflow-auto" id="canvasContainer">
                <svg id="hpoTreeSvg" class="w-full h-full absolute inset-0"></svg>
            </div>
        </section>

        <!-- Right Pane: Ontology Inspector -->
        <section class="w-2/3 flex flex-col bg-slate-900/20 overflow-hidden relative">
            <div class="p-6 bg-slate-900/40 border-b border-slate-800 shrink-0 flex items-center gap-4">
                <div class="w-14 h-14 rounded-xl bg-sky-500/10 text-sky-400 flex items-center justify-center border border-sky-500/20 text-3xl shadow-inner" id="headerIcon"><i class="fa-solid fa-sitemap"></i></div>
                <div>
                    <div class="text-[10px] uppercase font-bold tracking-widest text-sky-400 mb-1" id="headerType">ROOT ONTOLOGY</div>
                    <h3 class="text-2xl font-extrabold text-white tracking-tight" id="headerTitle">Select a Domain</h3>
                </div>
            </div>

            <div class="bg-slate-900/60 border-b border-slate-800 px-6 flex gap-6 shrink-0 z-10" id="tabsNavBar">
                <button onclick="switchTab('overview')" class="tab-btn py-3 text-sm font-bold text-slate-400 hover:text-white transition-all active" id="tab-overview"><i class="fa-solid fa-chart-pie"></i> Overview & Analysis</button>
                <button onclick="switchTab('variants')" class="tab-btn py-3 text-sm font-bold text-slate-400 hover:text-white transition-all" id="tab-variants"><i class="fa-solid fa-list-ul"></i> Genes & Variants <span class="bg-slate-800 text-slate-400 text-[10px] px-1.5 py-0.5 rounded-full ml-1" id="variantCountTag">0</span></button>
            </div>

            <div class="flex-1 overflow-y-auto p-6" id="tabContent">
                
                <div id="overviewContent" class="space-y-6">
                    <div class="grid grid-cols-3 gap-4 mb-6" id="statsGrid"></div>
                    
                    <div class="bg-slate-800/40 border border-slate-800/60 p-5 rounded-xl">
                        <h4 class="text-xs font-extrabold uppercase tracking-wider text-sky-400 mb-2"><i class="fa-solid fa-brain"></i> Domain Analysis</h4>
                        <p class="text-sm leading-relaxed text-slate-300 font-medium" id="domainAnalysisText"></p>
                    </div>

                    <div id="geneDescContainer" class="bg-slate-800/40 border border-slate-800/60 p-5 rounded-xl hidden">
                        <h4 class="text-xs font-extrabold uppercase tracking-wider text-sky-400 mb-2">NCBI Gene Summary</h4>
                        <p class="text-sm leading-relaxed text-slate-300 font-medium" id="ncbiSummary"></p>
                    </div>
                </div>

                <div id="variantsContent" class="hidden">
                    <div class="w-full overflow-x-auto bg-slate-800/30 border border-slate-700/60 rounded-lg shadow-lg">
                        <table class="w-full text-left text-sm text-slate-300">
                            <thead class="text-[10px] text-slate-400 uppercase bg-slate-900/80 border-b border-slate-700/80">
                                <tr>
                                    <th class="px-4 py-3">Gene / Description</th>
                                    <th class="px-4 py-3">Variant / Effect</th>
                                    <th class="px-4 py-3">ClinVar / Tier</th>
                                    <th class="px-4 py-3">CADD / REVEL</th>
                                    <th class="px-4 py-3">Depth / Qual</th>
                                </tr>
                            </thead>
                            <tbody id="variantTableBody" class="divide-y divide-slate-700/50"></tbody>
                        </table>
                    </div>
                </div>

            </div>
        </section>
    </main>

    <script>
        const monogenicFindings = {json.dumps(monogenic_findings)};
        let rootNode = null, svgGroup = null, treeLayout = null;
        const margin = {{top: 20, right: 120, bottom: 20, left: 120}};

        function buildHierarchy() {{
            const tree = {{ name: "All Patient Findings", type: "root", children: [], allGenes: monogenicFindings }};
            const l1Map = {{}};

            monogenicFindings.forEach(gene => {{
                const l1 = gene.domain_l1;
                const l2 = gene.domain_l2;
                
                if (!l1Map[l1.id]) {{
                    l1Map[l1.id] = {{ name: l1.label, type: "level1", color: l1.color, children: [], l2Map: {{}}, allGenes: [] }};
                    tree.children.push(l1Map[l1.id]);
                }}
                l1Map[l1.id].allGenes.push(gene);
                
                if (!l1Map[l1.id].l2Map[l2.id]) {{
                    l1Map[l1.id].l2Map[l2.id] = {{ name: l2.label, type: "level2", color: l2.color, children: [], allGenes: [] }};
                    l1Map[l1.id].children.push(l1Map[l1.id].l2Map[l2.id]);
                }}
                l1Map[l1.id].l2Map[l2.id].allGenes.push(gene);
                
                let tier = "Tier 3";
                if(gene.variants.length > 0) tier = gene.variants[0].tier;
                const gColor = tier === "Tier 1" ? "#ef4444" : (tier === "Tier 2" ? "#f59e0b" : "#a855f7");
                
                l1Map[l1.id].l2Map[l2.id].children.push({{ name: gene.gene_symbol, type: "gene", color: gColor, allGenes: [gene] }});
            }});
            return tree;
        }}

        function initializeD3Tree() {{
            const data = buildHierarchy();
            const container = document.getElementById('canvasContainer');
            const svg = d3.select("#hpoTreeSvg");
            
            svgGroup = svg.append("g").attr("transform", "translate(120,20)");
            const zoom = d3.zoom().scaleExtent([0.1, 4]).on("zoom", (e) => svgGroup.attr("transform", e.transform));
            svg.call(zoom).on("dblclick.zoom", null);

            rootNode = d3.hierarchy(data);
            rootNode.x0 = container.clientHeight / 2;
            rootNode.y0 = 0;

            rootNode.descendants().forEach(d => {{
                if (d.depth > 0) {{ d._children = d.children; d.children = null; }}
            }});

            treeLayout = d3.tree().nodeSize([30, 220]);
            updateTree(rootNode);
            inspectNode(rootNode.data);
        }}

        function updateTree(source) {{
            const treeData = treeLayout(rootNode);
            const nodes = treeData.descendants(), links = treeData.descendants().slice(1);
            nodes.forEach(d => d.y = d.depth * 200);

            const node = svgGroup.selectAll('g.node').data(nodes, d => d.id || (d.id = ++window.i));
            const nodeEnter = node.enter().append('g').attr('class', 'node')
                .attr("transform", d => "translate(" + source.y0 + "," + source.x0 + ")")
                .on('click', clickNode);

            nodeEnter.append('circle').attr('class', 'node-circle').attr('r', 1e-6)
                .style("fill", d => d._children ? "#1e293b" : d.data.color)
                .style("stroke", d => d.data.color);

            nodeEnter.append('text').attr("dy", ".35em")
                .attr("x", d => d.children || d._children ? -13 : 13)
                .attr("text-anchor", d => d.children || d._children ? "end" : "start")
                .text(d => d.data.name.length > 25 ? d.data.name.substring(0,23)+'...' : d.data.name)
                .attr("fill", "#e2e8f0").attr("font-size", "11px");

            const nodeUpdate = nodeEnter.merge(node);
            nodeUpdate.transition().duration(400).attr("transform", d => "translate(" + d.y + "," + d.x + ")");
            nodeUpdate.select('circle.node-circle').attr('r', 7)
                .style("fill", d => d._children ? "#1e293b" : d.data.color);

            const nodeExit = node.exit().transition().duration(400)
                .attr("transform", d => "translate(" + source.y + "," + source.x + ")").remove();
            nodeExit.select('circle').attr('r', 1e-6);

            const link = svgGroup.selectAll('path.node-link').data(links, d => d.id);
            const linkEnter = link.enter().insert('path', "g").attr("class", "node-link")
                .attr('d', d => "M " + source.y0 + " " + source.x0 + " C " + source.y0 + " " + source.x0 + ", " + source.y0 + " " + source.x0 + ", " + source.y0 + " " + source.x0);
                
            const linkUpdate = linkEnter.merge(link);
            linkUpdate.transition().duration(400)
                .attr('d', d => "M " + d.parent.y + " " + d.parent.x + " C " + (d.parent.y + d.y) / 2 + " " + d.parent.x + ", " + (d.parent.y + d.y) / 2 + " " + d.x + ", " + d.y + " " + d.x);

            link.exit().transition().duration(400)
                .attr('d', d => "M " + source.y + " " + source.x + " C " + source.y + " " + source.x + ", " + source.y + " " + source.x + ", " + source.y + " " + source.x).remove();
            nodes.forEach(d => {{ d.x0 = d.x; d.y0 = d.y; }});
        }}

        function clickNode(event, d) {{
            if (d.children) {{ d._children = d.children; d.children = null; }}
            else {{ d.children = d._children; d._children = null; }}
            updateTree(d);
            inspectNode(d.data);
        }}
        
        function inspectNode(nodeData) {{
            const genes = nodeData.allGenes || [];
            let allVariants = [];
            let tier1 = 0, tier2 = 0;
            genes.forEach(g => {{
                g.variants.forEach(v => {{
                    allVariants.push({{gene: g, variant: v}});
                    if(v.tier === 'Tier 1') tier1++;
                    if(v.tier === 'Tier 2') tier2++;
                }});
            }});

            document.getElementById("headerTitle").innerText = nodeData.name;
            document.getElementById("headerType").innerText = nodeData.type === 'root' ? "ALL FINDINGS" : (nodeData.type === 'gene' ? "GENE" : "ONTOLOGY DOMAIN");
            
            let tier1Cls = tier1 > 0 ? 'text-red-400' : 'text-white';
            let tier1Bg = tier1 > 0 ? 'border-red-500/50 bg-red-500/10' : 'border-slate-800';
            
            document.getElementById("statsGrid").innerHTML = 
                '<div class="bg-slate-900 border border-slate-800 p-4 rounded-lg text-center">' +
                    '<div class="text-3xl font-black text-white">' + genes.length + '</div>' +
                    '<div class="text-[10px] uppercase font-bold text-slate-400 mt-1">Total Genes</div>' +
                '</div>' +
                '<div class="bg-slate-900 border border-slate-800 p-4 rounded-lg text-center">' +
                    '<div class="text-3xl font-black text-white">' + allVariants.length + '</div>' +
                    '<div class="text-[10px] uppercase font-bold text-slate-400 mt-1">Total Variants</div>' +
                '</div>' +
                '<div class="bg-slate-900 border ' + tier1Bg + ' p-4 rounded-lg text-center">' +
                    '<div class="text-3xl font-black ' + tier1Cls + '">' + tier1 + '</div>' +
                    '<div class="text-[10px] uppercase font-bold text-slate-400 mt-1">Tier 1 Findings</div>' +
                '</div>';
            
            let analysis = "No variants found in this domain.";
            if (allVariants.length > 0) {{
                if (tier1 > 0) analysis = "<strong>High Risk Indicated:</strong> This domain contains " + tier1 + " pathogenic (Tier 1) finding(s). Immediate clinical correlation is recommended for the associated genes.";
                else if (tier2 > 0) analysis = "<strong>Monitoring Recommended:</strong> This domain contains " + tier2 + " Variant(s) of Uncertain Significance (Tier 2). Review clinical phenotypes to determine if these variants explain the patient's presentation.";
                else analysis = "<strong>Standard Risk:</strong> Found " + allVariants.length + " variants, mostly benign or Tier 3 regulatory variants. No major pathogenic drivers detected in this domain.";
            }}
            document.getElementById("domainAnalysisText").innerHTML = analysis;

            if (nodeData.type === 'gene') {{
                document.getElementById("geneDescContainer").classList.remove("hidden");
                document.getElementById("ncbiSummary").innerText = genes[0].ncbi_description;
            }} else {{
                document.getElementById("geneDescContainer").classList.add("hidden");
            }}

            document.getElementById("variantCountTag").innerText = allVariants.length;
            const vBody = document.getElementById("variantTableBody");
            let htmlBuf = "";

            allVariants.forEach((item, idx) => {{
                const g = item.gene;
                const v = item.variant;
                const sigCls = v.clinvar_significance.toLowerCase().includes("pathogenic") ? "text-red-400 border-red-500/30 bg-red-500/10" : "text-yellow-400 border-yellow-500/30 bg-yellow-500/10";
                
                let caddCls = "text-slate-300";
                if(v.cadd_phred !== "N/A" && parseFloat(v.cadd_phred) >= 20) caddCls = "text-amber-400 font-bold";
                if(v.cadd_phred !== "N/A" && parseFloat(v.cadd_phred) >= 30) caddCls = "text-red-400 font-bold";
                
                let revelCls = "text-slate-300";
                if(v.revel_score !== "N/A" && parseFloat(v.revel_score) >= 0.5) revelCls = "text-amber-400 font-bold";
                if(v.revel_score !== "N/A" && parseFloat(v.revel_score) >= 0.75) revelCls = "text-red-400 font-bold";

                htmlBuf += '<tr class="hover:bg-slate-800/40 transition-colors">' +
                    '<td class="px-4 py-4 w-1/4">' +
                        '<div class="font-bold text-white text-base">' + g.gene_symbol + '</div>' +
                        '<div class="text-[10px] text-slate-400 mt-1 line-clamp-2" title="' + g.ncbi_description + '">' + (g.ncbi_description || 'No description') + '</div>' +
                    '</td>' +
                    '<td class="px-4 py-4">' +
                        '<div class="font-mono font-bold text-sky-400">' + v.rsid + '</div>' +
                        '<div class="text-xs text-slate-300 mt-1">' + v.impact_consequence + '</div>' +
                        '<div class="text-[10px] text-slate-500 mt-0.5">' + v.genotype + '</div>' +
                    '</td>' +
                    '<td class="px-4 py-4">' +
                        '<div class="text-[10px] font-bold px-2 py-0.5 rounded uppercase border inline-block ' + sigCls + '">' + v.clinvar_significance + '</div>' +
                        '<div class="text-xs font-bold mt-1 text-slate-300">' + v.tier + '</div>' +
                    '</td>' +
                    '<td class="px-4 py-4">' +
                        '<div class="text-xs text-slate-400">CADD: <span class="' + caddCls + '">' + v.cadd_phred + '</span></div>' +
                        '<div class="text-xs text-slate-400 mt-1">REVEL: <span class="' + revelCls + '">' + v.revel_score + '</span></div>' +
                    '</td>' +
                    '<td class="px-4 py-4">' +
                        '<div class="text-xs text-slate-300 font-bold">' + v.read_depth + ' <span class="text-slate-500 font-normal ml-1">Reads</span></div>' +
                        '<div class="text-[10px] text-slate-500 mt-1">Qual: ' + v.read_quality + '</div>' +
                    '</td>' +
                '</tr>';
            }});
            vBody.innerHTML = htmlBuf;
        }}

        function switchTab(tabId) {{
            document.querySelectorAll('.tab-btn').forEach(btn => btn.classList.remove('active'));
            document.getElementById("tab-" + tabId).classList.add('active');
            ['overview', 'variants'].forEach(t => document.getElementById(t + "Content").classList.add("hidden"));
            document.getElementById(tabId + "Content").classList.remove("hidden");
        }}

        function onGlobalSearch(term) {{
            term = term.toLowerCase().trim();
            if (!term) return;
            
            let foundNode = null;
            rootNode.descendants().forEach(d => {{
                if (d.data.name.toLowerCase().includes(term)) {{
                    if(!foundNode) foundNode = d; 
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
                inspectNode(foundNode.data);
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
    print(f"Successfully generated Ontology Hub HTML at: {output_filepath}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input")
    parser.add_argument("-o", "--output", default="reports/visual_ontology_explorer.html")
    parser.add_argument("-d", "--demo", "--mock", action="store_true")
    args = parser.parse_args()

    out_dir = os.path.dirname(args.output)
    if out_dir: os.makedirs(out_dir, exist_ok=True)

    if args.demo:
        mel_path = Path(__file__).parent / "logs" / "mel_actionable.json"
        with open(mel_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        generate_upgraded_visual_report(data, args.output)
