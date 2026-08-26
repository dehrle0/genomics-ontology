#!/usr/bin/env python3
"""
Visual Ontology Explorer - Iteration 5: Claude-style Light UI & Deep Enrichment
Enriches missing data via MyGene.info and PubMed E-Utilities.
"""
import argparse, json, os, re, sys, yaml
import requests
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

CACHE_DIR = Path(__file__).parent / ".cache"
CACHE_DIR.mkdir(exist_ok=True)

def load_domain_registry():
    config_path = Path(__file__).parent / "config" / "ontology_domains.yaml"
    if config_path.exists():
        with open(config_path, "r", encoding="utf-8") as f:
            return yaml.safe_load(f).get("level1_systems", {})
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

def fetch_mygene_summary(hugo):
    cache_file = CACHE_DIR / f"{hugo}_summary.json"
    if cache_file.exists():
        with open(cache_file, "r") as f: return json.load(f).get("summary", "")
    try:
        res = requests.get(f"https://mygene.info/v3/query?q=symbol:{hugo}&fields=summary", timeout=5)
        if res.status_code == 200:
            hits = res.json().get("hits", [])
            summary = hits[0].get("summary", "") if hits and "summary" in hits[0] else ""
            with open(cache_file, "w") as f: json.dump({"summary": summary}, f)
            return summary
    except Exception:
        pass
    return ""

def fetch_pubmed_pubs(hugo):
    cache_file = CACHE_DIR / f"{hugo}_pubs.json"
    if cache_file.exists():
        with open(cache_file, "r") as f: return json.load(f).get("pubs", [])
    pubs = []
    try:
        url_search = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term={hugo}&retmode=json&retmax=3"
        search_res = requests.get(url_search, timeout=5).json()
        idlist = search_res.get("esearchresult", {}).get("idlist", [])
        if idlist:
            ids_str = ",".join(idlist)
            url_sum = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id={ids_str}&retmode=json"
            sum_res = requests.get(url_sum, timeout=5).json()
            for p_id in idlist:
                doc = sum_res.get("result", {}).get(p_id, {})
                title = doc.get("title", f"Publication {p_id}")
                pubdate = doc.get("pubdate", "Unknown Date")
                pubs.append({"id": p_id, "title": title, "date": pubdate, "link": f"https://pubmed.ncbi.nlm.nih.gov/{p_id}/"})
    except Exception:
        pass
    with open(cache_file, "w") as f: json.dump({"pubs": pubs}, f)
    return pubs

def enrich_gene(gene_obj):
    if not gene_obj["ncbi_description"] or gene_obj["ncbi_description"] == "No description":
        sum_text = fetch_mygene_summary(gene_obj["gene_symbol"])
        if sum_text: gene_obj["ncbi_description"] = sum_text
    
    gene_obj["publications"] = fetch_pubmed_pubs(gene_obj["gene_symbol"])
    return gene_obj

def build_data_model(raw_records, domain_reg):
    genes_map = {}
    for r in raw_records:
        hugo = r.get("hugo") or r.get("gene_symbol") or "Unknown"
        hpo_ids_raw = r.get("gene_hpo_id") or ""
        hpo_terms_raw = r.get("gene_hpo_term") or ""
        hpo_ids = [h.strip() for h in hpo_ids_raw.split(";") if h.strip()]
        hpo_terms = [h.strip() for h in hpo_terms_raw.split(";") if h.strip()]
        hpo_pairs = [{"id": hid, "label": hpo_terms[i] if i < len(hpo_terms) else hid} for i, hid in enumerate(hpo_ids)]

        ev = r.get("evidence") or {}
        hpo_ctx = ev.get("hpo_context") or []
        go_ctx = ev.get("go_context") or []

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
                "associated_hpos": hpo_pairs,
                "hpo_context": list(set(hpo_ctx)),
                "go_context": list(set(go_ctx)),
                "variants": [],
                "publications": []
            }
        else:
            genes_map[hugo]["hpo_context"] = list(set(genes_map[hugo]["hpo_context"] + hpo_ctx))
            genes_map[hugo]["go_context"] = list(set(genes_map[hugo]["go_context"] + go_ctx))
            
        zyg = r.get("zygosity") or ev.get("zygosity") or "Heterozygous"
        depth = ev.get("tot_reads") or r.get("tot_reads") or r.get("vcfinfo__tot_reads") or "N/A"
        alt_depth = ev.get("alt_reads") or r.get("alt_reads") or "N/A"
        quality = ev.get("qual") or r.get("phred") or r.get("qual") or "N/A"
        cadd = r.get("cadd_phred") or r.get("cadd") or "N/A"
        spliceai = get_max_spliceai(r)
        revel = r.get("revel") or r.get("revel_score") or "N/A"
        
        am_path = r.get("am_path") or "N/A"
        am_class = r.get("am_class") or "N/A"
        
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
            "read_depth_alt": alt_depth,
            "read_depth_total": depth,
            "read_quality": quality,
            "cadd_phred": cadd,
            "spliceai_max": spliceai,
            "am_path": am_path,
            "am_class": am_class,
            "revel_score": revel,
            "tier": tier,
            "gnomad_af": r.get("gnomad4_af") or r.get("gnomad_af") or "0.0",
            "allofus_af": r.get("allofus_af") or "0.0"
        })

    gene_list = list(genes_map.values())
    with ThreadPoolExecutor(max_workers=10) as executor:
        enriched_genes = list(executor.map(enrich_gene, gene_list))
    
    return enriched_genes

def generate_upgraded_visual_report(report_data: dict, output_filepath: str):
    domain_reg = load_domain_registry()
    raw_records = report_data.get("records") or report_data.get("monogenic_findings") or []
    monogenic_findings = build_data_model(raw_records, domain_reg)

    with open(Path(__file__).parent / "template_ontology_hub.html", "r", encoding="utf-8") as f:
        html_template = f.read()

    html_content = html_template.replace("{{MONOGENIC_FINDINGS_JSON}}", json.dumps(monogenic_findings))
    
    with open(output_filepath, "w", encoding="utf-8") as f:
        f.write(html_content)
    print(f"Successfully generated Ontology Hub HTML at: {output_filepath}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output", default="reports/visual_ontology_explorer.html")
    parser.add_argument("-d", "--demo", action="store_true")
    parser.add_argument("-i", "--input", help="Path to actionable JSON file")
    args = parser.parse_args()

    out_dir = os.path.dirname(args.output)
    if out_dir: os.makedirs(out_dir, exist_ok=True)

    input_path = args.input or (Path(__file__).parent / "logs" / "mel_actionable.json" if args.demo else None)
    if input_path and os.path.exists(input_path):
        with open(input_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        generate_upgraded_visual_report(data, args.output)
    else:
        print("Please provide --input <file.json> or --demo")
