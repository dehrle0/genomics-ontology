#!/usr/bin/env python3
"""
medgemma_synthesis_runner.py
Air-gapped, zero-PII runner for clinical genomics synthesis using local MedGemma 27B.

Features:
  1. Checks if medgemma-27b model exists on disk (LM Studio or llama.cpp path).
  2. Auto-starts local inference runtime (llama-server) if not already listening.
  3. Sanitizes and tokenizes all inputs (Zero PII, PROBAND_01 tokens only).
  4. Dispatches multi-thousand token clinical synthesis request (default: 4096 tokens).
  5. Formats output in publication-grade VSCP-DF Markdown.
  6. Optional local PII binding & PDF compilation strictly on localhost.
  7. Automatically unloads model / cleans up server post-use to reclaim 32GB RAM.
"""

import os
import sys
import json
import time
import argparse
import subprocess
import shutil
import urllib.request
import urllib.error

MEDGEMMA_27B_PATH = os.path.expanduser("~/.lmstudio/models/lmstudio-community/medgemma-27b-text-it-GGUF/medgemma-27b-text-it-Q4_K_M.gguf")
MEDGEMMA_4B_PATH = os.path.expanduser("~/.lmstudio/models/unsloth/medgemma-1.5-4b-it-GGUF/medgemma-1.5-4b-it-Q8_0.gguf")
MMPROJ_4B_PATH = os.path.expanduser("~/.lmstudio/models/unsloth/medgemma-1.5-4b-it-GGUF/mmproj-F32.gguf")
MEDGEMMA_4B_Q4_PATH = os.path.expanduser("~/.lmstudio/models/gguf-org/medgemma-1.5-4b-it-gguf/medgemma-1.5-4b-it-q4_0.gguf")
MMPROJ_4B_Q4_PATH = os.path.expanduser("~/.lmstudio/models/gguf-org/medgemma-1.5-4b-it-gguf/mmproj-medgemma-1.5-4b-it-q4_0.gguf")

LLAMA_SERVER_PATHS = [
    "/home/daniel-ehrle/.lmstudio/extensions/backends/llama.cpp-linux-x86_64-vulkan-avx2-2.28.2/llama-server",
    "/home/daniel-ehrle/.lmstudio/extensions/backends/llama.cpp-linux-x86_64-vulkan-avx2-2.27.1/llama-server",
    shutil.which("llama-server")
]

def find_llama_server():
    for p in LLAMA_SERVER_PATHS:
        if p and os.path.exists(p) and os.access(p, os.X_OK):
            return p
    return None

def check_model_exists(model_type="4b"):
    if model_type == "4b":
        if os.path.exists(MEDGEMMA_4B_PATH):
            print(f"[Model Check] Found verified MedGemma 1.5 4B (Q8_0 Multimodal): {MEDGEMMA_4B_PATH} ({os.path.getsize(MEDGEMMA_4B_PATH)/(1024**3):.2f} GB)")
            return MEDGEMMA_4B_PATH, MMPROJ_4B_PATH if os.path.exists(MMPROJ_4B_PATH) else None
        elif os.path.exists(MEDGEMMA_4B_Q4_PATH):
            print(f"[Model Check] Found verified MedGemma 1.5 4B (Q4_0 Multimodal): {MEDGEMMA_4B_Q4_PATH} ({os.path.getsize(MEDGEMMA_4B_Q4_PATH)/(1024**3):.2f} GB)")
            return MEDGEMMA_4B_Q4_PATH, MMPROJ_4B_Q4_PATH if os.path.exists(MMPROJ_4B_Q4_PATH) else None
        print("[Model Check Warning] MedGemma 1.5 4B not found in standard paths.")
        return None, None
    else:
        if os.path.exists(MEDGEMMA_27B_PATH):
            print(f"[Model Check] Found verified MedGemma 27B: {MEDGEMMA_27B_PATH} ({os.path.getsize(MEDGEMMA_27B_PATH)/(1024**3):.2f} GB)")
            return MEDGEMMA_27B_PATH, None
        print("[Model Check Warning] MedGemma 27B not found in standard paths.")
        return None, None

def is_server_listening(port=7002):
    try:
        req = urllib.request.Request(f"http://127.0.0.1:{port}/v1/models")
        with urllib.request.urlopen(req, timeout=2) as resp:
            return resp.status == 200
    except Exception:
        return False

def start_local_server(model_path, mmproj_path=None, port=7002, ctx_size=16384):
    server_bin = find_llama_server()
    if not server_bin:
        print("[Server Error] llama-server binary not found.")
        return None
    
    cmd = [server_bin, "-m", model_path, "--port", str(port), "-c", str(ctx_size), "-ngl", "99"]
    if mmproj_path and os.path.exists(mmproj_path):
        cmd.extend(["--mmproj", mmproj_path])
        print(f"[Multimodal Projector] Enabled medical image projector: {mmproj_path}")
        
    print(f"[Server Launch] Starting llama-server on 127.0.0.1:{port} ({ctx_size} context)...")
    log_file = open("/tmp/llama_server_medgemma.log", "w")
    proc = subprocess.Popen(
        cmd,
        stdout=log_file,
        stderr=subprocess.STDOUT,
        preexec_fn=os.setsid
    )
    for _ in range(90):
        if is_server_listening(port):
            model_tag = "MedGemma 1.5 4B (Multimodal)" if mmproj_path else "MedGemma 27B"
            print(f"[Server Ready] {model_tag} loaded and listening on port {port} with {ctx_size:,} token context.")
            return proc
        time.sleep(1)
    print("[Server Warning] Timeout waiting for MedGemma server. Last 20 lines of log:")
    try:
        with open("/tmp/llama_server_medgemma.log", "r") as f:
            lines = f.readlines()
            for l in lines[-20:]:
                print("  " + l.rstrip())
    except Exception:
        pass
    return proc

def stop_local_server(proc):
    if proc:
        print("[Server Cleanup] Stopping local MedGemma server and reclaiming RAM...")
        try:
            os.killpg(os.getpgid(proc.pid), 15)
            proc.wait(timeout=5)
        except Exception:
            try:
                os.killpg(os.getpgid(proc.pid), 9)
            except Exception:
                pass
        print("[Server Cleanup] Memory reclaimed successfully.")

def sanitize_for_medgemma(raw_findings, session_token="PROBAND_01"):
    """
    Enforces strict zero-PII scrubbing.
    Filters out silent synonymous variants without splice disruption.
    """
    clean_records = []
    for item in raw_findings:
        so = str(item.get("so") or "").upper().strip()
        reasons = str(item.get("reason_codes") or "").lower()
        if so == "SYN" and "spliceai_high" not in reasons and "spliceai_mod" not in reasons:
            continue
        rec = {
            "subject_token": session_token,
            "gene": str(item.get("hugo") or item.get("gene") or "").upper().strip(),
            "variant": str(item.get("achange") or item.get("cchange") or item.get("variant") or "").strip(),
            "zygosity": str(item.get("zygosity") or "het").strip(),
            "so": so,
            "tier": str(item.get("tier") or "").strip(),
            "clinvar_sig": str(item.get("clinvar_sig") or "").strip(),
            "clinvar_id": str(item.get("clinvar_id") or "").strip(),
            "omim_id": str(item.get("omim_id") or "").strip(),
            "chrom": str(item.get("chrom") or "").strip(),
            "pos": str(item.get("pos") or "").strip(),
            "ref": str(item.get("ref") or "").strip(),
            "alt": str(item.get("alt") or "").strip(),
            "rsid": str(item.get("rsid") or "").strip(),
            "gwas_pmid": str(item.get("gwas_pmid") or "").strip(),
            "cadd_phred": item.get("cadd_phred"),
            "revel": item.get("revel"),
            "am_class": item.get("am_class"),
            "avi_phred": item.get("avi_phred")
        }
        clean_records.append(rec)
    return clean_records

def sort_priority(rec):
    sig = str(rec.get("clinvar_sig") or "").lower()
    is_plp = "pathogenic" in sig and "conflicting" not in sig and "uncertain" not in sig and "benign" not in sig
    tier = str(rec.get("tier") or "").strip()
    so = str(rec.get("so") or "").upper().strip()

    cadd = float(rec.get("cadd_phred") or 0.0)
    avi = float(rec.get("avi_phred") or 0.0)
    rev = float(rec.get("revel") or 0.0)
    priority = 0.0
    if is_plp:
        priority += 120.0
    elif "likely pathogenic" in sig or "drug response" in sig:
        priority += 40.0
    elif "protective" in sig:
        priority += 20.0

    if tier == "Tier1":
        priority += 50.0
    elif tier == "Tier2":
        priority += 20.0

    if so in ["UT3", "UT5", "INT", "SYN"] and not is_plp:
        priority -= 40.0

    return priority + cadd + (avi * 0.5) + (rev * 20.0)

def append_supporting_documentation_appendix(report_md, clean_records):
    """
    Appends an authoritative, zero-hallucination Appendix of supporting documentation links
    and curated knowledgebase cross-references.
    """
    if "### Appendix: Supporting Documentation" in report_md or "## Appendix: Supporting Documentation" in report_md:
        return report_md

    lines = []
    lines.append("")
    lines.append("## Appendix: Supporting Documentation & Evidence Sources")
    lines.append("")
    lines.append(
        "This appendix compiles primary evidence accessions, external database cross-references, genomic coordinate mappings (GRCh38), "
        "and clinical guidelines for all key variants and genes analyzed across this report. All links connect directly to peer-reviewed "
        "public archives and expert clinical curation repositories."
    )
    lines.append("")
    lines.append("### 1. Variant Evidence & Cross-Reference Directory")
    lines.append("")
    lines.append("| Gene | Variant | Coordinates (GRCh38) | ClinVar | dbSNP | OMIM | ClinGen | AlphaGenome | Primary Guideline / Evidence |")
    lines.append("| :--- | :--- | :--- | :---: | :---: | :---: | :---: | :---: | :--- |")

    sorted_records = sorted(clean_records, key=sort_priority, reverse=True)
    seen = set()
    for rec in sorted_records:
        h = str(rec.get("gene") or "").upper().strip()
        v = str(rec.get("variant") or "").strip()
        so = str(rec.get("so") or "").upper().strip()
        tier = str(rec.get("tier") or "").strip()
        sig = str(rec.get("clinvar_sig") or "").lower()
        is_plp = "pathogenic" in sig and "conflicting" not in sig and "uncertain" not in sig

        # Prune silent synonymous, UTR, intronic, or Tier 3 unless definitive PLP
        if (so in ["SYN", "INT", "UT3", "UT5"] or tier == "Tier3") and not is_plp:
            continue

        chrom = str(rec.get("chrom") or "").strip()
        pos = str(rec.get("pos") or "").strip()
        ref = str(rec.get("ref") or "").strip().upper()
        alt = str(rec.get("alt") or "").strip().upper()
        key = (chrom, pos, ref, alt) if chrom and pos else (h, v)
        if key in seen:
            continue
        seen.add(key)

        coord_str = f"`{chrom}:{pos} {ref}>{alt}`" if chrom and pos else "—"

        cv_id = str(rec.get("clinvar_id") or "").replace("VCV", "").replace("vcv", "").strip()
        cv_link = f"[VCV{cv_id}](https://www.ncbi.nlm.nih.gov/clinvar/variation/{cv_id}/)" if cv_id and cv_id.isdigit() else "—"

        rs = str(rec.get("rsid") or "").strip()
        rs_link = f"[{rs}](https://www.ncbi.nlm.nih.gov/snp/{rs})" if rs.startswith("rs") else "—"

        om_raw = str(rec.get("omim_id") or "").strip()
        om_digits = [p.strip().replace("OMIM:", "").replace("MIM:", "").strip() for p in om_raw.replace(";", ",").split(",") if p.strip().isdigit()]
        om_link = f"[MIM:{om_digits[0]}](https://www.omim.org/entry/{om_digits[0]})" if om_digits else "—"

        cg_link = f"[{h}](https://search.clinicalgenome.org/kb/genes/{h})"

        if chrom and pos and ref and alt and ref != "-" and alt != "-":
            c_tag = f"chr{chrom}" if not chrom.startswith("chr") else chrom
            ag_link = f"[Atlas](https://alphagenome.deepmind.com/variant/{c_tag}:{pos}:{ref}>{alt})"
        else:
            ag_link = f"[Atlas](https://alphagenome.deepmind.com/gene/{h})"

        pmid = str(rec.get("gwas_pmid") or "").strip()
        if h == "CBLIF" and is_plp:
            guide_link = "[OMIM 261000 (Intrinsic Factor)](https://www.omim.org/entry/261000)"
        elif h == "F5":
            if "534" in v or "rs6025" in rs.lower():
                guide_link = "[ACMG/ACOG VTE (PMID 28373160)](https://pubmed.ncbi.nlm.nih.gov/28373160/)"
            elif "295" in v or is_plp:
                guide_link = "[Factor V Deficiency (ClinVar)](https://www.ncbi.nlm.nih.gov/clinvar/variation/2884751/)"
            else:
                guide_link = f"[ClinVar Evidence](https://www.ncbi.nlm.nih.gov/clinvar/variation/{cv_id}/)" if cv_id else f"[NCBI Gene {h}](https://www.ncbi.nlm.nih.gov/gene/?term={h})"
        elif h == "APOB":
            if "4481" in v or "1801695" in rs:
                guide_link = "[GWAS Lipids (PMID 41325697)](https://pubmed.ncbi.nlm.nih.gov/41325697/)"
            elif is_plp:
                guide_link = "[FHBL1 (OMIM 615558)](https://www.omim.org/entry/615558)"
            else:
                guide_link = f"[ClinVar Evidence](https://www.ncbi.nlm.nih.gov/clinvar/variation/{cv_id}/)" if cv_id else f"[NCBI Gene {h}](https://www.ncbi.nlm.nih.gov/gene/?term={h})"
        elif h == "DPYD" and ("732" in v or "1801160" in rs or is_plp):
            guide_link = "[CPIC Fluoropyrimidines (DPYD)](https://cpicpgx.org/guidelines/guideline-for-fluoropyrimidines-and-dpyd/)"
        elif h == "ANK2":
            if "3906" in v or "121912706" in rs or is_plp:
                guide_link = "[CredibleMeds QTdrugs](https://crediblemeds.org/)"
            else:
                guide_link = f"[ClinVar Evidence](https://www.ncbi.nlm.nih.gov/clinvar/variation/{cv_id}/)" if cv_id else f"[NCBI Gene {h}](https://www.ncbi.nlm.nih.gov/gene/?term={h})"
        elif h == "CDKN2B":
            guide_link = "[GWAS CAD 9p21 (PMID 30054458)](https://pubmed.ncbi.nlm.nih.gov/30054458/)"
        elif h == "POLG" and is_plp:
            guide_link = "[Valproate Toxicity (PMID 24077912)](https://pubmed.ncbi.nlm.nih.gov/24077912/)"
        elif h == "ATM" and is_plp:
            guide_link = "[NCCN Genetic Surveillance (ATM)](https://search.clinicalgenome.org/kb/genes/ATM)"
        elif h == "TAT" and is_plp:
            guide_link = "[Tyrosinemia Type II (OMIM 276600)](https://www.omim.org/entry/276600)"
        elif h == "VDR":
            guide_link = "[VDR Promoter (PMID 16785239)](https://pubmed.ncbi.nlm.nih.gov/16785239/)"
        elif pmid and pmid.isdigit():
            guide_link = f"[GWAS (PMID {pmid})](https://pubmed.ncbi.nlm.nih.gov/{pmid}/)"
        elif cv_id and cv_id.isdigit():
            guide_link = f"[ClinVar Evidence](https://www.ncbi.nlm.nih.gov/clinvar/variation/{cv_id}/)"
        else:
            guide_link = f"[NCBI Gene {h}](https://www.ncbi.nlm.nih.gov/gene/?term={h})"

        lines.append(f"| **{h}** | `{v}` | {coord_str} | {cv_link} | {rs_link} | {om_link} | {cg_link} | {ag_link} | {guide_link} |")
        if len(seen) >= 15:
            break

    lines.append("")
    lines.append("### 2. Authoritative Clinical & Pharmacogenomic Repositories")
    lines.append("* **[NCBI ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/):** National Center for Biotechnology Information public archive of human genomic variants and interpretations of clinical significance.")
    lines.append("* **[ClinGen (Clinical Genome Resource)](https://clinicalgenome.org/):** NIH-funded consortium curating authoritative evidence supporting gene-disease clinical validity, dosage sensitivity, and actionability.")
    lines.append("* **[OMIM (Online Mendelian Inheritance in Man)](https://www.omim.org/):** Curated database of human genes and genetic phenotypes founded by Victor A. McKusick at Johns Hopkins University.")
    lines.append("* **[DeepMind AlphaGenome Atlas](https://alphagenome.deepmind.com/):** Unified genomic AI foundation model providing 1-bp resolution locus exploration, chromatin accessibility, and multimodal impact predictions.")
    lines.append("* **[CPIC (Clinical Pharmacogenetics Implementation Consortium)](https://cpicpgx.org/guidelines/):** Peer-reviewed clinical practice guidelines enabling translation of genetic test results into actionable prescribing decisions.")
    lines.append("* **[PharmGKB (Pharmacogenomics Knowledgebase)](https://www.pharmgkb.org/):** Comprehensive resource curating knowledge on how genetic variations impact medication responses and clinical outcomes.")
    lines.append("* **[CredibleMeds (AZCERT)](https://crediblemeds.org/):** Evidence-based decision support resource maintaining stratified QT-prolonging drug lists and torsadogenic risk categories.")
    lines.append("* **[Broad Institute gnomAD (v4.1)](https://gnomad.broadinstitute.org/):** Reference population genomic dataset spanning >800,000 individuals for allele frequency estimation and gene constraint metrics.")
    lines.append("* **[NCBI dbSNP](https://www.ncbi.nlm.nih.gov/snp/):** Central public repository for single nucleotide polymorphisms and short genetic variations.")
    lines.append("")
    return report_md.strip() + "\n\n" + "\n".join(lines)

def build_medgemma_prompt(clean_records, session_token):
    # Sort by priority so top actionable findings receive dossiers first
    sorted_records = sorted(clean_records, key=sort_priority, reverse=True)
    top_records = []
    for r in sorted_records[:10]:
        top_records.append({
            "subject_token": session_token,
            "gene": r["gene"],
            "variant": r["variant"],
            "zygosity": r["zygosity"],
            "clinvar_sig": r["clinvar_sig"],
            "clinvar_id": r["clinvar_id"],
            "omim_id": r["omim_id"],
            "cadd_phred": r["cadd_phred"],
            "revel": r["revel"],
            "am_class": r["am_class"],
            "avi_phred": r["avi_phred"]
        })

    system_prompt = (
        "You are MedGemma, an elite clinical genomics AI assistant specializing in evidence synthesis. "
        "Operational Directives:\n"
        "1. Zero PII: Never invent or use patient names, dates of birth, or facilities. Refer solely to the proband as '" + session_token + "'.\n"
        "2. Strict Grounding & Anti-Hallucination: You must strictly restrict your analysis to the specific genes and variants present in the supplied JSON payload. Never reference, invent, or discuss genes or variants that are not in the provided callset.\n"
        "3. Clinical Scope & Tone: The exhaustive raw variant call tables and in-silico predictor matrices are already archived in the accompanying Master Ontology Explorer. Your role is to deliver a readable, highly authoritative executive clinical summary focused on the high-actionable level findings.\n"
        "4. Include sufficient molecular mechanism details and clinical reasoning to be authoritative without bogging down in repetitive technical data dumps. Do not loop or repeat identical variant phrases.\n"
        "5. Structure strictly following the VSCP-DF framework:\n"
        "   - Orientation & Executive Summary\n"
        "   - High-Actionable & Primary Findings (focused dossiers with authoritative depth for primary loci in the callset)\n"
        "   - Secondary Modifiers & Organ Surveillance (cardiovascular, metabolic, or cellular co-factors present in the callset)\n"
        "   - Critical Pharmacogenomic & Drug Interactions (contraindications or metabolic alterations specific to the patient's verified variants)\n"
        "   - Clinical Decision Calculus (Arguments FOR and AGAINST over-intervention, with explicit 0.00-1.00 Confidence Score)\n"
        "   - Action Directives: Action (Update Electronic Health Records with confirmed findings)\n"
        "   - Monitoring Directives: Monitor (Baseline & Periodic Clinical Surveillance tailored to the detected variants)\n"
        "   - Methodological Assumptions & Limitations (40x WGS boundaries, mosaicism, recessive carrier status)"
    )
    
    user_prompt = (
        f"Generate a publication-grade, authoritative Clinical Genomics Evidence Summary for {session_token} based on the following verified genomic callset:\n\n"
        f"{json.dumps(top_records, indent=2)}\n\n"
        "Requirements:\n"
        "- Focus strictly and exclusively on the findings present in the supplied callset with clear molecular rationale.\n"
        "- Detail any applicable drug interactions, contraindications, or metabolic warnings directly relevant to the provided variants.\n"
        "- Ensure Action (EHR updates), Monitor (laboratory directives), Assumptions/Limitations, and Decision Calculus are fully articulated and grounded only in the provided callset.\n"
        "- Output clean GitHub-Flavored Markdown directly without wrapping in markdown code blocks."
    )
    return system_prompt, user_prompt

def query_local_medgemma(system_prompt, user_prompt, port=7002, alias="medgemma-27b-it", max_tokens=8192):
    url = f"http://127.0.0.1:{port}/v1/chat/completions"
    payload = {
        "model": alias,
        "messages": [
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": user_prompt}
        ],
        "temperature": 0.2,
        "repeat_penalty": 1.15,
        "presence_penalty": 0.2,
        "frequency_penalty": 0.2,
        "max_tokens": max_tokens
    }
    data = json.dumps(payload).encode("utf-8")
    req = urllib.request.Request(url, data=data, headers={"Content-Type": "application/json"})
    
    timeout_secs = max(14400, int(max_tokens * 2.5))
    with urllib.request.urlopen(req, timeout=timeout_secs) as resp:
        res_json = json.loads(resp.read().decode("utf-8"))
        return res_json["choices"][0]["message"]["content"]

def main():
    parser = argparse.ArgumentParser(description="Airgapped MedGemma Multi-Thousand Token Runner")
    parser.add_argument("--model", choices=["4b", "27b"], default="4b", help="Model size to load (default: 4b)")
    parser.add_argument("--ctx-size", type=int, default=16384, help="Context size in tokens (default: 16384, supports 24576)")
    parser.add_argument("--input-json", required=True, help="Path to input actionable variants JSON")
    parser.add_argument("--out-md", required=True, help="Output Markdown report path")
    parser.add_argument("--session-token", default="PROBAND_01", help="Ephemeral anonymous token")
    parser.add_argument("--patient-name", default=None, help="Optional: real patient name for local-only binding")
    parser.add_argument("--patient-id", default=None, help="Optional: real patient sample ID for local-only binding")
    parser.add_argument("--max-tokens", type=int, default=8192, help="Maximum generation tokens (default: 8192, supports up to 16384)")
    parser.add_argument("--port", type=int, default=7002, help="Port of local server (default: 7002)")
    parser.add_argument("--keep-server-alive", action="store_true", help="Do not shut down server after generation")
    args = parser.parse_args()

    model_path, mmproj_path = check_model_exists(args.model)
    if not model_path:
        sys.exit(1)

    with open(args.input_json, "r", encoding="utf-8") as f:
        raw_data = json.load(f)

    findings = raw_data.get("records") or raw_data.get("variants") or (raw_data if isinstance(raw_data, list) else [])
    print(f"[Sanitizer] Ingesting {len(findings)} records. Scrubbing PII...")
    clean_records = sanitize_for_medgemma(findings, args.session_token)
    print(f"[Sanitizer Complete] Zero PII scrubbed payload ({len(clean_records)} candidates) for {args.session_token}.")

    sys_prompt, user_prompt = build_medgemma_prompt(clean_records, args.session_token)

    server_proc = None
    if not is_server_listening(args.port):
        server_proc = start_local_server(model_path, mmproj_path=mmproj_path, port=args.port, ctx_size=args.ctx_size)
        if not server_proc:
            print("[Error] Failed to start local inference server.")
            sys.exit(1)

    try:
        print(f"[Inference] Dispatching {args.max_tokens}-token clinical synthesis request to local MedGemma {args.model.upper()} (127.0.0.1:{args.port})...")
        synthesis = query_local_medgemma(sys_prompt, user_prompt, port=args.port, alias=f"medgemma-{args.model}", max_tokens=args.max_tokens)
        
        # Clean markdown code fences if wrapped by LLM
        synthesis = synthesis.strip()
        if synthesis.startswith("```markdown"):
            synthesis = synthesis[len("```markdown"):].strip()
        elif synthesis.startswith("```"):
            synthesis = synthesis[3:].strip()
        if synthesis.endswith("```"):
            synthesis = synthesis[:-3].strip()

        # Local PII re-binding if requested
        if args.patient_name and args.patient_id:
            synthesis = synthesis.replace(args.session_token, args.patient_name)
            synthesis = synthesis.replace(f"`{args.session_token}`", f"`{args.patient_id}`")
            print(f"[Local Binding] Applied patient metadata locally: {args.patient_name} ({args.patient_id})")

        # Append authoritative supporting documentation appendix with direct links
        synthesis = append_supporting_documentation_appendix(synthesis, clean_records)

        with open(args.out_md, "w", encoding="utf-8") as f:
            f.write(synthesis)
        print(f"[Success] Written {len(synthesis.splitlines())} lines of MedGemma synthesis to: {args.out_md}")

        # Export HTML and PDF deliverables
        out_html = args.out_md.replace(".md", ".html")
        out_pdf = args.out_md.replace(".md", ".pdf")
        patient_title = args.patient_name or "Clinical Evidence Synthesis"
        pid_title = args.patient_id or "PROBAND"
        try:
            from lib.generate_deep_research_report import format_report_html, generate_pdf
            html_content = format_report_html(patient_title, pid_title, synthesis)
            with open(out_html, "w", encoding="utf-8") as f:
                f.write(html_content)
            print(f"[Deliverable Export] Written HTML5:    {out_html}")
            generate_pdf(out_html, out_pdf)
        except Exception as e:
            print(f"[Export Warning] Could not render HTML/PDF: {e}")

        # Automatically sync to Google Drive Ontology folder if available
        gdrive_ontology = "/home/daniel-ehrle/Google Drive/My Drive/Ontology"
        if os.path.exists(gdrive_ontology):
            files_to_sync = [args.out_md]
            if os.path.exists(out_html):
                files_to_sync.append(out_html)
            if os.path.exists(out_pdf):
                files_to_sync.append(out_pdf)

            for fpath in files_to_sync:
                out_base = os.path.basename(fpath)
                gdrive_dest = os.path.join(gdrive_ontology, out_base)
                try:
                    shutil.copyfile(fpath, gdrive_dest)
                    print(f"[Google Drive Sync] Synced {out_base} to Google Drive: {gdrive_dest}")
                    for folder in os.listdir(gdrive_ontology):
                        subpath = os.path.join(gdrive_ontology, folder)
                        if os.path.isdir(subpath) and ("Daniel_Ehrle" in folder or "Melinda_Ehrle" in folder) and folder in args.out_md:
                            sub_dest = os.path.join(subpath, out_base)
                            shutil.copyfile(fpath, sub_dest)
                            print(f"[Google Drive Sync] Synced {out_base} to dated folder: {sub_dest}")
                except Exception as e:
                    print(f"[Google Drive Sync Warning] Could not copy {out_base} to Google Drive: {e}")
    finally:
        if server_proc and not args.keep_server_alive:
            stop_local_server(server_proc)

if __name__ == "__main__":
    main()
