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

LM_STUDIO_MODEL_PATH = os.path.expanduser("~/.lmstudio/models/lmstudio-community/medgemma-27b-text-it-GGUF/medgemma-27b-text-it-Q4_K_M.gguf")
LLAMA_MODEL_PATH = os.path.expanduser("~/ai-infrastructure/models/medgemma-27b-it.gguf")
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

def check_model_exists():
    paths = [LM_STUDIO_MODEL_PATH, LLAMA_MODEL_PATH]
    for p in paths:
        if os.path.exists(p) and os.path.getsize(p) > 10 * 1024 * 1024 * 1024:
            print(f"[Model Check] Found verified MedGemma 27B: {p} ({os.path.getsize(p)/(1024**3):.2f} GB)")
            return p
    print("[Model Check Warning] MedGemma 27B not found in standard paths.")
    return None

def is_server_listening(port=7002):
    try:
        req = urllib.request.Request(f"http://127.0.0.1:{port}/v1/models")
        with urllib.request.urlopen(req, timeout=2) as resp:
            return resp.status == 200
    except Exception:
        return False

def start_local_server(model_path, port=7002):
    server_bin = find_llama_server()
    if not server_bin:
        print("[Server Error] llama-server binary not found.")
        return None
    print(f"[Server Launch] Starting llama-server with Vulkan/AVX2 on 127.0.0.1:{port}...")
    proc = subprocess.Popen(
        [server_bin, "-m", model_path, "--port", str(port), "-c", "8192", "-ngl", "99"],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        preexec_fn=os.setsid
    )
    for _ in range(30):
        if is_server_listening(port):
            print(f"[Server Ready] MedGemma 27B loaded and listening on port {port}.")
            return proc
        time.sleep(1)
    print("[Server Warning] Timeout waiting for MedGemma 27B server.")
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
            "clinvar_sig": str(item.get("clinvar_sig") or "").strip(),
            "clinvar_id": str(item.get("clinvar_id") or "").strip(),
            "omim_id": str(item.get("omim_id") or "").strip(),
            "cadd_phred": item.get("cadd_phred"),
            "revel": item.get("revel"),
            "am_class": item.get("am_class"),
            "avi_phred": item.get("avi_phred")
        }
        clean_records.append(rec)
    return clean_records

def sort_priority(rec):
    sig = str(rec.get("clinvar_sig") or "").lower()
    is_plp = "pathogenic" in sig and "conflicting" not in sig and "uncertain" not in sig
    cadd = float(rec.get("cadd_phred") or 0.0)
    avi = float(rec.get("avi_phred") or 0.0)
    hugo = str(rec.get("gene") or "").upper()
    priority = 0.0
    if is_plp: priority += 100.0
    if hugo in ["CBLIF", "F5", "APOB", "DPYD", "CDKN2B", "POLG"]: priority += 50.0
    return priority + cadd + (avi * 0.5)

def build_medgemma_prompt(clean_records, session_token):
    # Sort by priority so top actionable findings receive dossiers first
    sorted_records = sorted(clean_records, key=sort_priority, reverse=True)
    top_records = sorted_records[:5]

    system_prompt = (
        "You are MedGemma, an elite clinical genomics AI assistant specializing in evidence synthesis. "
        "Strict Operational Directives:\n"
        "1. Zero PII: Never invent or use patient names, dates of birth, or facilities. Refer solely to the proband as '" + session_token + "'.\n"
        "2. Zero Extrapolation: Ground all clinical evaluations strictly in the provided variant data and peer-reviewed ACMG/ClinVar/CPIC evidence.\n"
        "3. Output Format: Produce a deeply researched, multi-page clinical genomics evidence dossier in clean GitHub-Flavored Markdown adhering strictly to the VSCP-DF framework:\n"
        "   - Orientation (What We Are Covering)\n"
        "   - Body: Comprehensive Evidence Dossiers for the prioritized findings (molecular mechanism, ClinVar/OMIM disease correlation, multi-engine in silico consensus, ACMG pathogenicity calculus, actionable surveillance, and drug contraindications)\n"
        "   - Conclusions: Diagnostic & Decision Calculus (Arguments FOR and AGAINST clinical intervention; explicit Confidence Score 0.00 to 1.00)\n"
        "   - Opportunities: High-Yield Clinical Next Steps"
    )
    
    user_prompt = (
        f"Generate a publication-grade, deeply researched Clinical Genomics Evidence Synthesis for {session_token} based on the following verified genomic callset:\n\n"
        f"{json.dumps(top_records, indent=2)}\n\n"
        "Requirements:\n"
        "- Deliver complete narrative dossiers for all prioritized loci without cutting off.\n"
        "- Disallow monogenic disease claims for common alleles (gnomAD AF > 1%) or benign in-silico variants.\n"
        "- Provide explicit laboratory surveillance and pharmacogenomic contraindications.\n"
        "- Ensure Conclusions (FOR/AGAINST over-intervention with Confidence Score) and Opportunities are fully written.\n"
        "- Output clean Markdown directly without surrounding code blocks."
    )
    return system_prompt, user_prompt

def query_local_medgemma(system_prompt, user_prompt, port=7002, alias="medgemma-27b-it", max_tokens=4096):
    url = f"http://127.0.0.1:{port}/v1/chat/completions"
    payload = {
        "model": alias,
        "messages": [
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": user_prompt}
        ],
        "temperature": 0.15,
        "max_tokens": max_tokens
    }
    data = json.dumps(payload).encode("utf-8")
    req = urllib.request.Request(url, data=data, headers={"Content-Type": "application/json"})
    
    timeout_secs = max(3600, int(max_tokens * 1.5))
    with urllib.request.urlopen(req, timeout=timeout_secs) as resp:
        res_json = json.loads(resp.read().decode("utf-8"))
        return res_json["choices"][0]["message"]["content"]

def main():
    parser = argparse.ArgumentParser(description="Airgapped MedGemma Multi-Thousand Token Runner")
    parser.add_argument("--input-json", required=True, help="Path to input actionable variants JSON")
    parser.add_argument("--out-md", required=True, help="Output Markdown report path")
    parser.add_argument("--session-token", default="PROBAND_01", help="Ephemeral anonymous token")
    parser.add_argument("--patient-name", default=None, help="Optional: real patient name for local-only binding")
    parser.add_argument("--patient-id", default=None, help="Optional: real patient sample ID for local-only binding")
    parser.add_argument("--max-tokens", type=int, default=4096, help="Maximum generation tokens (default: 4096)")
    parser.add_argument("--port", type=int, default=7002, help="Port of local server (default: 7002)")
    parser.add_argument("--keep-server-alive", action="store_true", help="Do not shut down server after generation")
    args = parser.parse_args()

    model_path = check_model_exists()
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
        server_proc = start_local_server(model_path, port=args.port)
        if not server_proc:
            print("[Error] Failed to start local inference server.")
            sys.exit(1)

    try:
        print(f"[Inference] Dispatching {args.max_tokens}-token clinical synthesis request to local MedGemma (127.0.0.1:{args.port})...")
        synthesis = query_local_medgemma(sys_prompt, user_prompt, port=args.port, max_tokens=args.max_tokens)
        
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

        with open(args.out_md, "w", encoding="utf-8") as f:
            f.write(synthesis)
        print(f"[Success] Written {len(synthesis.splitlines())} lines of MedGemma synthesis to: {args.out_md}")
    finally:
        if server_proc and not args.keep_server_alive:
            stop_local_server(server_proc)

if __name__ == "__main__":
    main()
