#!/usr/bin/env python3
"""
lib/enrich_alphagenome.py
Enriches actionable variants with Google DeepMind AlphaGenome Variant Impact (AVI) scores.
Enforces offline-first caching via data/alphagenome_cache.json and strict query capping.
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
from typing import Dict, Any, List

AVI_SCRIPT = "/home/daniel-ehrle/.gemini/config/plugins/science/skills/alphagenome_variant_impact_score/scripts/alphagenome_atlas_avi.py"
UV_BIN = os.path.expanduser("~/.local/bin/uv")
GLOBAL_CACHE_FILE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "alphagenome_cache.json")


def parse_args():
    parser = argparse.ArgumentParser(description="Enrich actionable variants with AlphaGenome AVI predictions.")
    parser.add_argument("--in-json", required=True, help="Path to master_actionable.json")
    parser.add_argument("--cache", default=GLOBAL_CACHE_FILE, help="Path to global cache JSON")
    parser.add_argument("--sample-cache", default=None, help="Optional sample-specific cache JSON copy")
    parser.add_argument("--max-queries", type=int, default=1000, help="Max API queries allowed in a single run")
    parser.add_argument("--force", action="store_true", help="Force re-query of cached items")
    return parser.parse_args()


def load_cache(cache_path: str) -> Dict[str, Any]:
    if os.path.exists(cache_path):
        try:
            with open(cache_path, "r", encoding="utf-8") as f:
                return json.load(f)
        except Exception as e:
            print(f"[AlphaGenome Enrich Warning] Failed to load cache from {cache_path}: {e}")
    return {}


def save_cache(cache_data: Dict[str, Any], cache_path: str):
    os.makedirs(os.path.dirname(os.path.abspath(cache_path)), exist_ok=True)
    temp_path = cache_path + ".tmp"
    with open(temp_path, "w", encoding="utf-8") as f:
        json.dump(cache_data, f, indent=2)
    os.replace(temp_path, cache_path)


def format_var_key(chrom: str, pos: Any, ref: str, alt: str) -> str:
    c = str(chrom)
    if not c.startswith("chr"):
        c = f"chr{c}"
    return f"{c}:{pos}:{ref}>{alt}"


def is_snv(ref: str, alt: str) -> bool:
    r = str(ref or "").strip().upper()
    a = str(alt or "").strip().upper()
    return len(r) == 1 and len(a) == 1 and r in "ACGT" and a in "ACGT"


def query_alphagenome_batch(variants: List[str], chunk_size: int = 100) -> List[Dict[str, Any]]:
    if not variants:
        return []

    uv_path = shutil.which("uv") or (UV_BIN if os.path.exists(UV_BIN) else None)
    if not uv_path:
        print("[AlphaGenome Enrich Error] uv binary not found. Cannot invoke AlphaGenome CLI.")
        return []

    if not os.path.exists(AVI_SCRIPT):
        print(f"[AlphaGenome Enrich Error] AVI script not found at: {AVI_SCRIPT}")
        return []

    all_results = []
    total_variants = len(variants)

    for start_idx in range(0, total_variants, chunk_size):
        chunk = variants[start_idx:start_idx + chunk_size]
        chunk_num = (start_idx // chunk_size) + 1
        total_chunks = (total_variants + chunk_size - 1) // chunk_size

        with tempfile.TemporaryDirectory() as tmpdir:
            input_file = os.path.join(tmpdir, "variants_to_query.txt")
            output_file = os.path.join(tmpdir, "avi_results.json")

            with open(input_file, "w", encoding="utf-8") as f:
                for v in chunk:
                    f.write(f"{v}\n")

            cmd = [
                uv_path, "run",
                "--with", "alphagenome",
                "--with", "python-dotenv",
                "--with", "anndata",
                "--with", "pandas",
                "--with", "polars",
                "--with", "tabulate",
                "python3", AVI_SCRIPT,
                "query",
                "-i", input_file,
                "--format", "json",
                "-o", output_file
            ]

            print(f"[AlphaGenome Enrich] Executing AlphaGenome Atlas query (Chunk {chunk_num}/{total_chunks}: {len(chunk)} variants) via uv...")
            res = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
            if res.returncode != 0:
                print(f"[AlphaGenome Enrich Warning] CLI returned non-zero code {res.returncode}:\n{res.stderr}")
                continue

            if not os.path.exists(output_file):
                print(f"[AlphaGenome Enrich Warning] Output file not found: {output_file}")
                continue

            try:
                with open(output_file, "r", encoding="utf-8") as f:
                    chunk_data = json.load(f)
                    if isinstance(chunk_data, list):
                        all_results.extend(chunk_data)
                    elif isinstance(chunk_data, dict):
                        all_results.append(chunk_data)
            except Exception as e:
                print(f"[AlphaGenome Enrich Error] Failed parsing JSON output: {e}")

    return all_results


def main():
    args = parse_args()

    if not os.path.exists(args.in_json):
        print(f"[AlphaGenome Enrich Error] File not found: {args.in_json}", file=sys.stderr)
        sys.exit(1)

    with open(args.in_json, "r", encoding="utf-8") as f:
        data = json.load(f)

    records = data.get("records", [])
    cache = load_cache(args.cache)

    candidates = []
    for r in records:
        tier = r.get("tier")
        rcodes = r.get("reason_codes") or []
        ev = r.get("evidence", {}) or {}
        is_t1_t2 = str(tier) in ("Tier1", "Tier2")
        is_cand = ("RESCUE_ALPHAGENOME_TARGET" in rcodes) or ev.get("is_alphagenome_candidate", False) or is_t1_t2
        if is_cand:
            candidates.append(r)

    print(f"[AlphaGenome Enrich] Total Actionable Records : {len(records)}")
    print(f"[AlphaGenome Enrich] Tier1/Tier2/Triage Targets : {len(candidates)}")
    print(f"[AlphaGenome Enrich] Cached Entries in Memory  : {len(cache)}")

    to_query = []
    skipped_indels = 0
    cached_hits = 0

    for r in candidates:
        chrom = r.get("chrom")
        pos = r.get("pos")
        ref = r.get("ref")
        alt = r.get("alt")
        var_key = format_var_key(chrom, pos, ref, alt)

        if not args.force and var_key in cache:
            cached_hits += 1
            continue

        if not is_snv(ref, alt):
            skipped_indels += 1
            cache[var_key] = {
                "supported": False,
                "reason": "INDEL_NOT_SUPPORTED",
                "avi_phred": None,
                "avi_raw": None,
                "top_percentile": None,
                "top_modality": None
            }
            continue

        if len(to_query) < args.max_queries:
            to_query.append(var_key)

    print(f"[AlphaGenome Enrich] Existing Cache Hits       : {cached_hits}")
    print(f"[AlphaGenome Enrich] Indels/Non-SNV Skipped    : {skipped_indels}")
    print(f"[AlphaGenome Enrich] Variants to Query (Live)  : {len(to_query)}")

    if to_query:
        batch_results = query_alphagenome_batch(to_query)
        for res in batch_results:
            v_key = res.get("variant")
            if not v_key:
                continue
            features = {k: v for k, v in res.items() if k.startswith("fi_")}
            cache[v_key] = {
                "supported": True,
                "avi_phred": res.get("avi_phred"),
                "avi_raw": res.get("avi_raw"),
                "avi_quantile": res.get("avi_quantile"),
                "top_percentile": res.get("top_percentile"),
                "top_modality": res.get("top_modality"),
                "top_feature_importance": res.get("top_feature_importance"),
                "features": features
            }
        print(f"[AlphaGenome Enrich] Ingested {len(batch_results)} live predictions into cache.")

    # Apply cache to actionable records
    annotated_count = 0
    for r in records:
        chrom = r.get("chrom")
        pos = r.get("pos")
        ref = r.get("ref")
        alt = r.get("alt")
        var_key = format_var_key(chrom, pos, ref, alt)

        if var_key in cache:
            entry = cache[var_key]
            ev = r.setdefault("evidence", {})
            ev["alphagenome_data"] = entry
            ev["avi_phred"] = entry.get("avi_phred")
            ev["avi_percentile"] = entry.get("top_percentile")
            ev["avi_modality"] = entry.get("top_modality")
            ev["avi_raw"] = entry.get("avi_raw")
            ev["avi_status"] = entry.get("reason") or ("SCORED" if entry.get("avi_phred") is not None else "UNSCORED")
            annotated_count += 1
        elif not is_snv(ref, alt):
            ev = r.setdefault("evidence", {})
            ev["avi_status"] = "INDEL_NOT_SUPPORTED"

    # Save actionable JSON
    with open(args.in_json, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)

    # Save cache
    save_cache(cache, args.cache)
    if args.sample_cache:
        save_cache(cache, args.sample_cache)

    print(f"[AlphaGenome Enrich] Total variants annotated in JSON: {annotated_count}")
    print(f"[AlphaGenome Enrich] Master cache saved to: {args.cache}")


if __name__ == "__main__":
    main()
