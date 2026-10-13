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

AVI_SCRIPT = os.path.expanduser("~/.gemini/config/plugins/science/skills/alphagenome_variant_impact_score/scripts/alphagenome_atlas_avi.py")
UV_BIN = os.path.expanduser("~/.local/bin/uv")
GLOBAL_CACHE_FILE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "alphagenome_cache.json")


def parse_args():
    parser = argparse.ArgumentParser(description="Enrich actionable variants with AlphaGenome AVI predictions.")
    parser.add_argument("--in-json", required=True, help="Path to master_actionable.json")
    parser.add_argument("--cache", default=GLOBAL_CACHE_FILE, help="Path to global cache JSON")
    parser.add_argument("--sample-cache", default=None, help="Optional sample-specific cache JSON copy")
    parser.add_argument("--max-queries", type=int, default=5000, help="Max API queries allowed in a single run")
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


OFFLINE_TABIX_TABLE = "/data/vep/data/alphagenome/avi_scores.tsv.gz"


def is_snv(ref: str, alt: str) -> bool:
    r = str(ref or "").strip().upper()
    a = str(alt or "").strip().upper()
    return len(r) == 1 and len(a) == 1 and r in "ACGT" and a in "ACGT"


def query_offline_tabix(variants: List[str], table_path: str = OFFLINE_TABIX_TABLE) -> Dict[str, Dict[str, Any]]:
    if not variants or not os.path.exists(table_path):
        return {}

    tabix_bin = shutil.which("tabix") or "/data/vep/env/bin/tabix"
    if not (shutil.which(tabix_bin) or os.path.exists(tabix_bin)):
        return {}

    results = {}
    chunk_size = 500
    for i in range(0, len(variants), chunk_size):
        chunk = variants[i:i + chunk_size]
        regions = []
        lookup = {}
        for v in chunk:
            parts = v.split(":")
            if len(parts) != 3:
                continue
            c = parts[0]
            try:
                p = int(parts[1])
            except ValueError:
                continue
            alts = parts[2].split(">")
            if len(alts) != 2:
                continue
            ref, alt = alts[0], alts[1]
            regions.append(f"{c}:{p}-{p}")
            lookup[(c, p, ref, alt)] = v

        if not regions:
            continue

        cmd = [tabix_bin, table_path] + regions
        try:
            res = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
            if res.returncode == 0:
                for line in res.stdout.splitlines():
                    fields = line.strip().split("\t")
                    if len(fields) >= 6:
                        c, p, r, a = fields[0], int(fields[1]), fields[2], fields[3]
                        raw = float(fields[4]) if fields[4] != "." else None
                        phred = float(fields[5]) if fields[5] != "." else None
                        if (c, p, r, a) in lookup:
                            vk = lookup[(c, p, r, a)]
                            results[vk] = {
                                "supported": True,
                                "avi_phred": phred,
                                "avi_raw": raw,
                                "source": "offline_tabix"
                            }
        except Exception as e:
            print(f"[AlphaGenome Enrich Warning] Tabix offline lookup error: {e}")

    return results


def query_alphagenome_batch(variants: List[str], chunk_size: int = 100, cache_ref: Dict[str, Any] = None, cache_path: str = None) -> List[Dict[str, Any]]:
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
                    chunk_items = []
                    if isinstance(chunk_data, list):
                        chunk_items = chunk_data
                    elif isinstance(chunk_data, dict):
                        chunk_items = [chunk_data]
                    all_results.extend(chunk_items)

                    # Incrementally populate and persist cache per completed chunk
                    if cache_ref is not None:
                        for item in chunk_items:
                            vk = item.get("variant")
                            if not vk:
                                continue
                            features = {k: v for k, v in item.items() if k.startswith("fi_")}
                            cache_ref[vk] = {
                                "supported": True,
                                "avi_phred": item.get("avi_phred"),
                                "avi_raw": item.get("avi_raw"),
                                "avi_quantile": item.get("avi_quantile"),
                                "top_percentile": item.get("top_percentile"),
                                "top_modality": item.get("top_modality"),
                                "top_feature_importance": item.get("top_feature_importance"),
                                "features": features
                            }
                        if cache_path:
                            save_cache(cache_ref, cache_path)
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

    # Evaluate 100% of actionable records across all tiers and categories
    candidates = records

    print(f"[AlphaGenome Enrich] Total Actionable Records : {len(records)}")
    print(f"[AlphaGenome Enrich] Candidates for Evaluation: {len(candidates)}")
    print(f"[AlphaGenome Enrich] Cached Entries in Memory  : {len(cache)}")

    to_query = []
    to_query_set = set()
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

        if var_key not in to_query_set:
            to_query_set.add(var_key)
            if len(to_query) < args.max_queries:
                to_query.append(var_key)

    print(f"[AlphaGenome Enrich] Existing Cache Hits       : {cached_hits}")
    print(f"[AlphaGenome Enrich] Indels/Non-SNV Skipped    : {skipped_indels}")
    print(f"[AlphaGenome Enrich] Unique Variants to Query  : {len(to_query)}")

    if to_query and os.path.exists(OFFLINE_TABIX_TABLE):
        print(f"[AlphaGenome Enrich] Querying offline Tabix database for {len(to_query)} variants...")
        tabix_hits = query_offline_tabix(to_query, OFFLINE_TABIX_TABLE)
        print(f"[AlphaGenome Enrich] Offline Tabix Hits: {len(tabix_hits)}/{len(to_query)}")
        for vk, hit in tabix_hits.items():
            cache[vk] = hit
        save_cache(cache, args.cache)
        to_query = [v for v in to_query if v not in cache]
        print(f"[AlphaGenome Enrich] Remaining to query via API: {len(to_query)}")

    if to_query:
        batch_results = query_alphagenome_batch(to_query, cache_ref=cache, cache_path=args.cache)
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
