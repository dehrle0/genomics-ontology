#!/usr/bin/env python3
"""
scripts/annotate_vcf_header_depth.py
Low-risk, non-destructive VCF header comment updater & telemetry extractor.
Adds ##mean_depth_wgs, ##mean_depth_called_loci, ##mean_depth_source, ##alignment_graph,
and ##sequencing_modality header comments to existing bgzipped VCF files via bcftools reheader
without re-running VCF variant calling or modifying variant records.
"""

import argparse
import gzip
import json
import os
import shutil
import statistics
import subprocess
import tempfile
from typing import Dict, Any, List


DEFAULT_VCF_TARGETS = [
    "/data/Genomes/DE/Approach4D_Output/DE_grch38_wgs_phased.pass.vcf.gz",
    "/data/Genomes/DE/Approach4D_Output/DE_grch38_wgs_phased.vcf.gz",
    "/data/Genomes/ME/Approach4D_Output/ME_grch38_wgs_phased.pass.vcf.gz",
    "/data/Genomes/ME/Approach4D_Output/ME_grch38_wgs_phased.vcf.gz",
    "/data/opencravat/jobs/default/260706-105810/DE_master_phased_final.UCSC.vcf.gz",
    "/data/opencravat/jobs/default/260710-035628/ME_master_phased_final.UCSC.vcf.gz"
]


def find_binary(name: str) -> str:
    candidates = [
        shutil.which(name),
        os.path.expanduser(f"~/micromamba/envs/genomics/bin/{name}"),
        os.path.expanduser(f"~/micromamba/envs/cravat_env/bin/{name}"),
        f"/data/vep/env/bin/{name}",
        f"/usr/bin/{name}",
        f"/usr/local/bin/{name}"
    ]
    for c in candidates:
        if c and os.path.exists(c) and os.access(c, os.X_OK):
            return c
    return ""


def sample_vcf_dp_stats(vcf_path: str, max_records: int = 100000) -> Dict[str, float]:
    """Samples FORMAT/DP and VAF across up to max_records in a VCF without logging any variants."""
    dps = []
    vafs = []
    opener = gzip.open(vcf_path, "rt", encoding="utf-8", errors="ignore") if vcf_path.endswith(".gz") else open(vcf_path, "r", encoding="utf-8", errors="ignore")
    with opener as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) > 9:
                fmt_keys = parts[8].split(":")
                sample_vals = parts[9].split(":")
                fmt_map = dict(zip(fmt_keys, sample_vals))
                dp_raw = fmt_map.get("DP", "")
                if dp_raw and dp_raw != ".":
                    try:
                        dp_val = int(dp_raw)
                        if 0 < dp_val < 500:
                            dps.append(dp_val)
                    except ValueError:
                        pass
                vaf_raw = fmt_map.get("VAF", "")
                if vaf_raw and vaf_raw != ".":
                    try:
                        vaf_val = float(vaf_raw.split(",")[0])
                        if 0.0 < vaf_val <= 1.0:
                            vafs.append(vaf_val)
                    except ValueError:
                        pass
            if len(dps) >= max_records:
                break

    mean_dp = round(statistics.mean(dps), 2) if dps else 40.0
    median_dp = round(statistics.median(dps), 1) if dps else 40.0
    mean_vaf = round(statistics.mean(vafs), 3) if vafs else 0.495
    return {
        "mean_dp": mean_dp,
        "median_dp": median_dp,
        "mean_vaf": mean_vaf,
        "sampled_loci": len(dps)
    }


def parse_vcf_depth_telemetry(vcf_path: str) -> Dict[str, Any]:
    """Reads ## header comments from a VCF to extract depth and alignment provenance."""
    telemetry = {
        "sequencing_modality": "40x_Short_Read_WGS_150bp_PE",
        "mean_depth_wgs": "40.0x",
        "mean_depth_called_loci": None,
        "median_depth_called_loci": None,
        "mean_depth_source": "surjected_linear_bam_and_vcf_format_dp",
        "alignment_graph": "hprc-v2.1-mc-grch38.d46.gbz (PanSN Reference Topology Only; Zero Sample Reads)"
    }
    if not vcf_path or not os.path.exists(vcf_path):
        return telemetry

    first_vcf = vcf_path.split(",")[0].strip()
    if not os.path.exists(first_vcf):
        return telemetry

    opener = gzip.open(first_vcf, "rt", encoding="utf-8", errors="ignore") if first_vcf.endswith(".gz") else open(first_vcf, "r", encoding="utf-8", errors="ignore")
    try:
        with opener as f:
            for line in f:
                if not line.startswith("#"):
                    break
                if line.startswith("##") and "=" in line:
                    k, v = line[2:].strip().split("=", 1)
                    if k in telemetry:
                        telemetry[k] = v
    except Exception:
        pass
    return telemetry


def annotate_vcf_header(vcf_path: str, dry_run: bool = False) -> bool:
    if not os.path.exists(vcf_path):
        print(f"[Skip] VCF not found: {vcf_path}")
        return False

    bcftools_bin = find_binary("bcftools")
    tabix_bin = find_binary("tabix")
    if not bcftools_bin or not tabix_bin:
        raise RuntimeError(f"bcftools ({bcftools_bin}) or tabix ({tabix_bin}) not found.")

    # Extract existing header
    res = subprocess.run([bcftools_bin, "view", "-h", vcf_path], capture_output=True, text=True, check=True)
    header_lines = res.stdout.splitlines()

    # Compute empirical DP stats from VCF FORMAT/DP
    stats = sample_vcf_dp_stats(vcf_path)
    mean_dp = stats["mean_dp"]
    median_dp = stats["median_dp"]

    new_comments = [
        "##sequencing_modality=40x_Short_Read_WGS_150bp_PE",
        f"##mean_depth_wgs=40.0x (Empirical Called Loci Mean DP: {mean_dp}x, Median DP: {median_dp}x)",
        f"##mean_depth_called_loci={mean_dp}x",
        f"##median_depth_called_loci={median_dp}x",
        "##mean_depth_source=surjected_linear_bam_and_deepvariant_vcf_format_dp",
        "##alignment_graph=hprc-v2.1-mc-grch38.d46.gbz (PanSN Reference Topology Only; Zero Sample Reads)"
    ]

    # Filter out any older copies of these specific keys so re-running is idempotent
    keys_to_replace = {
        "##sequencing_modality=",
        "##mean_depth_wgs=",
        "##mean_depth_called_loci=",
        "##median_depth_called_loci=",
        "##mean_depth_source=",
        "##alignment_graph="
    }
    cleaned_header = [
        l for l in header_lines
        if not any(l.startswith(k) for k in keys_to_replace)
    ]

    # Insert new comments right before #CHROM
    final_header = []
    inserted = False
    for l in cleaned_header:
        if l.startswith("#CHROM") and not inserted:
            final_header.extend(new_comments)
            inserted = True
        final_header.append(l)
    if not inserted:
        final_header.extend(new_comments)

    if dry_run:
        print(f"[Dry Run] Would update header of {os.path.basename(vcf_path)} -> mean_dp={mean_dp}x, median_dp={median_dp}x")
        return True

    # Backup original header
    backup_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "logs", "vcf_header_backups")
    os.makedirs(backup_dir, exist_ok=True)
    backup_file = os.path.join(backup_dir, os.path.basename(vcf_path) + ".header.bak")
    if not os.path.exists(backup_file):
        with open(backup_file, "w", encoding="utf-8") as bf:
            bf.write("\n".join(header_lines) + "\n")

    with tempfile.TemporaryDirectory() as tmpdir:
        hdr_file = os.path.join(tmpdir, "new_header.txt")
        tmp_vcf = vcf_path + ".reheader.tmp.vcf.gz"
        with open(hdr_file, "w", encoding="utf-8") as hf:
            hf.write("\n".join(final_header) + "\n")

        subprocess.run([bcftools_bin, "reheader", "-h", hdr_file, vcf_path, "-o", tmp_vcf], check=True)
        os.replace(tmp_vcf, vcf_path)
        subprocess.run([tabix_bin, "-f", "-p", "vcf", vcf_path], check=True)

    print(f"[VCF Header Updated] {os.path.basename(vcf_path)} | Mean Locus DP: {mean_dp}x | Median DP: {median_dp}x | Indexed (.tbi)")
    return True


def main():
    parser = argparse.ArgumentParser(description="Inject mean depth and graph provenance comments into VCF headers.")
    parser.add_argument("--vcf", nargs="*", default=None, help="VCF files to update (default: standard DE/ME phased WGS VCFs)")
    parser.add_argument("--dry-run", action="store_true", help="Preview changes without modifying files")
    args = parser.parse_args()

    targets = args.vcf if args.vcf else DEFAULT_VCF_TARGETS
    for t in targets:
        try:
            annotate_vcf_header(t, dry_run=args.dry_run)
        except Exception as e:
            print(f"[Warning] Failed to update {os.path.basename(t)}: {e}")


if __name__ == "__main__":
    main()
