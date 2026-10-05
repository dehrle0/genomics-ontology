#!/usr/bin/env python3
"""
phase_duo_parental_orientation.py
=================================
Single-parent pedigree phase set orientation & isolated variant resolver.

Given:
1. Child WGS variant callset phased by WhatsHap using read pairs.
2. Single parent (maternal or paternal) microarray callset (e.g. 23andMe or AncestryDNA).

Function:
- Identifies parental homozygous markers (0/0 or 1/1).
- Assesses parental phase orientation across every WhatsHap phase set (PS block).
- Uniformly flips blocks to maintain standard orientation:
    - Haplotype 1 = Paternal allele
    - Haplotype 2 = Maternal allele
- Resolves and directly phases isolated heterozygous variants where the single parent is homozygous.
- Outputs standardized, parentally oriented VCF 4.2.

Author: Genomics Contributor <developer@users.noreply.github.com>
License: MIT / BSD-3
"""

import sys
import os
import argparse
import pysam

def phase_duo(child_vcf_path, parent_vcf_path, output_vcf_path, child_sample, parent_sample, parent_role="maternal", chrom=None):
    print(f"[duo_phase] Child VCF: {child_vcf_path} (Sample: {child_sample})")
    print(f"[duo_phase] Parent VCF: {parent_vcf_path} (Sample: {parent_sample}, Role: {parent_role})")
    print(f"[duo_phase] Target Contig: {chrom if chrom else 'ALL'}")

    parent_in = pysam.VariantFile(parent_vcf_path)
    parent_markers = {}  # (chrom, pos) -> (ref, alts, allele)

    fetch_args = [chrom] if chrom else []
    try:
        for rec in parent_in.fetch(*fetch_args):
            gt = rec.samples[parent_sample]['GT']
            if gt in [(0, 0), (1, 1)]:
                parent_markers[(rec.chrom, rec.pos)] = (rec.ref, rec.alts, gt[0])
    except KeyError:
        print(f"[duo_phase] Warning: sample {parent_sample} not found or contig missing in parent VCF.")
    parent_in.close()

    print(f"[duo_phase] Loaded {len(parent_markers)} informative homozygous parent markers.")

    # Pass 1: Tally orientation per phase set (PS)
    child_in = pysam.VariantFile(child_vcf_path)
    block_stats = {}  # (chrom, ps) -> {'h1_matches': 0, 'h2_matches': 0}

    for rec in child_in.fetch(*fetch_args):
        gt = rec.samples[child_sample]['GT']
        ps = rec.samples[child_sample].get('PS')

        if gt in [(0, 1), (1, 0)] and ps is not None:
            key = (rec.chrom, rec.pos)
            if key in parent_markers:
                p_ref, p_alts, p_allele = parent_markers[key]
                if p_ref == rec.ref and p_alts == rec.alts:
                    ps_key = (rec.chrom, ps)
                    if ps_key not in block_stats:
                        block_stats[ps_key] = {'h1_matches': 0, 'h2_matches': 0}

                    # gt[0] = H1, gt[1] = H2
                    if gt[0] == p_allele:
                        block_stats[ps_key]['h1_matches'] += 1
                    elif gt[1] == p_allele:
                        block_stats[ps_key]['h2_matches'] += 1

    child_in.close()

    # Determine which blocks to invert:
    # Standard: H1 = Paternal, H2 = Maternal
    # If parent is MATERNAL: we want p_allele to match H2.
    # If h1_matches > h2_matches, we invert the block.
    # If parent is PATERNAL: we want p_allele to match H1.
    # If h2_matches > h1_matches, we invert the block.
    blocks_to_invert = set()
    for ps_key, counts in block_stats.items():
        if parent_role.lower() == "maternal":
            if counts['h1_matches'] > counts['h2_matches']:
                blocks_to_invert.add(ps_key)
        else:  # paternal
            if counts['h2_matches'] > counts['h1_matches']:
                blocks_to_invert.add(ps_key)

    print(f"[duo_phase] Assessed {len(block_stats)} phase sets. Inverting {len(blocks_to_invert)} sets for canonical parental alignment.")

    # Pass 2: Write normalized, parentally oriented callset
    child_in = pysam.VariantFile(child_vcf_path)
    header = child_in.header.copy()

    if 'PS' not in header.formats:
        header.formats.add('PS', 1, 'Integer', 'Phase set identifier')

    out_vcf = pysam.VariantFile(output_vcf_path, 'w', header=header)

    phased_in_blocks = 0
    direct_rescued = 0

    for rec in child_in.fetch(*fetch_args):
        new_rec = rec.copy()
        gt = rec.samples[child_sample]['GT']
        ps = rec.samples[child_sample].get('PS')

        if gt in [(0, 1), (1, 0)]:
            ps_key = (rec.chrom, ps)
            if ps is not None:
                # Inside an existing read-backed phase set
                if ps_key in blocks_to_invert:
                    new_rec.samples[child_sample]['GT'] = (gt[1], gt[0])
                new_rec.samples[child_sample].phased = True
                phased_in_blocks += 1
            else:
                # Isolated heterozygous site: check single parent
                key = (rec.chrom, rec.pos)
                if key in parent_markers:
                    p_ref, p_alts, p_allele = parent_markers[key]
                    if p_ref == rec.ref and p_alts == rec.alts:
                        # Parent transmitted p_allele.
                        # Other allele must come from unsequenced parent.
                        alt_allele = 1 if p_allele == 0 else 0
                        if parent_role.lower() == "maternal":
                            # H1 = Paternal (alt_allele), H2 = Maternal (p_allele)
                            new_rec.samples[child_sample]['GT'] = (alt_allele, p_allele)
                        else:
                            # H1 = Paternal (p_allele), H2 = Maternal (alt_allele)
                            new_rec.samples[child_sample]['GT'] = (p_allele, alt_allele)

                        new_rec.samples[child_sample]['PS'] = rec.pos
                        new_rec.samples[child_sample].phased = True
                        direct_rescued += 1

        out_vcf.write(new_rec)

    child_in.close()
    out_vcf.close()

    print(f"[duo_phase] Complete: {phased_in_blocks} variants phased in oriented blocks, {direct_rescued} isolated sites parentally resolved.")


def main():
    parser = argparse.ArgumentParser(description="Orient phase blocks and phase isolated variants using a single parent microarray.")
    parser.add_argument("--child-vcf", required=True, help="Child VCF (WhatsHap read-phased)")
    parser.add_argument("--parent-vcf", required=True, help="Parent VCF (23andMe or AncestryDNA converted)")
    parser.add_argument("--output-vcf", "-o", required=True, help="Output parentally phased VCF")
    parser.add_argument("--child-sample", default="PROBAND", help="Child sample column name (default: PROBAND)")
    parser.add_argument("--parent-sample", default="PARENT", help="Parent sample column name (default: PARENT)")
    parser.add_argument("--parent-role", choices=["maternal", "paternal"], default="maternal", help="Role of sequenced single parent")
    parser.add_argument("--chrom", default=None, help="Target chromosome (e.g. chr20). Omit for whole genome.")

    args = parser.parse_args()
    phase_duo(args.child_vcf, args.parent_vcf, args.output_vcf, args.child_sample, args.parent_sample, args.parent_role, args.chrom)

if __name__ == "__main__":
    main()
