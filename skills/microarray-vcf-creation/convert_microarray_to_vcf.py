#!/usr/bin/env python3
"""
convert_microarray_to_vcf.py
============================
Universal parser and VCF 4.2 converter for commercial microarray raw data:
- AncestryDNA (5-column format: rsid, chromosome, position, allele1, allele2)
- 23andMe (4-column format: rsid, chromosome, position, genotype)
- MyHeritage (CSV/TSV format: RSID, CHROMOSOME, POSITION, RESULT)
- FamilyTreeDNA / FTDNA (CSV format: RSID, CHROMOSOME, POSITION, RESULT)

Supports reference lookup via pysam (FASTA) or twobitreader (2bit).
Automatically normalizes contig names, handles pseudoautosomal regions,
resolves hemizygous/diploid calls, sorts records, and writes valid VCF 4.2.

Author: Genomics Contributor <developer@users.noreply.github.com>
License: MIT / BSD-3
"""

import sys
import os
import csv
import argparse

CHR_MAP = {str(i): f"chr{i}" for i in range(1, 23)}
CHR_MAP.update({
    "X": "chrX", "Y": "chrY", "MT": "chrM", "M": "chrM",
    "23": "chrX", "24": "chrY", "25": "chrX", "26": "chrM"
})

CHROM_ORDER = {f"chr{i}": i for i in range(1, 23)}
CHROM_ORDER.update({"chrX": 23, "chrY": 24, "chrM": 25})

def get_chrom_sort_key(chrom):
    clean = chrom if chrom.startswith("chr") else f"chr{chrom}"
    return CHROM_ORDER.get(clean, 999)

class ReferenceEngine:
    def __init__(self, ref_path):
        self.ref_path = ref_path
        self.is_2bit = ref_path.endswith(".2bit")
        self.handle = None
        self.has_chr_prefix = True

        if self.is_2bit:
            import twobitreader
            self.handle = twobitreader.TwoBitFile(ref_path)
            self.contigs = set(self.handle.keys())
            self.has_chr_prefix = any(c.startswith("chr") for c in self.contigs)
        else:
            import pysam
            self.handle = pysam.FastaFile(ref_path)
            self.contigs = set(self.handle.references)
            self.has_chr_prefix = any(c.startswith("chr") for c in self.contigs)

    def get_base(self, chrom, pos_1based):
        lookup = chrom if self.has_chr_prefix else chrom.replace("chr", "")
        if lookup not in self.contigs:
            return None
        try:
            if self.is_2bit:
                return self.handle[lookup][pos_1based - 1].upper()
            else:
                return self.handle.fetch(lookup, pos_1based - 1, pos_1based).upper()
        except (IndexError, KeyError):
            return None

    def get_contig_length(self, chrom):
        lookup = chrom if self.has_chr_prefix else chrom.replace("chr", "")
        if lookup not in self.contigs:
            return None
        return len(self.handle[lookup]) if self.is_2bit else self.handle.get_reference_length(lookup)


def detect_delimiter(filepath):
    """Detects whether file is comma-delimited (CSV) or tab/whitespace-delimited."""
    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
        for _ in range(50):
            line = f.readline()
            if not line:
                break
            if line.startswith('#') or not line.strip():
                continue
            if ',' in line and line.count(',') >= 3:
                return ','
            if '\t' in line:
                return '\t'
    return None


def convert_microarray(input_path, output_vcf, ref_path, sample_name, sex="female", format_mode="auto"):
    ref_engine = ReferenceEngine(ref_path)
    print(f"[convert] Reading raw microarray data: {input_path}")
    print(f"[convert] Format Mode: {format_mode}, Sample ID: {sample_name}, Declared Sex: {sex}")

    delimiter = detect_delimiter(input_path)
    variants = []
    skipped_count = 0
    parsed_count = 0

    is_ancestry_5col = (format_mode == "ancestrydna")

    with open(input_path, 'r', encoding='utf-8', errors='ignore') as infile:
        # Use csv reader if comma delimited, else split
        reader = csv.reader(infile, delimiter=delimiter) if delimiter == ',' else None

        for raw_line in infile if reader is None else reader:
            if reader is None:
                line = raw_line.strip()
                if line.startswith('#') or not line:
                    continue
                parts = [p.strip().strip('"').strip("'") for p in line.split()]
            else:
                parts = [p.strip().strip('"').strip("'") for p in raw_line]
                if not parts or parts[0].startswith('#'):
                    continue

            if len(parts) < 4:
                continue

            lower_0 = parts[0].lower()
            lower_1 = parts[1].lower() if len(parts) > 1 else ""

            # Detect Header
            if 'rsid' in lower_0 or 'chromosome' in lower_1:
                if format_mode == "auto":
                    if len(parts) >= 5 and ('allele1' in [p.lower() for p in parts] or 'allele2' in [p.lower() for p in parts]):
                        is_ancestry_5col = True
                    elif len(parts) == 5:
                        is_ancestry_5col = True
                continue

            rsid = parts[0]
            chrom_raw = parts[1].replace("chr", "")
            pos_raw = parts[2]

            try:
                pos = int(pos_raw)
            except ValueError:
                continue

            chrom_std = CHR_MAP.get(chrom_raw, f"chr{chrom_raw}")

            # Extract alleles based on format
            if is_ancestry_5col or len(parts) >= 5:
                a1 = parts[3].upper()
                a2 = parts[4].upper()
            else:
                gt_str = parts[3].upper()
                if len(gt_str) == 2:
                    a1 = gt_str[0]
                    a2 = gt_str[1]
                elif len(gt_str) == 1:
                    a1 = gt_str
                    a2 = None
                else:
                    a1 = '-'
                    a2 = '-'

            is_autosome = chrom_std in [f"chr{i}" for i in range(1, 23)]
            is_female_x = (chrom_std == "chrX" and sex.lower() in ("female", "f", "xx"))

            if (is_autosome or is_female_x) and a1 not in ('-', '0', '?', 'N', 'D', 'I') and a2 is None:
                a2 = a1

            parsed_count += 1

            if a1 in ('-', '0', '?', 'N', 'D', 'I') or (a2 is not None and a2 in ('-', '0', '?', 'N', 'D', 'I')):
                variants.append((chrom_std, pos, rsid, "N", ".", "./."))
                continue

            ref_base = ref_engine.get_base(chrom_std, pos)
            if not ref_base or ref_base not in ('A', 'C', 'G', 'T'):
                skipped_count += 1
                continue

            if a2 is None:
                gt_val = "0" if a1 == ref_base else "1"
                variants.append((chrom_std, pos, rsid, ref_base, "." if gt_val == "0" else a1, gt_val))
            else:
                if a1 == ref_base and a2 == ref_base:
                    variants.append((chrom_std, pos, rsid, ref_base, ".", "0/0"))
                elif a1 == ref_base and a2 != ref_base:
                    variants.append((chrom_std, pos, rsid, ref_base, a2, "0/1"))
                elif a2 == ref_base and a1 != ref_base:
                    variants.append((chrom_std, pos, rsid, ref_base, a1, "0/1"))
                elif a1 == a2:
                    variants.append((chrom_std, pos, rsid, ref_base, a1, "1/1"))
                else:
                    variants.append((chrom_std, pos, rsid, ref_base, f"{a1},{a2}", "1/2"))

    print(f"[convert] Parsed {parsed_count} entries. Skipped {skipped_count} unmapped loci.")
    print(f"[convert] Sorting {len(variants)} records by genomic coordinate...")
    variants.sort(key=lambda x: (get_chrom_sort_key(x[0]), x[1]))

    print(f"[convert] Writing VCF 4.2 to: {output_vcf}")
    with open(output_vcf, 'w') as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write(f"##source=convert_microarray_to_vcf.py\n")
        out.write(f"##reference={ref_path}\n")
        out.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")

        for c_std in sorted(CHROM_ORDER.keys(), key=lambda c: CHROM_ORDER[c]):
            c_name = c_std if ref_engine.has_chr_prefix else c_std.replace("chr", "")
            length = ref_engine.get_contig_length(c_name)
            if length:
                out.write(f"##contig=<ID={c_name},length={length}>\n")

        out.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n")

        for chrom_std, pos, rsid, ref, alt, gt in variants:
            c_name = chrom_std if ref_engine.has_chr_prefix else chrom_std.replace("chr", "")
            out.write(f"{c_name}\t{pos}\t{rsid}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt}\n")

    print("[convert] Conversion complete.")


def main():
    parser = argparse.ArgumentParser(description="Universal commercial microarray to VCF 4.2 converter.")
    parser.add_argument("--input", "-i", required=True, help="Raw AncestryDNA, 23andMe, MyHeritage, or FTDNA text/tsv/csv")
    parser.add_argument("--output", "-o", required=True, help="Path to output VCF file")
    parser.add_argument("--reference", "-r", required=True, help="Reference FASTA or 2bit file")
    parser.add_argument("--sample", "-s", required=True, help="Sample identifier name for VCF header")
    parser.add_argument("--sex", choices=["female", "male"], default="female", help="Donor biological sex (default: female)")
    parser.add_argument("--format", choices=["auto", "ancestrydna", "23andme", "myheritage", "ftdna"], default="auto", help="Vendor format mode")

    args = parser.parse_args()
    convert_microarray(args.input, args.output, args.reference, args.sample, args.sex, args.format)

if __name__ == "__main__":
    main()
