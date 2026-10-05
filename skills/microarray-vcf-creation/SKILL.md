---
name: microarray-vcf-creation
description: >-
  Parse, normalize, and convert commercial microarray raw data (AncestryDNA 5-column TSV and
  23andMe 4-column TSV) into standards-compliant, sorted, bgzipped, and tabix-indexed VCF 4.2.
  Handles reference genome coordinate lookups (GRCh38 / hg19 / T2T-CHM13), sex chromosome hemizygosity,
  and CrossMap / Picard liftover transformations.
---

# Microarray to VCF Creation Skill (AncestryDNA & 23andMe)

This skill provides an authoritative, self-contained reference and execution engine for transforming raw consumer microarray text dumps (AncestryDNA and 23andMe) into standardized, indexed VCF 4.2 files suitable for downstream variant calling, imputation, and pedigree phasing.

---

## 1. Input Formats & Idiosyncrasies

Commercial genotyping arrays produce flat TSV/CSV tables without reference alleles or standard VCF headers:

| Provider | Column Layout | Typical Coordinates | Hemizygous Calling | Missing Values |
| :--- | :--- | :--- | :--- | :--- |
| **AncestryDNA** | `rsid`, `chromosome`, `position`, `allele1`, `allele2` (5 columns) | GRCh37 / hg19 | Split into two columns (`allele1`, `allele2`) | `0 0`, `- -`, `? ?` |
| **23andMe** | `rsid`, `chromosome`, `position`, `genotype` (4 columns) | GRCh37 / hg19 (v3/v4/v5) | Concatenated single string (e.g. `AG`, `AA`, `A`) | `--`, `DD`, `II`, `??` |

### Critical Normalization Rules:
1. **Contig Renaming:** Microarrays typically use numeric strings for chromosomes. Normalize:
   - `23` $\rightarrow$ `chrX`
   - `24` $\rightarrow$ `chrY`
   - `25` $\rightarrow$ `chrX` (PAR: Pseudoautosomal Region)
   - `26` $\rightarrow$ `chrM`
2. **Reference Allele Retrieval:** Microarrays report only observed patient alleles, NOT the reference allele. A reference FASTA or `.2bit` file (hg19 or GRCh38) is mandatory to look up `REF` at `POS` (1-based index).
3. **Ploidy and Hemizygosity:**
   - Autosomes (`chr1`–`chr22`): Always diploid. If a single base is reported for homozygous calls, expand to homozygous pair (e.g., `A` $\rightarrow$ `A/A`).
   - Female `chrX`: Diploid (`0/0`, `0/1`, `1/1`).
   - Male `chrX` (non-PAR) and `chrY`: Hemizygous (`0` or `1`).
   - Mitochondrial (`chrM`): Hemizygous (`0` or `1`).
4. **Indel Markers (`D` / `I`):** Commercial arrays denote micro-indels as `D` (deletion) or `I` (insertion). Because these lack explicit nucleotide sequences, record them as missing (`./.`) or filter to strictly preserve high-confidence SNV anchors.

---

## 2. Embedded Production Converter: `convert_microarray_to_vcf.py`

Execute the following standalone Python script to convert raw AncestryDNA or 23andMe TSV data into VCF 4.2 format:

```python
#!/usr/bin/env python3
import sys
import os
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

def convert_microarray(input_path, output_vcf, ref_path, sample_name, sex="female"):
    ref = ReferenceEngine(ref_path)
    is_ancestry = False
    variants = []

    with open(input_path, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            lower = line.lower()
            if 'rsid' in lower or 'chromosome' in lower:
                if 'allele1' in lower:
                    is_ancestry = True
                continue

            parts = line.strip().split()
            if len(parts) < 4:
                continue

            rsid, chrom_raw, pos_raw = parts[0], parts[1].replace("chr", ""), parts[2]
            try:
                pos = int(pos_raw)
            except ValueError:
                continue

            chrom_std = CHR_MAP.get(chrom_raw, f"chr{chrom_raw}")
            if is_ancestry or len(parts) == 5:
                a1, a2 = parts[3].upper(), parts[4].upper()
            else:
                gt = parts[3].upper()
                a1, a2 = (gt[0], gt[1]) if len(gt) == 2 else (gt, None) if len(gt) == 1 else ('-', '-')

            is_autosome = chrom_std in [f"chr{i}" for i in range(1, 23)]
            is_female_x = (chrom_std == "chrX" and sex.lower() in ("female", "f", "xx"))
            if (is_autosome or is_female_x) and a1 not in ('-', '0', '?', 'N', 'D', 'I') and a2 is None:
                a2 = a1

            if a1 in ('-', '0', '?', 'N', 'D', 'I') or (a2 is not None and a2 in ('-', '0', '?', 'N', 'D', 'I')):
                variants.append((chrom_std, pos, rsid, "N", ".", "./."))
                continue

            ref_base = ref.get_base(chrom_std, pos)
            if not ref_base or ref_base not in ('A', 'C', 'G', 'T'):
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

    variants.sort(key=lambda x: (get_chrom_sort_key(x[0]), x[1]))

    with open(output_vcf, 'w') as out:
        out.write("##fileformat=VCFv4.2\n##source=convert_microarray_to_vcf.py\n")
        out.write(f"##reference={ref_path}\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        for c in sorted(CHROM_ORDER.keys(), key=lambda k: CHROM_ORDER[k]):
            c_name = c if ref.has_chr_prefix else c.replace("chr", "")
            length = ref.get_contig_length(c_name)
            if length:
                out.write(f"##contig=<ID={c_name},length={length}>\n")
        out.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n")
        for chrom_std, pos, rsid, ref_base, alt, gt in variants:
            c_name = chrom_std if ref.has_chr_prefix else chrom_std.replace("chr", "")
            out.write(f"{c_name}\t{pos}\t{rsid}\t{ref_base}\t{alt}\t.\tPASS\t.\tGT\t{gt}\n")

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--input", "-i", required=True)
    p.add_argument("--output", "-o", required=True)
    p.add_argument("--reference", "-r", required=True)
    p.add_argument("--sample", "-s", required=True)
    p.add_argument("--sex", default="female", choices=["female", "male"])
    args = p.parse_args()
    convert_microarray(args.input, args.output, args.reference, args.sample, args.sex)
```

---

## 3. Coordinate Liftover Workflow (hg19 $\rightarrow$ GRCh38)

If raw microarray data is in hg19/b37 coordinates, lift over the generated VCF to GRCh38 using CrossMap or Picard LiftoverVcf:

### Option A: CrossMap Liftover (Fast & Python-Native)
```bash
# 1. Compress and index initial hg19 VCF
bgzip -c sample_hg19.vcf > sample_hg19.vcf.gz
tabix -p vcf sample_hg19.vcf.gz

# 2. Lift over to GRCh38 using UCSC chain file
CrossMap vcf \
  hg19ToHg38.over.chain.gz \
  sample_hg19.vcf.gz \
  GRCh38_reference.fa \
  sample_hg38.vcf

# 3. Sort, compress, and index
bcftools sort -Oz -o sample_hg38.vcf.gz sample_hg38.vcf
tabix -p vcf sample_hg38.vcf.gz
```

### Option B: Picard LiftoverVcf (GATK / Java)
```bash
java -jar picard.jar LiftoverVcf \
  I=sample_hg19.vcf.gz \
  O=sample_hg38.vcf.gz \
  CHAIN=hg19ToHg38.over.chain.gz \
  REJECT=rejected_variants.vcf.gz \
  R=GRCh38_reference.fa \
  WARN_ON_MISSING_CONTIG=true
```

---

## 4. Verification & Quality Gates

Run `bcftools stats` on the final converted VCF to ensure integrity:
```bash
bcftools stats sample_hg38.vcf.gz | grep "^SN"
```
Expectations:
- **Total Records:** 500,000 to 720,000 loci (typical Illumina OmniExpress / Global Screening Array content).
- **SNPs:** >99.5% of total records.
- **Ti/Tv Ratio:** Approximately 2.0 to 2.3 for human autosomes.
