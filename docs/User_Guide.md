# User Guide: Genomic Ontology Reporting & Clinical Pipeline (v5.1)

Welcome to the **Genomic Ontology Reporting & Clinical Pipeline**. This platform translates Whole Genome Sequencing (WGS), Whole Exome Sequencing (WES), and targeted panel variant calls into formal multi-level Directed Acyclic Graph (DAG) biomedical ontologies (**HPO**, **GO**, and **Anatomical Organ/System** views), integrating in-silico predictors, polygenic risk scores (PRS), pharmacogenomic recommendations (CPIC/DPWG), and peer-reviewed literature.

---

## 1. Environment & Prerequisites

Activate your micromamba or conda environment containing Python 3.10+ and OpenCRAVAT:

```bash
micromamba activate cravat_env
```

### 1.1 Python Dependencies
Install required packages:
```bash
pip install pyyaml openpyxl pydantic linkml-runtime requests
```

---

## 2. Running OpenCRAVAT (OC) Annotation

To generate the rich multi-annotator SQLite database from raw VCF files, run OpenCRAVAT with the full suite of clinical and functional annotators.

### 2.1 Complete OpenCRAVAT Annotators List
The pipeline leverages the following **25 essential annotators**:

* **Clinical & Pathogenicity**: `clinvar`, `clingen`, `civic`, `clinvar_acmg`
* **Ontologies**: `hpo`, `go`
* **Gene Models & Disease Databases**: `ncbigene`, `omim`, `interpro`
* **In-Silico Missense & Splicing Predictors**: `revel`, `alphamissense`, `cadd`, `spliceai`, `bayesdel`, `metarnn`, `esm1b`, `varity`
* **Non-Coding & Regulatory**: `linsight`, `ncer`, `regulomedb`, `ccre_screen`
* **Population & Frequency**: `gnomad4`, `dbsnp`, `vcfinfo`, `gtex`
* **Pharmacogenomics & GWAS**: `pharmgkb`, `gwas_catalog`

### 2.2 Bash Command to Install All Annotators
```bash
oc module install hpo go clinvar clingen omim ncbigene revel alphamissense \
  bayesdel metarnn esm1b varity spliceai cadd linsight ncer regulomedb \
  ccre_screen gtex dbsnp vcfinfo gwas_catalog pharmgkb civic interpro
```

### 2.3 Bash Command to Run OpenCRAVAT on a VCF
```bash
oc run /path/to/input_phased.vcf.gz \
  -l hg38 \
  -a hpo go clinvar clingen omim ncbigene revel alphamissense bayesdel metarnn esm1b varity spliceai cadd linsight ncer regulomedb ccre_screen gtex dbsnp vcfinfo gwas_catalog pharmgkb civic interpro \
  -d /path/to/output_dir \
  --mp $(nproc) \
  -n Sample_Prefix
```

---

## 3. Running the Unified Ontology Pipeline Engine

You can execute the entire end-to-end report generation using either the Python runner (`run_ontology_pipeline.py`) or the Bash wrapper (`run_ontology_master_pipeline.sh`).

### 3.1 Accepted Input Sources
1. **Raw or Phased VCF (`.vcf`, `.vcf.gz`)**: Automatically runs `oc run` with all 25 annotators if an annotated database does not already exist.
2. **Pre-Annotated SQLite Database (`.sqlite`)**: Immediately executes filtering, enrichment, and report rendering.
3. **OpenCRAVAT Job ID (e.g. `260706-105810`)**: Automatically locates the SQLite database and phased VCF in `/data/opencravat/jobs/default/<JOB_ID>/`.

### 3.2 Command Examples

#### Option A: Running from OpenCRAVAT Job ID
```bash
# Python Engine
python3 run_ontology_pipeline.py --sample DE_master --input 260706-105810

# Bash Wrapper
./run_ontology_master_pipeline.sh DE_master 260706-105810
```

#### Option B: Running directly from a VCF / VCF.GZ file
```bash
# Python Engine
python3 run_ontology_pipeline.py --sample HG003 --input /data/genomes/HG003.phased.vcf.gz

# Bash Wrapper
./run_ontology_master_pipeline.sh HG003 /data/genomes/HG003.phased.vcf.gz
```

#### Option C: Running from an annotated SQLite database
```bash
python3 run_ontology_pipeline.py --sample Patient_101 --input /path/to/Patient_101.sqlite
```

---

## 4. Output Deliverables & Google Drive Sync

The pipeline automatically creates a dated subfolder under your local `reports/` directory and synchronizes all deliverables to your Google Drive `"Ontology"` folder:

### 4.1 Output Path Structure
* **Local Workspace**: `./reports/{Sample_ID}-{DD-MM-YYYY}/`
* **Google Drive Target**: `~/Google Drive/My Drive/Ontology/{Sample_ID}-{DD-MM-YYYY}/`

### 4.2 Generated Deliverables
| File Name Pattern | Description |
| :--- | :--- |
| **`{Sample_ID}_visual_explorer.html`** | Standalone interactive Visual Ontology Explorer (D3 DAG hierarchy, non-redundant variant drawers, zoomable graph, interactive exons). |
| **`{Sample_ID}_master_ontology_report.html`** | Comprehensive Universal Master Hub Report with full domain breakdowns. |
| **`{Sample_ID}_master_actionable.json`** | Filtered, tiered actionable variants dataset with protective frequency bypass. |
| **`{Sample_ID}_variants.tsv`** | Tab-separated matrix of all clinical variant calls, in-silico scores, and phasing. |
| **`{Sample_ID}_summary.txt`** | Text summary for clinical review. |
| **`{Sample_ID}_report.pdf`** | Headless Chrome auto-generated high-resolution PDF document. |
| **`{Sample_ID}_iOS_bundle.zip`** | Portable offline archive containing all HTML, TSV, and JSON assets. |

---

## 5. Visual Ontology Explorer Interactive Features

1. **Expanding Variant Details**:
   * Click on the character/symbol (**`▸`** / **`▾`**) to the left of the **`rs`** identifier to expand or collapse the 4-box clinical drawer.
2. **Active rs Links**:
   * Clicking directly on any **`rsID`** opens the official NCBI dbSNP (or ClinVar) page in a new browser tab.
3. **Variant-Specific Literature & Studies**:
   * The expanded variant drawer displays direct links to **LitVar2**, **GWAS Catalog**, **UCSC Genome Browser (GRCh38)**, and associated PubMed articles referenced by that specific variant.
4. **Protective & Drug Response Filtering**:
   * Protective variants ($MAF \approx 0.10 - 0.76$) and CPIC/DPWG pharmacogenomic drug responses bypass the 1% rare disease ceiling and are categorized in green under **Protective Associations**.
5. **Phasing Badges**:
   * Direct extraction from VCF piped genotypes displays **Maternal** (`0|1`), **Paternal** (`1|0`), **Homozygous** (`1|1`), or **Unphased** badges.
