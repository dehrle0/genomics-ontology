# User Guide: Genomic Ontology Reporting & Clinical Pipeline (v5.3)

Welcome to the **Genomic Ontology Reporting & Clinical Pipeline**. This platform translates Whole Genome Sequencing (WGS), Whole Exome Sequencing (WES), and targeted panel variant calls into formal multi-level Directed Acyclic Graph (DAG) biomedical ontologies (**HPO**, **GO**, and **Anatomical Organ/System** views), integrating in-silico predictors, polygenic risk scores (PRS), pharmacogenomic recommendations (CPIC/DPWG), Google DeepMind AlphaGenome Atlas annotations, and publication-grade **Clinician Action Briefs**.

---

## 1. Environment & Prerequisites

Activate your conda, mamba, or virtual environment containing Python 3.10+:

```bash
micromamba activate cravat_env
```

### 1.1 Python Dependencies
Install all required dependencies using `requirements.txt`:
```bash
pip install -r requirements.txt
```

### 1.2 Automated Zero-Curation Guard
Before executing pipeline scripts, verify that no static PII or variant curations exist in code or skills:
```bash
python3 lib/zero_curation_guard.py
# Must output: ZERO_CURATION_GUARD_PASSED
```

---

## 2. Standard Reporting Pipeline (`run_standard_reporting.py`)

The primary orchestrator for generating clinical deliverables is `run_standard_reporting.py`. It reads single-sample annotated JSONs, applies ACMG BA1 population frequency and synonymous gating, generates the canonical 3-Page Minimalist Clinician Action Brief, compiles the Deep Clinical Research Dossier, renders PDFs via headless Chrome, and synchronizes deliverables to `export/` and Google Drive.

### 2.1 Standard Execution
```bash
# Generate deliverables for sample DE
python3 run_standard_reporting.py de

# Generate deliverables for sample ME
python3 run_standard_reporting.py me

# Explicit dated sample directory override
python3 run_standard_reporting.py de --sample-dir reports/Daniel_Ehrle-09-10-2026

# Validate inputs without writing files (dry-run)
python3 run_standard_reporting.py de --dry-run
```

---

## 3. Output Deliverables & Hierarchy

The reporting pipeline generates synchronized deliverables stored in `./reports/{Sample_ID}-{Date}/`, mirrored to `./export/`, and synced to `~/Google Drive/My Drive/Ontology/`:

### 3.1 Core Deliverables Matrix
| Deliverable | File Pattern | Pages / Format | Clinical Purpose |
| :--- | :--- | :--- | :--- |
| **Canonical Clinician Action Brief** | `{Sample_ID}_clinical_brief.html` / `.pdf` | **Strictly 3 Pages** (Letter Portrait) | **Primary Clinical Brief:** Minimalist layout covering Priority Findings, Medication Precautions (fatal-only Life-Threatening badges), Preventive Care, PGx Reference, and EHR Note. |
| **Deep Clinical Research Dossier** | `{Sample_ID}_deep_research_report.md` / `.html` / `.pdf` | Multi-Page Research Synthesis | **In-Depth Evidence:** Mechanistic causality vs. statistical correlation, multi-engine AI scores (AlphaGenome, AlphaMissense, CADD, REVEL), and cis/trans phasing. |
| **FHIR Patient EHR Import** | `{Sample_ID}_ehr_import.json` | JSON Data Bundle | **EHR Integration:** Machine-readable clinical decision support alerts, structured phenotypes, and surveillance recommendations. |
| **Standalone Landscape Dashboard** | `{Sample_ID}_landscape_dashboard.html` | Single-File HTML5 (Offline) | **Interactive Explorer:** Offline visual explorer with organ system ontology navigation, AlphaGenome candidate heatmaps, and SVG trees. |
| **Deprecated Baseline Brief** *(Archival)* | `{Sample_ID}_clinical_brief (old).html` / `.pdf` | 3 Pages (Letter Portrait) | **Legacy Baseline:** Older verbose baseline brief preserved for historical comparisons and audit trail. Deprecated in favor of the canonical brief. |

---

## 4. Canonical 3-Page Clinician Action Brief Structure

The canonical brief is optimized for clinical readability, immediate comprehension, and zero printing truncation:

* **Page 1: Priority Genomic Findings & Carrier Status**
  - Minimalist header with patient identifier, date, and audience routing.
  - Section 1 table: *Gene & Biological System*, *Variant & Zygosity*, *Classification & Condition*, and *Clinical Impact & Guidance*.
  - Strict text limit: 1–2 sentences (~80–140 chars) per row with explicit single-copy carrier reassurance. Capped at 7 prioritized rows.
* **Page 2: Medication Precautions & Preventive Care**
  - Section 2 table: High-Priority Medication & Exposure Precautions.
    - Badges labeled `🚨 LIFE-THREATENING RISK` are strictly gated to fatal/life-threatening hazards (e.g. fatal hepatic necrosis, lethal myelotoxicity, severe thromboembolism).
    - Other high risks are styled with bold amber badges (`⚠️ HIGH RISK · AVOID`, `⬇️ HIGH RISK · REDUCE DOSE`). Capped at 5 rows.
  - Section 3 table: High-Level Preventive Care & Annual Checkup Discussion Points. Capped at 5 rows.
* **Page 3: Pharmacogenomic Guide, Clinical Perspective & EHR Note**
  - Section 4 table: Actionable Pharmacogenomic Prescribing Guide (CPIC / DPWG) in 4 columns, capped at 5 key actionable pharmacogenes.
  - Section 5: Clinical Perspective side-by-side cards (Proactive Care Utility vs. Avoiding Over-Intervention).
  - Section 6: Copy-Paste Physician Handoff / EHR Summary Note.
  - Technical Telemetry Footnote: Displays empirical BAM autosomal mean alignment depth and median locus DP.

---

## 5. Running OpenCRAVAT (OC) Annotation (Pre-computation)

To generate the rich multi-annotator SQLite database from raw VCF files:

### 5.1 Required Annotators (25 Essential Annotators)
* **Clinical & Pathogenicity:** `clinvar`, `clingen`, `civic`, `clinvar_acmg`
* **Ontologies:** `hpo`, `go`
* **Gene Models & Disease Databases:** `ncbigene`, `omim`, `interpro`
* **In-Silico Predictors:** `revel`, `alphamissense`, `cadd`, `spliceai`, `bayesdel`, `metarnn`, `esm1b`, `varity`
* **Non-Coding & Regulatory:** `linsight`, `ncer`, `regulomedb`, `ccre_screen`
* **Population & Frequency:** `gnomad4`, `dbsnp`, `vcfinfo`, `gtex`
* **Pharmacogenomics & GWAS:** `pharmgkb`, `gwas_catalog`

### 5.2 Command to Run OpenCRAVAT on a VCF
```bash
oc run /path/to/input_phased.vcf.gz \
  -l hg38 \
  -a hpo go clinvar clingen omim ncbigene revel alphamissense bayesdel metarnn esm1b varity spliceai cadd linsight ncer regulomedb ccre_screen gtex dbsnp vcfinfo gwas_catalog pharmgkb civic interpro \
  -d /path/to/output_dir \
  --mp $(nproc) \
  -n Sample_Prefix
```

---

## 6. Pre-Commit Verification Checklist

Before pushing changes to GitHub:
1. `python3 lib/zero_curation_guard.py` $\rightarrow$ Must output `ZERO_CURATION_GUARD_PASSED`.
2. `pdfinfo reports/<Sample_Dir>/<Sample>_clinical_brief.pdf` $\rightarrow$ Must report `Pages: 3`.
3. Check `export/` and Google Drive directories to confirm all files are synchronized.
