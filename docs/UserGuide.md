# UserGuide.md: Genomic & Pharmacogenomic Ontology Reporting

Welcome to the user guide for running the Unified Ontology and Pharmacogenomics Reporting Tool. This tool maps Whole Genome Sequencing (WGS) variant calls to core biomedical ontologies (HPO, MONDO, EFO/PGS) and translates star-allele diplotypes to clinical prescribing recommendations (PharmCAT).

---

## 1. Setup & Environment
Ensure you have activated your micromamba `cravat_env` environment containing the necessary local OpenCRAVAT databases:
```bash
micromamba activate cravat_env
```

### 1.1 Python Dependencies
Install the required packages locally:
```bash
pip install pyyaml openpyxl pydantic linkml-runtime
```

---

## 2. Running the Core Pipeline
The tool runs dynamically based on domain configurations (e.g. `config/cardiology.yaml` or `config/autoimmunity.yaml`).

### 2.1 Standard Command Line Interface (CLI)
```bash
# General usage
./run_ontology_report.sh -c <config.yaml> <input_sqlite_or_vcf> <output_dir> <patient_prefix>

# Running cardiology report on an annotated SQLite database:
./run_ontology_report.sh -c config/cardiology.yaml \
  /data/Genomes/TEST/Data/Final/2026-03-22/ontology_reports/TEST.sqlite \
  /data/Genomes/TEST/Data/Final/2026-03-22/ontology_reports \
  TEST
```

---

## 3. Polygenic Trait Roll-ups vs. Monogenic Tiers
- **Monogenic Variants**: ClinVar and REVEL annotated variants are tiered (Tier 1/2/3) and mapped under HPO Level 1 (organ system) nodes.
- **Polygenic Roll-ups**: Polygenic variants are mathematically integrated (using raw logs or odds-ratio multipliers) to compute a single **trait percentile** matching your demographic background. Trait percentiles are visualized under HPO physiological subcategories using horizontal bar charts.

---

## 4. Google Drive MCP Integration
Automated delivery utilizes a Model Context Protocol (MCP) server connected to your personal Google Drive account.
The script `cloud_delivery_service.py` is called at the end of the pipeline run to securely transfer files directly to the `"Ontology"` (or `"Pharmacogenomics"`) drive folder.
To configure, provide your OAuth2 client secret file at `config/gdrive_credentials.json` and authenticate on your first run.
