# Genomic Ontology Reporting Engine & Visual Explorer (v5.2)

An ontology-driven clinical genomics interpretation engine powered by **OpenCRAVAT 3.1.1**, **LinkML**, **Pydantic v2**, and modern **vanilla Web Standards**.

This system bridges raw genomic variants with biomedical ontologies (**HPO**, **GO**, **Anatomical Organ Systems**), multi-predictor in-silico scores (**AlphaMissense**, **CADD**, **SpliceAI**, **REVEL**, **LINSIGHT**), whole-genome pedigree phasing (**Approach 4D Trio/Duo**), polygenic risk scores (**PRS**), pharmacogenomics (**CPIC / DPWG**), and automated **Google DeepMind AlphaGenome Atlas** candidate triage.

---

## 📚 Core Documentation Links
- 📋 **[Pull Request Document](docs/PULL_REQUEST.md)**: Full clinical and engineering specifications, requirements, architecture, and verification results.
- 🚀 **[Brainstorming & Strategic Opportunities](docs/OPPORTUNITIES_BRAINSTORMING.md)**: Visionary blueprint covering Multimodal Audio/Visual AI, Agent Skills, 3D AlphaFold mapping, and single-cell epigenomics.
- 📖 **[User Guide & Pipelines](docs/User_Guide.md)**: Detailed instructions on running the 7-stage orchestrator against OpenCRAVAT SQLite databases and phased VCFs.

---

## ⚡ Highlights & Key Capabilities (v5.2)

### 1. Unified 7-Stage End-to-End Orchestrator (`run_ontology_pipeline.py`)
- Automated single-command CLI executing the entire clinical reporting lifecycle:
  1. **Input Resolution & Validation**: Direct inspection of OpenCRAVAT SQLite and phased VCF inputs.
  2. **Multi-Domain Panel Construction**: HPO + GO + Organ panel synthesis (9,038 panel genes).
  3. **Database Schema & Column Probing**: Dynamic mapping of 239 variant annotation fields.
  4. **Actionable Variant Filtering & Phased VCF Integration**: High-performance multi-tier filtration with pedigree haplotype phasing.
  5. **Gene & Literature Enrichment**: Cached fallback PubMed/LitVar and OMIM gene summaries.
  6. **Multi-Portal Rendering**: Standalone single-file HTML5 Visual Explorer with embedded data, Universal Master Portal, and D3 graph engine.
  7. **Deliverables Packaging & Dual Cloud Sync**: Automated generation of PDF reports, candidate matrices, and synchronization to local Google Drive and cloud remote (`rclone`).

### 2. Whole-Genome Pedigree Phasing Integration (Approach 4D)
- Directly ingests and correlates whole-genome phased VCFs (SNV, Indel, SV, CNV, and STR):
  - **Trio / Duo Informative Anchors**: Differentiates maternal allele inheritance (`0|1`, e.g., `(SE Anchor)` or `(MI Anchor)`) and paternal allele inheritance (`1|0`).
  - **Haplotype Block Context**: Assigns physical block IDs and strand orientations for phase-set-aware compound heterozygosity detection.

### 3. DeepMind AlphaGenome Atlas Candidate Triage & Rescue
- Automatically isolates variants requiring 9-billion-variant precomputed AlphaGenome resolution:
  - **ClinVar Discordance (`CLINVAR_CONFLICT`)**: Targets variants with conflicting classifications of pathogenicity among submitters.
  - **Predictor Discordance (`PREDICTOR_DISCORDANCE`)**: Highlights variants where protein-language models (AlphaMissense, ESM1b) conflict with ensemble predictors (REVEL, BayesDel) or splice tools.
  - **Ultra-Conserved Non-Coding Rescue (`ULTRA_CONSERVED_NONCODING`)**: Rescues ultra-rare non-coding variants ($AF < 10^{-4}$ or unobserved) in disease-panel genes displaying extreme vertebrate negative selection (`LINSIGHT >= 0.80`) or top 1% deleteriousness (`CADD Phred >= 20.0`).
- **Dedicated Candidate Matrix TSV**: Automatically exports `*_alphagenome_candidates.tsv` complete with pre-formatted, 1-click DeepMind Atlas exploration deep-links.

### 4. Standalone Single-File HTML5 Visual Explorer
- Self-contained, zero-dependency HTML5 application (10+ MB) with inlined styles, scripts, and complete JSON datasets:
  - **Offline Portability**: Runs natively in any browser without requiring an active backend or node process.
  - **4-Level Ontology Hierarchy**: Interactive navigation across HPO, GO (Biological Process, Molecular Function, Cellular Component), and 9 Anatomical Organ Systems.
  - **Integrated Pan & Zoom Tree View**: Interactive SVG graph engine with smooth cubic Bezier curves and full matrix transformation controls.

### 5. Multi-System Genomic Risk Profile & Pharmacogenomics
- 8-system risk matrix summarizing risk tiers (**HIGH**, **MODERATE**, **TYPICAL**), primary affected biological pathways, polygenic risk percentiles, and high-risk concern genes.
- Actionable CPIC / DPWG pharmacogenomics guidance table with diplotypes and dosing recommendations.

---

## 📂 Repository Structure

```
ontology_report/
├── run_ontology_pipeline.py            # Master 7-stage Python pipeline orchestrator (v5.2)
├── generate_claude_v2_report.py        # Python ETL pipeline for generating DAG JSON data
├── index.html                          # 5-view web application shell
├── js/
│   └── app.js                          # Client-side reactive router, tree, and graph engine
├── css/
│   └── style.css                       # Clinical design tokens, responsive grids, and print CSS
├── lib/
│   ├── ontology_filter.py              # Clinical multi-tier filtering & AlphaGenome triage engine
│   ├── phased_vcf_integrator.py        # Whole-genome phased VCF integration engine
│   ├── render_master_hub.py            # HTML5 Master Hub generator
│   ├── render_report.py                # TSV and HTML clinical report renderer
│   ├── render_autoimmune.py            # Specialized SVG trait-burden chart & clinical tables
│   ├── enrich_report.py                # Gene annotation & PubMed literature enrichment
│   └── schema_probe.py                 # OpenCRAVAT SQLite schema probe
├── data/
│   └── mock-data.js                    # Verified data payload
├── docs/
│   ├── PULL_REQUEST.md                 # Formal Pull Request & technical specification
│   ├── OPPORTUNITIES_BRAINSTORMING.md  # Visionary roadmap (Audio AI, Skills, 3D AlphaFold)
│   ├── PLAN.md                         # Architecture notes and initial designs
│   └── User_Guide.md                   # Operational guide
└── package.json                        # Testing dependencies
```

---

## 🚀 Quick Start

### 1. Run the Unified End-to-End Pipeline

Execute the full 7-stage pipeline against an OpenCRAVAT SQLite database and phased VCFs:

```bash
# Daniel Ehrle (with SE Maternal Anchor):
python3 run_ontology_pipeline.py \
  --sample Daniel_Ehrle \
  --input reports/Daniel_Ehrle-02-10-2026/Daniel_Ehrle.sqlite \
  --phased-vcf /data/Genomes/DE/Approach4D_Output/DE_grch38_wgs_phased.pass.vcf.gz,/data/Genomes/DE/Approach4D_Output/DE_grch38_sv_phased.vcf.gz,/data/Genomes/DE/Approach4D_Output/DE_grch38_cnv_phased.vcf.gz,/data/Genomes/DE/Approach4D_Output/DE_grch38_str.vcf.gz

# Melinda Ehrle (with MI Maternal Anchor):
python3 run_ontology_pipeline.py \
  --sample Melinda_Ehrle \
  --input reports/Melinda_Ehrle-03-10-2026/Melinda_Ehrle.sqlite \
  --phased-vcf /data/Genomes/ME/Approach4D_Output/ME_grch38_wgs_phased.pass.vcf.gz,/data/Genomes/ME/Approach4D_Output/ME_grch38_sv_phased.vcf.gz,/data/Genomes/ME/Approach4D_Output/ME_grch38_cnv_phased.vcf.gz,/data/Genomes/ME/Approach4D_Output/ME_grch38_str.vcf.gz
```

### 2. Launch Local Web Server & Inspect Reports

```bash
# Serve reports directory:
python3 -m http.server 8080 --directory reports

# Open in your browser:
# http://localhost:8080/Daniel_Ehrle-03-10-2026/Daniel_Ehrle_visual_explorer.html
# http://localhost:8080/Melinda_Ehrle-03-10-2026/Melinda_Ehrle_visual_explorer.html
```

### 3. Review Prioritized AlphaGenome Atlas Candidates

Each run outputs a dedicated TSV containing coordinates, gene symbols, consequences, in-silico scores, and direct 1-click AlphaGenome Atlas links:

```bash
# View top prioritized candidate variants:
head -n 20 reports/Daniel_Ehrle-03-10-2026/Daniel_Ehrle_alphagenome_candidates.tsv
head -n 20 reports/Melinda_Ehrle-03-10-2026/Melinda_Ehrle_alphagenome_candidates.tsv
```

---

## 🔒 Automated Cloud Delivery

Deliverables are automatically mirrored and synchronized to Google Drive:
- **Local Mirror**: `/home/daniel-ehrle/Google Drive/My Drive/Ontology/<Sample>-<Date>/`
- **Cloud Remote**: `drive:Ontology/<Sample>-<Date>/` via automated background `rclone` sync.
