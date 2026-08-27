# Genomics Ontology Reporting Engine & Visual Explorer (v4.5)

An ontology-driven clinical genomics interpretation engine powered by **OpenCRAVAT 3.1.1**, **LinkML**, **Pydantic v2**, and modern **vanilla Web Standards**. 

This system bridges raw genomic variants with biomedical ontologies (**HPO**, **GO**, **Anatomical Organ Systems**), multi-predictor in-silico scores (**AlphaMissense**, **CADD**, **SpliceAI**, **REVEL**), polygenic risk scores (**PRS**), pharmacogenomics (**CPIC / DPWG**), and curated peer-reviewed literature (**PubMed / LitVar**).

---

## 📚 Core Documentation Links
- 📋 **[Pull Request Document](docs/PULL_REQUEST.md)**: Full clinical and engineering specifications, requirements, architecture, and verification results.
- 🚀 **[Brainstorming & Strategic Opportunities](docs/OPPORTUNITIES_BRAINSTORMING.md)**: Visionary blueprint covering Multimodal Audio/Visual AI, Agent Skills, 3D AlphaFold mapping, and single-cell epigenomics.
- 📖 **[User Guide & Pipelines](docs/User_Guide.md)**: Detailed instructions on running pipelines against OpenCRAVAT SQLite databases.

---

## ⚡ Highlights & Key Capabilities

### 1. Multi-Level Formal Ontology Lineage (4 Levels Deep)
- **Human Phenotype Ontology (HPO)**: Formal hierarchy (e.g. `HP:0001626` Cardiovascular System $\rightarrow$ `HP:0001627` Heart Morphology / `HP:0011025` Physiology $\rightarrow$ `HP:0001638` Cardiomyopathy / `HP:0011675` Arrhythmias $\rightarrow$ `HP:0001657` Long QT / `HP:0001663` Brugada syndrome $\rightarrow$ Associated Genes).
- **Gene Ontology (GO)**: Distinct roots for **Biological Process** (`GO:0008150`), **Molecular Function** (`GO:0003674`), and **Cellular Component** (`GO:0005575`).
- **Organ / System View**: True anatomical breakdown across 9 organ systems (Heart, Brain, Lungs, Skeleton, Immune, Kidneys, Digestive, Blood, Sensory Organs).

### 2. Gene-Grouped Variants with NCBI Summaries
- Under each clinical tier (*Potential Concerns*, *Protective*, *Uncertain*, *Baseline*), variants are organized into **Gene Cards** featuring:
  - Official gene symbols and full names.
  - NCBI / OMIM biological summary paragraphs.
  - Clickable reference links (**NCBI Gene**, **OMIM**, **GeneCards**, **ClinVar**).
  - Clean nested tables for all variants identified in that gene.

### 3. Non-Redundant Comprehensive Variant Detail Drawer
- Expanding a variant displays a 4-column clinical card:
  - **Genomic & Transcript**: Chromosome locus, cDNA change (`c.`), protein change (`p.`), transcript ID, and Variant Allele Fraction (VAF %).
  - **In-Silico Predictions**: AlphaMissense score & class, CADD Phred, REVEL ensemble score, and all four SpliceAI delta scores (**AG, AL, DG, DL**).
  - **ClinVar & ACMG Evidence**: Classification, review status, and ACMG PM5 / PS1 hotspot codes.
  - **Research & GWAS Studies**: Study title, statistical findings, odds ratio/beta, p-value, and risk allele.

### 4. Interactive Scalable Graph View with Zoom & Pan
- High-contrast SVG tree graph with smooth cubic curves.
- Zoom toolbar (**Zoom In `+`**, **Zoom Out `−`**, **Reset `↺`**) with CSS matrix scaling to inspect all levels without clipping.

### 5. Multi-System Genomic Risk Profile & Pharmacogenomics
- 8-system risk matrix summarizing risk tiers (**HIGH**, **MODERATE**, **TYPICAL**), primary affected biological pathways, polygenic risk percentiles, and high-risk concern genes.
- Actionable CPIC / DPWG pharmacogenomics guidance table with diplotypes and dosing recommendations.

### 6. Peer-Reviewed Bibliography Integration
- 67+ curated PubMed citations across clinical loci with direct NCBI links.

### 7. Artistic Epigenetic Watermark Background
- Faint SVG background watermark featuring glowing DNA helices, histone octamers, and methylation tags.

---

## 📂 Repository Structure

```
genomics-ontology/
├── docs/
│   ├── PULL_REQUEST.md                 # Formal Pull Request & technical specification
│   ├── OPPORTUNITIES_BRAINSTORMING.md  # Visionary roadmap (Audio AI, Skills, 3D AlphaFold)
│   ├── PLAN.md                         # Architecture notes and initial designs
│   └── User_Guide.md                   # Operational guide
├── generate_claude_v2_report.py        # Python ETL pipeline for generating DAG JSON data
├── index.html                          # 5-view web application shell
├── js/
│   └── app.js                          # Client-side reactive router, tree, and graph engine
├── css/
│   └── style.css                       # Clinical design tokens, responsive grids, and print CSS
├── data/
│   └── mock-data.js                    # Verified data payload (702 variants, 577 genes)
├── reports/                            # Generated HTML clinical reports
│   └── DE_master_260706/
│       ├── DE_master_master_actionable.json
│       ├── DE_master_master_ontology_report.html
│       └── claude_v2_explorer.html
├── test_claude_v2.js                   # Automated JSDOM validation suite
├── run_ontology_master_pipeline.sh     # Master bash pipeline runner
├── serve_reports.sh                    # Local HTTP server launcher
└── package.json                        # Testing dependencies
```

---

## 🚀 Quick Start

### 1. Launch the Visual Explorer
```bash
# Start local HTTP server
python3 -m http.server 8082 --directory .

# Open in your browser:
# http://localhost:8082/index.html
```

### 2. Run Automated Verification Tests
```bash
NODE_PATH=node_modules node test_claude_v2.js
```

### 3. Re-generate Report Data from OpenCRAVAT SQLite/JSON
```bash
python3 generate_claude_v2_report.py \
  reports/DE_master_260706/DE_master_master_actionable.json \
  data/mock-data.js
```

---
