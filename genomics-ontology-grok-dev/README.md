# Genomics Ontology Reporting Engine

An ontology-driven, actionable variant-to-phenotype mapping engine powered by **OpenCRAVAT**, **LinkML**, and **Pydantic v2**. This project dynamically builds clinical panels and roll-ups by querying biomedical ontologies (HPO, MONDO, EFO) and pharmacogenomics databases (PharmGKB / CPIC / PharmCAT).

---

## 📂 Project Repository Layout
```
ontology_report/
├── run_ontology_report.sh          # Orchestrator bash script
├── genomic_ontology_schema.yaml    # Formal LinkML schema specification
├── render_new_ontology_report.py   # Publication-grade interactive HTML/JS report renderer
├── cloud_delivery_service.py       # Google Drive upload and email notification webhook script
├── genomics_ontology_io/
│   ├── __init__.py
│   └── models.py                   # Pydantic v2 classes generated from LinkML schema
├── config/
│   ├── cardiology.yaml             # Cardiology domain seeds & rules
│   ├── autoimmunity.yaml           # Autoimmunity polygenic seeds with GWAS bypass
│   └── template.yaml               # scaffold template for custom domains
├── docs/
│   ├── PLAN.md                     # Architecture, design decisions, and Agent logs
│   └── UserGuide.md                # Environment setup and run instructions
└── tests/
    └── test_end_to_end.py          # Simulated mock SQLite validation suite
```

---

## 🚀 Highlights & Features
1. **Clinical Morphology vs. Physiology (HPO Level 2)**: Reports are grouped logically by Level 1 Organ Systems (e.g. Cardiovascular, Immune, Nervous) and partitioned into structural anomalies (Morphology) or functional defects (Physiology/Immunodeficiency).
2. **Polygenic Score Roll-ups**: Multi-locus additive polygenic scores are calculated and rolled up into clinical trait percentiles (rather than listed as raw disconnected SNPs).
3. **Interactive Visual Dashboard**: Incorporates horizontal bar charts, search/filter bars, and bold NCBI Gene descriptions above nested variant lists.
4. **Offline-First cached Enrichment**: Automatically retrieves live NCBI Gene summaries and GWAS study correlations, keeping runs offline-safe with local file caches.
