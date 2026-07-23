# PLAN — Domain-Agnostic Ontology-Driven Actionable Variant Report

**Project:** `ontology_report` (formerly `cardio_ontology_report`)
**Location:** `~/My-Projects/genomics/development/ontology_report`
**Status:** Living document — updated at each increment (plan → implement → validate → refactor → validate).

> **Scope note (I7):** The engine is now domain-agnostic. Cardiology is one
> example config; any phenotype domain is defined by a config file under
> `config/` (see `hereditary_cancer.yaml`, `template.yaml`). Nothing below is
> cardiology-specific except the worked cardiology example.

---

## 1. Motivation & Difference From the Previous Pipeline

The previous engine (`../cardiology_pipeline`) selects variants with a **hard-coded gene list**
(~60 genes) embedded in `run_sql_filter.py`. This has three problems:

1. **Static** — new gene–disease relationships require editing SQL.
2. **Opaque** — the clinical rationale for each gene is not encoded; a gene is either "in" or "out".
3. **Not phenotype-driven** — it cannot answer "show me variants in genes linked to *arrhythmia*"
   without a human curating the list.

The new engine selects genes **dynamically from ontologies already installed in OpenCRAVAT**:

- **HPO** (Human Phenotype Ontology, `hpo` annotator) — gene → clinical phenotype terms.
- **GO** (Gene Ontology, `go` annotator) — gene → molecular function / biological process / cellular component.

We define a small **seed set of cardiac ontology terms** (HPO IDs + GO IDs + keyword patterns).
The gene panel is *derived* by asking the ontologies "which genes are annotated with these terms?".
Editing the panel now means editing a phenotype/function list, not a gene list.

## 2. Actionability Definition

The report focuses only on **actionable variants** — variants for which there is real
phenotype and/or genotype information. A variant is *actionable* if it passes the ontology
gene gate **and** has at least one line of evidence:

**Phenotype (P) evidence — the gene/variant is tied to disease**
- ClinVar significance is Pathogenic / Likely pathogenic / VUS (`clinvar__sig`).
- ClinGen gene–disease validity present (`clingen__classification`).
- OMIM disease entry present (`omim__*`).
- Gene has ≥1 HPO disease-phenotype term (from the ontology gate itself).

**Genotype (G) evidence — the variant is plausibly functional**
- Coding consequence: missense / nonsense / frameshift / start-loss / canonical splice (`base__so`).
- Predictor support above threshold: REVEL, AlphaMissense, CardioBoost, BayesDel, MetaRNN, ESM1b, VARITY.
- SpliceAI Δ ≥ 0.20 (splice-altering).
- Non-coding regulatory support in a panel gene: CADD phred, LINSIGHT, ncER, RegulomeDB, cCRE-SCREEN.
- Rare in population: gnomAD4 AF and All-of-Us AF below tier thresholds.

Variants with **no** phenotype and **no** genotype evidence are dropped (they are not actionable).

## 3. Tiering

| Tier | Meaning | Rule (simplified) |
|------|---------|-------------------|
| **Tier 1 — Reportable** | Pathogenic-grade | ClinVar P/LP, OR PVS1-type LoF in a haploinsufficient panel gene, OR strong multi-predictor consensus + rare |
| **Tier 2 — VUS of interest** | Uncertain but supported | ClinVar VUS, OR ≥1 predictor pathogenic + rare, OR splice-altering |
| **Tier 3 — Monitor** | Weak/regulatory | Non-coding regulatory support only, or single weak signal |
| **Filtered** | Not actionable | Fails ontology gate or has no evidence |

Each kept variant carries **reason codes** (e.g. `HPO_CARDIOMYOPATHY`, `GO_HEART_CONTRACTION`,
`CLINVAR_PLP`, `PVS1_HAPLOINSUFFICIENT`, `PP3_MULTI_PREDICTOR`, `PM2_RARE`, `SPLICEAI_HIGH`)
so the selection is fully explainable.

## 4. Architecture

```
                 TEST_master_final.vcf.gz
                          │
              ┌───────────▼─────────────┐
              │ Stage 1: oc run          │  OpenCRAVAT annotation
              │ (go hpo clinvar clingen  │  → TEST.sqlite (raw)
              │  omim gnomad4 predictors)│
              └───────────┬─────────────┘
                          │
              ┌───────────▼─────────────┐
              │ Stage 2: build panel     │  lib/build_ontology_panel.py
              │ HPO seeds + GO seeds     │  → panel.json (gene set + why)
              │ (reads hpo.sqlite,       │
              │  go.sqlite ontology data)│
              └───────────┬─────────────┘
                          │
              ┌───────────▼─────────────┐
              │ Stage 3: schema probe    │  lib/schema_probe.py
              │ detect present columns   │  → schema.json
              └───────────┬─────────────┘
                          │
              ┌───────────▼─────────────┐
              │ Stage 4: actionable      │  lib/ontology_filter.py
              │ filter + tier + reasons  │  → <prefix>_actionable.sqlite
              └───────────┬─────────────┘
                          │
              ┌───────────▼─────────────┐
              │ Stage 5: enrich          │  lib/enrich_report.py (cached,
              │ NCBI Gene desc + GWAS     │  offline-safe) → enriched JSON
              │ Catalog study evidence   │
              └───────────┬─────────────┘
                          │
              ┌───────────▼─────────────┐
              │ Stage 6: render          │  lib/render_report.py (generic) or
              │ HTML + TSV + text +      │  lib/render_autoimmune.py (viz + live
              │ native Excel/VCF export  │  studies) + oc report
              └──────────────────────────┘
```

### Why this leverages OpenCRAVAT modules
- **Gene selection** uses the installed `hpo` and `go` annotator data SQLite files directly
  (`/data/opencravat/modules/annotators/{hpo,go}/data/*.sqlite`) — no external downloads.
- **Gene–disease validity** uses `clingen` + `omim` + `clinvar_acmg` annotators.
- **Report rendering** uses the installed OpenCRAVAT **reporters** (`tsvreporter`,
  `excelreporter`, `vcfreporter`, `textreporter`) via `oc report ... --filtersql`,
  in addition to our own styled HTML.

## 5. Data Model Notes (discovered during exploration)

- `hpo` and `go` are **gene-level** annotators → their columns live in the **`gene`** table,
  not `variant`. Report-time OpenCRAVAT joins gene-level onto variant-level.
- Relevant `gene` columns: `hpo__id`, `hpo__term`, `hpo__all`, `go__bpo_name/_id`,
  `go__cco_name/_id`, `go__mfo_name/_id`.
- Ontology source tables:
  - `hpo.sqlite` → `genes(hugo, hpo_id, hpo_term)` (329k rows).
  - `go.sqlite` → `go_annotation(hugo, go_id, go_aspect, qualifier, evidence_code)` (400k rows)
    and `go_name(go_id, name)`.
- Variant table join key: `base__uid`; gene join key: `base__hugo` ↔ `gene.base__hugo`.

## 6. Incremental Roadmap

| Increment | Goal | Validation |
|-----------|------|-----------|
| **I0** | Explore env, modules, schema (DONE) | `oc module ls`, schema dumps |
| **I1** | Ontology panel builder from HPO+GO seeds | panel gene count sane, contains known cardiac genes (MYH7, KCNQ1, TTR…) |
| **I2** | Schema probe + actionable filter + tiering | run on Melinda DB, non-zero tiered variants, reason codes present |
| **I3** | Report renderer (HTML/TSV/text) + orchestrator | reports open, counts match filter |
| **I4** | Offline test suite (mock sqlite) | `python3 tests/test_flow.py` green |
| **I5** | Full validation on TEST DB (after annotation) | end-to-end run, deliverables produced |
| **I6** | Refactor for clarity/perf; re-validate | tests still green, TEST run reproducible |

## 7. Iteration Strategy for Cost Control

Full annotation of the 226k-variant TEST VCF takes ~4.5 h (CADD/ncER/gnomAD are the bottleneck).
To keep the plan→implement→validate loop fast, tooling is developed and validated against an
**existing fully-annotated database** (`Melinda_Ehrle.sqlite`, 5.2M variant rows, same annotator
set) while the TEST annotation runs once in the background. Final validation (I5) uses the
freshly annotated `TEST.sqlite`.

## 8. Deliverables

- `config/<domain>.yaml` — per-domain seed HPO/GO terms + thresholds (the only file to edit to retarget); ships `cardiology.yaml`, `hereditary_cancer.yaml`, `autoimmunity.yaml`, `template.yaml`.
- `lib/*.py` — panel builder, schema probe, filter, enrichment (NCBI + GWAS), generic renderer, autoimmune renderer.
- `run_ontology_report.sh` — one-command orchestrator (7 stages; `-c/-o/-E`).
- `tests/` — offline mock validation (`test_flow.py`, `test_autoimmune.py`) + a runnable sample-report demo (`demo_autoimmune.py`).
- `README.md`, `docs/User_Guide.md`, `docs/PLAN.md` (this file), plus per-stage change notes.

## 9. Change Log
- **2026-07-07 I0:** Environment/module/schema exploration complete. TEST annotation started in
  background. Directory scaffolded. Plan authored.
- **2026-07-07 I1:** `build_ontology_panel.py` complete. Panel = 2,362 genes from HPO+GO;
  all classic cardiac genes verified present.
- **2026-07-07 I2:** `schema_probe.py` + `ontology_filter.py` complete. Actionability gate
  tightened after first pass was too loose (Melinda 481→378 kept, sensible tiers).
- **2026-07-07 I3:** `render_report.py` + `run_ontology_report.sh` complete. Fixed multi-kb
  allele blow-up in text report. Native OC Excel/VCF export wired in (Stage 6).
- **2026-07-07 I4:** `tests/test_flow.py` — 11/11 decision-path checks pass.
- **2026-07-07 I5:** Full TEST run validated (19 actionable; HTML verified in-browser).
  Documented SV-skip converter limitation (152,834 symbolic-ALT lines dropped).
- **2026-07-07 I6:** Refactored gene-level column handling to be fully schema-aware;
  re-validated on both Melinda (no clingen/omim) and TEST (with clingen/omim). Tests green.
  Detailed evidence in [`STEP_NOTES.md`](STEP_NOTES.md).
- **2026-07-09 I8 (gene/genotype context + live-evidence autoimmunity report):**
  Two report-wide field additions and a new domain that showcases "live"
  reporting.
  - **Zygosity + rsID:** `schema_probe` now detects genotype/zygosity (single-
    sample `vcfinfo__*` and the multi-sample `sample` table) and dbSNP rsID;
    `ontology_filter` normalizes zygosity (Het/Hom/Hemi), derives variant allele
    fraction from read depths, and carries rsID through. Shown on every card, in
    TSV and text (all domains).
  - **NCBI Gene description:** new `lib/enrich_report.py` fetches per-gene
    NCBI Gene descriptions/summaries (E-utilities) and current GWAS Catalog study
    associations (EBI REST). Cached, rate-limited, offline-safe; a new orchestrator
    stage (5/7) with `-o` (offline) / `-E` (skip) flags and config-driven toggles.
  - **`bypass_frequency` domain evidence:** a `domain_evidence` entry can now keep
    population-common variants (e.g. GWAS risk alleles), tagged `RISK_ALLELE_COMMON`.
  - **Autoimmunity domain:** `config/autoimmunity.yaml` (autoimmune HPO + immune GO
    seeds, HLA/PTPN22/CTLA4/… force-includes, `gwas_catalog` bypass, live enrichment
    on) + `lib/render_autoimmune.py` (inline SVG trait-burden visualization + per-
    variant live GWAS study-evidence tables + NCBI gene context). Renderer is chosen
    by the config's `report.renderer`.
  - **Validation:** `tests/test_autoimmune.py` (offline, pre-seeded cache) green;
    `tests/demo_autoimmune.py` builds a full sample report without OpenCRAVAT.
    End-to-end online demo confirmed 12/12 records enriched with real NCBI
    descriptions + GWAS associations (e.g. PTPN22/rs2476601 → T1D p≈2e-80). Generic
    `tests/test_flow.py` still 12/12.

- **2026-07-07 I7 (generalization):** De-coupled the engine from cardiology. Renamed
  project `cardio_ontology_report` → `ontology_report`. Config is now per-domain
  (`config/cardiology.yaml`, `config/hereditary_cancer.yaml`, `config/template.yaml`);
  orchestrator takes `-c CONFIG`; output filenames namespaced by domain. Removed
  hard-coded CardioBoost/ArrVars from the filter — they are now generic, optional
  `domain_predictors` / `domain_evidence` config entries (absent columns ignored).
  Renamed the output column `cardio_tier` → `tier`; renderer title is config-driven.
  **Validation:** 12/12 offline checks pass (incl. a config-driven domain-evidence
  case); on the TEST DB, cardiology kept 19 (panel 2,362) and hereditary cancer kept
  8 (panel 2,303, BRCA1/BRCA2/MLH1/… present) from the *same* database — confirming
  the domain is fully config-driven.
