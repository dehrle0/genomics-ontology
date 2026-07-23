# ontology_report

**A domain-agnostic, ontology-driven, actionable variant report built on OpenCRAVAT.**

The engine knows nothing about any specific disease. A **domain** (cardiology,
hereditary cancer, neurodevelopmental disorders, …) is defined entirely by a
**config file** that lists the HPO/GO ontology terms used to *derive* the gene
panel and the thresholds that define "actionable". Swap the config → retarget
the whole report. No gene lists are hard-coded anywhere.

```
Genes are NOT listed by hand.
They are SELECTED by asking the ontologies:
  "Which genes does HPO link to <these phenotypes>?"
  "Which genes does GO link to <these functions>?"
```

---

## Quick start

```bash
micromamba activate cravat_env
cd ~/My-Projects/genomics/development/ontology_report

# Cardiology (default config) from an already-annotated OpenCRAVAT DB:
./run_ontology_report.sh -c config/cardiology.yaml \
  /data/Genomes/TEST/Data/Final/2026-03-22/ontology_reports/TEST.sqlite \
  /data/Genomes/TEST/Data/Final/2026-03-22/ontology_reports \
  TEST

# Hereditary cancer — SAME data, SAME engine, different config:
./run_ontology_report.sh -c config/hereditary_cancer.yaml \
  /data/Genomes/TEST/Data/Final/2026-03-22/ontology_reports/TEST.sqlite \
  /data/Genomes/TEST/Data/Final/2026-03-22/ontology_reports \
  TEST

# From a raw VCF (runs OpenCRAVAT annotation first — hours on WGS):
./run_ontology_report.sh -c config/cardiology.yaml \
  /data/Genomes/TEST/Data/Final/2026-03-22/TEST_master_final.vcf.gz \
  /data/Genomes/TEST/Data/Final/2026-03-22/ontology_reports \
  TEST
```

`-c` defaults to `config/cardiology.yaml`. Input ending in `.sqlite` skips
annotation; anything else is treated as a VCF. Output filenames are namespaced
by domain (`<prefix>_<domain>_report.html`, …) so multiple domains coexist.

---

## Define a new domain in one file

Copy [`config/template.yaml`](config/template.yaml) → `config/<your_domain>.yaml`,
fill in HPO/GO seeds + thresholds, and run with `-c config/<your_domain>.yaml`.
Shipped examples:

| Config | Domain | Panel from |
|--------|--------|-----------|
| `config/cardiology.yaml` | Cardiovascular / cardiometabolic | cardiac HPO + heart/ion-channel GO |
| `config/hereditary_cancer.yaml` | Hereditary cancer predisposition | neoplasm HPO + DNA-repair GO |
| `config/autoimmunity.yaml` | Autoimmune / inflammatory (polygenic) | autoimmune HPO + immune-system GO |
| `config/template.yaml` | (blank scaffold) | — |

Each config can also declare **optional** `domain_predictors` (disease-tuned
models like CardioBoost) and `domain_evidence` columns (e.g. ArrVars,
Cancer-Hotspots, CIViC, GWAS-Catalog). Missing annotator columns are silently
ignored, so a config is portable across databases annotated with different
module sets. A `domain_evidence` entry may set `bypass_frequency: true` to keep
population-**common** variants (e.g. established GWAS risk alleles) that the
usual rare-disease frequency ceiling would otherwise drop.

### Every card now carries gene + genotype context

Two fields were added to **all** domains: the **NCBI Gene description** for each
gene and the **zygosity** (plus variant allele fraction and dbSNP rsID) for each
variant. Zygosity is read straight from the database; the NCBI description is
filled in by an optional **enrichment** stage (below).

### Live enrichment (NCBI Gene + GWAS Catalog)

Stage 5 (`lib/enrich_report.py`) annotates the kept variants on the fly from
public web services, then hands the enriched JSON to the renderer:

- **NCBI Gene** (E-utilities) → per-gene description, summary, cytogenetic band.
- **GWAS Catalog** (EBI REST) → *current* study associations per variant/gene:
  trait, p-value, odds ratio, risk allele, and the PubMed citation.

Enrichment is **cached** (a JSON file next to the outputs), **rate-limited**,
and **offline-safe**: pass `-o` (offline) or run with no network and the report
still builds from cache — it just omits the live layers. Which enrichments run
is set per-config under `report.enrichment`.

### The autoimmunity report (a worked "have fun" domain)

Autoimmune disease is polygenic, so `config/autoimmunity.yaml` turns on the
`gwas_catalog` `bypass_frequency` rule (keep common risk alleles) and selects a
dedicated renderer (`lib/render_autoimmune.py`) that adds an inline **SVG
trait-burden visualization** (how many catalogued GWAS associations point at
each autoimmune trait, coloured by strongest p-value) and per-variant **live
study-evidence tables**. Build a full sample report with no OpenCRAVAT install:

```bash
python3 tests/demo_autoimmune.py /tmp/ai_demo            # live enrichment
python3 tests/demo_autoimmune.py /tmp/ai_demo --offline  # cache-only
# -> /tmp/ai_demo/DEMO_autoimmunity_report.html
```

---

## What it produces (per domain)

| File | Description |
|------|-------------|
| `<prefix>_<domain>_report.html` | **Primary deliverable** — styled, filterable, printable report grouped by tier |
| `<prefix>_<domain>_report.tsv` / `.txt` | Flat table / plain-text summary |
| `<prefix>_<domain>_actionable.xlsx` / `.vcf` | Native OpenCRAVAT exports (full annotation columns) |
| `<prefix>_<domain>_actionable.sqlite` / `.json` | Filtered/tiered variants + structured records |
| `<prefix>_<domain>_panel.json` | The derived gene panel + per-gene ontology rationale |
| `<prefix>_<domain>_enrich_cache.json` | Cached NCBI/GWAS responses (reproducible re-runs, offline mode) |
| `<prefix>_schema.json` | Detected column map for this database |

---

## How selection works (in one paragraph)

The config lists ontology seeds; the panel builder queries the installed `hpo`
and `go` annotator databases and collects every gene annotated with those terms
→ that is the panel. Variants are kept only if they fall in a panel gene **and**
carry a variant-specific signal: ClinVar significance, a damaging coding/splice
consequence, multi-predictor consensus, a configured domain signal, or (for
non-coding) rare regulatory evidence. Everything kept is tiered (Tier 1
reportable / Tier 2 VUS-of-interest / Tier 3 monitor) and tagged with
explainable **reason codes** (`HPO_…`, `GO_…`, `CLINVAR_PLP`, `PP3_CONSENSUS`,
`PVS1_HAPLOINSUFFICIENT`, `SPLICEAI_HIGH`, `PM2_RARE`, …).

See [`docs/PLAN.md`](docs/PLAN.md) for design and [`docs/User_Guide.md`](docs/User_Guide.md)
for usage, tuning, and interpretation.

---

## Repository layout

```
ontology_report/
├── run_ontology_report.sh        # orchestrator: ./run_ontology_report.sh [-c CONFIG] [-o] [-E] <in> <out> <prefix>
├── config/
│   ├── cardiology.yaml           # example domain
│   ├── hereditary_cancer.yaml    # example domain (proves generality)
│   ├── autoimmunity.yaml         # polygenic domain: GWAS bypass + live evidence + viz
│   └── template.yaml             # copy this to define a new domain
├── lib/
│   ├── build_ontology_panel.py   # HPO + GO  -> gene panel
│   ├── schema_probe.py           # detect columns (adds zygosity + rsID)
│   ├── ontology_filter.py        # actionable selection + tiering + reasons (domain-agnostic)
│   ├── enrich_report.py          # NCBI Gene + GWAS Catalog enrichment (cached, offline-safe)
│   ├── make_filtersql.py         # emit --filtersql for native OC reporters
│   ├── render_report.py          # generic HTML / TSV / text renderer
│   └── render_autoimmune.py      # autoimmune renderer (SVG viz + live study cards)
├── tests/
│   ├── test_flow.py              # generic offline mock-DB validation (12 checks)
│   ├── test_autoimmune.py        # autoimmune + zygosity + enrichment validation (offline)
│   └── demo_autoimmune.py        # build a full sample autoimmune report w/o OpenCRAVAT
├── docs/{PLAN,User_Guide,STEP_NOTES}.md
└── logs/                         # dev scratch (gitignored)
```

---

## Requirements

- `micromamba` env **`cravat_env`** (OpenCRAVAT 3.1.1) with modules under
  `/data/opencravat/modules` (`oc module ls`). Set `OC_MODULES_DIR` to override.
- Python: `pyyaml`, `openpyxl` (already in `cravat_env`).

Enrichment additionally uses outbound HTTPS to NCBI E-utilities and the EBI GWAS
Catalog. It is optional: run with `-o`/offline and the report still builds. Set
`NCBI_API_KEY` to raise NCBI rate limits.

## Testing

```bash
python3 tests/test_flow.py        # generic engine: mock SQLite, every decision path
python3 tests/test_autoimmune.py  # autoimmunity + zygosity + enrichment (offline, deterministic)
```

## Scope / disclaimer

Research and personal-screening use. Structural variants (symbolic `<...>` ALT
alleles) are skipped by the OpenCRAVAT VCF converter; this covers SNV/indel and
OpenCRAVAT-resolved deletions. Not a substitute for accredited clinical
diagnostic interpretation.
