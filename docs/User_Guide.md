# User Guide — ontology_report

This guide covers running the report, reading the output, and tuning the
ontology selection. For the design rationale see [`PLAN.md`](PLAN.md).

---

## 1. Prerequisites

```bash
micromamba activate cravat_env
oc module ls | grep -E "hpo|go|clinvar|clingen|omim|gnomad4"   # confirm modules present
```

The report reads two ontology databases directly:

- `/data/opencravat/modules/annotators/hpo/data/hpo.sqlite`
- `/data/opencravat/modules/annotators/go/data/go.sqlite`

If your modules live elsewhere, set `OC_MODULES_DIR` before running:

```bash
export OC_MODULES_DIR=/path/to/opencravat/modules
```

---

## 2. Running

```bash
cd ~/My-Projects/genomics/development/ontology_report
./run_ontology_report.sh [-c CONFIG] [-o] [-E] <INPUT> <OUTPUT_DIR> <PREFIX>
```

- `-c CONFIG` — the **domain** config (default `config/cardiology.yaml`). This is
  what makes the report cardiology, cancer, autoimmunity, neuro, etc.
- `-o` — **offline** enrichment: never hit the network; use the enrichment cache
  only (the report still builds, minus any un-cached live layers).
- `-E` — skip the enrichment stage entirely.
- `<INPUT>` — a `.vcf`/`.vcf.gz` (annotation runs first, hours on WGS) **or** an
  existing OpenCRAVAT `.sqlite` (annotation skipped, seconds to minutes).
- `<OUTPUT_DIR>` — where all deliverables are written.
- `<PREFIX>` — patient/sample label used in filenames and the report header.

Output files are namespaced by domain (`<prefix>_<domain>_report.html`, …) so you
can run several domains into the same directory without collisions.

### Example (two domains, same data)

```bash
OUT=/data/Genomes/TEST/Data/Final/2026-03-22/ontology_reports
./run_ontology_report.sh -c config/cardiology.yaml        "$OUT/TEST.sqlite" "$OUT" TEST
./run_ontology_report.sh -c config/hereditary_cancer.yaml "$OUT/TEST.sqlite" "$OUT" TEST
```

### The seven stages

1. **Annotate** (skipped for `.sqlite`) — `oc run` with the full annotator set.
2. **Build panel** — HPO + GO → `<prefix>_<domain>_panel.json`.
3. **Probe schema** — detect which annotator columns exist (incl. genotype /
   zygosity and dbSNP rsID) → `<prefix>_schema.json`.
4. **Filter** — actionable selection + tiering → `<prefix>_<domain>_actionable.{sqlite,json}`.
5. **Enrich** — NCBI Gene descriptions and (if the config asks) live GWAS Catalog
   study evidence; cached to `<prefix>_<domain>_enrich_cache.json`. Skipped by
   `-E`; `-o` forces cache-only.
6. **Render** — HTML + TSV + text (renderer chosen by the config's
   `report.renderer`: `generic` or `autoimmune`).
7. **Native export** — OpenCRAVAT Excel + VCF of exactly the actionable set.

---

## 3. Reading the report

Open `<prefix>_<domain>_report.html` in a browser. Layout:

- **Header stats** — actionable count, per-tier counts, panel size, variants scanned.
- **Filter box** — type a gene symbol or a reason code to live-filter cards.
- **Tier sections**
  - **Tier 1 — Reportable / Pathogenic-grade**
  - **Tier 2 — VUS of interest**
  - **Tier 3 — Monitor / regulatory**
- **Each card** shows location, consequence, **zygosity** (+ variant allele
  fraction and dbSNP rsID), population AF (gnomAD4 + All of Us), ClinVar (linked),
  key predictor scores, the **NCBI Gene description** (after enrichment), the
  gene's **HPO phenotype context** and **GO function context** (with outbound
  HPO/GenCC links), and the full list of **reason codes**.
- **Print / Export PDF** button uses the browser print dialog.

### Zygosity

Zygosity is normalized to `Heterozygous` / `Homozygous` / `Hemizygous` from
whatever the database carries (single-sample `vcfinfo__zygosity`, or the
per-sample `sample` table on multi-sample databases). When only read depths are
present the variant allele fraction is derived as `alt_reads / tot_reads`.

### NCBI Gene description

Added by the enrichment stage from NCBI Gene (E-utilities): the official gene
description and cytogenetic band, with a link to the gene's NCBI page. Absent if
enrichment was skipped or the gene could not be resolved.

### The autoimmunity report

`config/autoimmunity.yaml` uses the `autoimmune` renderer, which adds two things
on top of the standard tiered cards:

- a **trait-burden chart** (inline SVG) summarising how many catalogued GWAS
  risk-allele associations point at each autoimmune trait across your genome,
  coloured by the strongest reported p-value; and
- per-variant **current GWAS Catalog study evidence** (trait, p-value, odds
  ratio, risk allele, PubMed link), fetched live and cached.

Because autoimmune risk is polygenic, this domain keeps **common** risk alleles
via the `gwas_catalog` `domain_evidence` rule (`bypass_frequency: true`); such
variants are tagged `GWAS_RISK_ALLELE` and `RISK_ALLELE_COMMON`.

### Reason codes

| Code | Meaning | Kind |
|------|---------|------|
| `CLINVAR_PLP` | ClinVar Pathogenic / Likely pathogenic | phenotype |
| `CLINVAR_VUS` | ClinVar Uncertain significance | phenotype |
| `CLINVAR_CONFLICT` | ClinVar conflicting classifications | phenotype |
| `CLINGEN_VALIDITY` | ClinGen gene–disease validity present | phenotype |
| `OMIM_DISEASE` | OMIM disease entry present | phenotype |
| `HPO_<TERM>` | Gene's HPO annotation matched a configured phenotype term | phenotype |
| `GO_<TERM>` | Gene's GO annotation matched a configured function term | phenotype |
| `<CODE>` (from `domain_evidence`) | Configured domain evidence column present (e.g. `ARRVARS_KNOWN`, `CANCER_HOTSPOT`) | phenotype/genotype |
| `LOF_<so>` | Loss-of-function consequence (stop/frameshift/splice/…) | genotype |
| `MISSENSE` / `INFRAME_INDEL` | coding consequence | genotype |
| `PP3_<predictor>` | Predictor above threshold (REVEL/AM/… + configured `domain_predictors`) | genotype |
| `PP3_CONSENSUS` | ≥ N predictors agree (config `pp3_consensus_n`) | genotype |
| `SPLICEAI_HIGH` / `SPLICEAI_MOD` | SpliceAI Δ ≥ tier1 / ≥ min | genotype |
| `PM2_RARE` / `RARE_TIER2` | rare at Tier 1 / Tier 2 frequency threshold | genotype |
| `NONCODING_CADD/LINSIGHT/NCER/REGULOMEDB` | regulatory support (non-coding) | genotype |
| `PVS1_HAPLOINSUFFICIENT` | LoF in a curated haploinsufficient gene | genotype |
| `GWAS_RISK_ALLELE` | catalogued in `gwas_catalog` (config `domain_evidence`) | phenotype |
| `COMMON_AF_FLAG` | kept because ClinVar P/LP but common in population | caution |
| `RISK_ALLELE_COMMON` | kept via a `bypass_frequency` domain evidence though common (e.g. GWAS risk allele) | caution |

---

## 4. Defining / tuning a domain — the only file you edit

Pick a domain config under `config/` (or copy `config/template.yaml` to make a
new one). Everything below lives in that file.

### Add/remove phenotypes (HPO)
Add exact term IDs (high precision) and/or keyword substrings (broader recall):

```yaml
hpo:
  term_ids:
    - HP:0001644   # Dilated cardiomyopathy
  term_keywords:
    - cardiomyopathy
```

Find term IDs with:

```bash
python3 - <<'PY'
import sqlite3
c=sqlite3.connect("/data/opencravat/modules/annotators/hpo/data/hpo.sqlite").cursor()
for r in c.execute("select distinct hpo_id, hpo_term from genes where lower(hpo_term) like '%arrhythmia%' limit 20"):
    print(r)
PY
```

### Add/remove functions (GO)
```yaml
go:
  term_keywords:
    - heart contraction
    - potassium channel activity
```

### Precision vs recall
```yaml
panel:
  min_ontology_support: 1   # 1 = union (recall). 2 = gene must hit BOTH HPO and GO (precision).
```

### Thresholds
- `frequency:` — `max_af` is the hard actionability ceiling; `tier1_af`/`tier2_af`
  drive rarity-based tier upgrades.
- `predictors:` — generic pan-disease cutoffs and `pp3_consensus_n`.
- `noncoding:` — CADD/LINSIGHT/ncER/RegulomeDB cutoffs.
- `haploinsufficient_genes:` — genes where LoF → Tier 1.

### Optional domain-specific signals (portable, no code changes)
- `domain_predictors:` — disease-tuned model columns + a threshold. Each firing
  adds a `PP3_<code>` reason and counts toward the predictor consensus.
  ```yaml
  domain_predictors:
    - code: CARDIOBOOST
      score_cols: [cardioboost__cardiomyopathy, cardioboost__arrhythmias]
      text_cols:  [cardioboost__cardiomyopathy1, cardioboost__arrhythmias1]
      min: 0.90
  ```
- `domain_evidence:` — any OpenCRAVAT column(s) whose presence makes a variant
  actionable and adds a reason code:
  ```yaml
  domain_evidence:
    - code: CANCER_HOTSPOT
      kind: phenotype
      columns: [cancer_hotspots__samples]
    - code: GWAS_RISK_ALLELE
      kind: phenotype
      bypass_frequency: true          # keep this variant even if population-common
      columns: [gwas_catalog__trait, gwas_catalog__pmid]
  ```
  Columns that are not present in a given database are ignored automatically, so
  the same config runs against databases with different annotator sets.
  `bypass_frequency: true` (default false) lets an evidence hit keep a variant
  that exceeds the frequency ceiling — essential for polygenic domains whose risk
  alleles are common (see `config/autoimmunity.yaml`).

### Enrichment (config `report:` block)
```yaml
report:
  renderer: autoimmune       # or omit / "generic"
  enrichment:
    genes: true              # NCBI Gene descriptions (default on)
    studies: true            # live GWAS Catalog study evidence
```
The orchestrator reads this block to decide which enrichments to run and which
renderer to use. Run with `-o` for offline (cache-only) or `-E` to skip.

After editing, just re-run stages 2–6 (they are fast); re-annotation is only
needed if you add a **new annotator** that must be present in the SQLite.

---

## 5. Rebuilding just the panel (to preview gene selection)

```bash
python3 lib/build_ontology_panel.py --config config/hereditary_cancer.yaml --out /tmp/panel.json
python3 - <<'PY'
import json; p=json.load(open("/tmp/panel.json"))
print(p["counts"])
print("BRCA1 in panel:", "BRCA1" in p["genes"])
PY
```

---

## 6. Troubleshooting

| Symptom | Cause / fix |
|---------|-------------|
| `HPO data not found` | Wrong modules dir — set `OC_MODULES_DIR`. |
| Panel is huge/tiny | Adjust `min_ontology_support` or trim keywords in the config. |
| Many `COMMON_AF_FLAG` Tier 1s | These are ClinVar P/LP but population-common; review manually. |
| Native Excel/VCF missing | Non-fatal; the orchestrator warns and continues. Re-run stage 7 manually with `oc report`. |
| Few actionable variants on test data | The TEST VCF is SV-heavy; symbolic `<...>` ALTs are dropped by the OC converter (see PLAN §Scope). |
| No gene descriptions / study evidence | Enrichment was skipped (`-E`), offline (`-o`) with an empty cache, or the network was unreachable. The report still builds; re-run online to populate. |
| Zygosity shows `-` | The database has no genotype columns (`vcfinfo__zygosity` or a `sample` table). Annotate with `vcfinfo`, or use a multi-sample DB. |
| Enrichment slow | Each new gene/variant is one or more API calls (cached afterwards). Set `NCBI_API_KEY`; results cache to `<prefix>_<domain>_enrich_cache.json`. |

---

## 7. Validating changes

Always run the offline suites after editing `lib/`:

```bash
python3 tests/test_flow.py        # generic engine, every decision path
python3 tests/test_autoimmune.py  # autoimmunity, zygosity, enrichment, renderer
```

`test_flow.py` builds a mock database exercising every decision path and asserts
the expected tiers and reason codes (including the config-driven
`domain_evidence` path), so it is domain-independent. `test_autoimmune.py`
additionally checks zygosity normalization, the `bypass_frequency` mechanism,
offline enrichment injection (from a pre-seeded cache — no network), and the
autoimmune renderer's SVG chart + study tables. To eyeball a full sample report:

```bash
python3 tests/demo_autoimmune.py /tmp/ai_demo            # live enrichment
python3 tests/demo_autoimmune.py /tmp/ai_demo --offline  # cache-only
```

---

## 8. Defining a brand-new domain (worked recipe)

1. `cp config/template.yaml config/neuro.yaml`
2. Set `domain: neuro`, `report_title: …`.
3. Find relevant HPO IDs/keywords (`Seizure`, `Intellectual disability`, …) and
   GO terms (`neuron differentiation`, `synaptic signaling`, …) using the
   snippet in §4; paste them into `hpo:`/`go:`.
4. Optionally add `force_include`, `haploinsufficient_genes`, `domain_predictors`,
   `domain_evidence`.
5. Preview the panel (§5). Adjust `min_ontology_support`/keywords until the gene
   count and known genes look right.
6. Run: `./run_ontology_report.sh -c config/neuro.yaml <db> <outdir> <prefix>`.

No Python changes are ever required to add a domain.
