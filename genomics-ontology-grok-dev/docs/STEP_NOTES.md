# STEP_NOTES — Per-increment implementation log

Chronological record of the *plan → implement → validate → refactor → validate*
loop, with the concrete evidence gathered at each step.

---

## I0 — Environment & schema discovery

**Findings**
- OpenCRAVAT 3.1.1 in `cravat_env`; modules under `/data/opencravat/modules`.
- `hpo` and `go` are **gene-level** annotators → their columns land in the
  `gene` table, not `variant`. Confirmed:
  - `gene.hpo__id`, `gene.hpo__term`, `gene.hpo__all`
  - `gene.go__bpo_name/_id`, `gene.go__cco_name/_id`, `gene.go__mfo_name/_id`
- Ontology source DBs:
  - `hpo.sqlite → genes(hugo, hpo_id, hpo_term)` (329,339 rows)
  - `go.sqlite → go_annotation(hugo, go_id, go_aspect, qualifier, evidence_code)`
    (400,135 rows) + `go_name(go_id, name)`
- `clingen` is gene-level (`gene.clingen__classification/disease/mondo`);
  `omim__omim_id` and `clinvar_acmg__{ps1,pm5}_id` are variant-level.
- hg38 mapper SO short codes captured (LoF = STG/FSD/FSI/SPL/MLO/STL/EXL/TAB,
  missense = MIS, inframe = IND/INI, complex = CSS).
- 16 cores, 45 GiB RAM.

**Cost control decision**
- Full annotation of the 226k-variant TEST VCF was started in the background
  (finished in ~5 min here because SV lines are skipped — see below). Tooling
  was developed against the pre-existing `Melinda_Ehrle.sqlite` (5.2M rows, same
  annotators) for fast iteration.

---

## I1 — Ontology panel builder (`build_ontology_panel.py`)

- Verified every seed HPO ID resolves to a sensible gene count
  (e.g. HP:0001639 Hypertrophic cardiomyopathy → 290 genes).
- Dropped two mis-labelled/irrelevant seeds (HP:0001636 = Tetralogy of Fallot,
  HP:0011947 = respiratory) and non-matching keywords (`long qt`, `brugada` →
  0 hits; replaced with `prolonged qt`).
- **Validation:** panel = 2,362 genes (hpo_only 1,182 / go_only 944 / both 236 /
  forced 4). Spot-check confirmed all classic cardiac genes present
  (MYH7, MYBPC3, KCNQ1, KCNH2, SCN5A, TTR, LDLR, APOB, PCSK9, RYR2, LMNA, TNNT2,
  DSP, FBN1, PKP2).

---

## I2 — Schema probe + actionable filter (`schema_probe.py`, `ontology_filter.py`)

- Schema probe records the actual column name for each logical field, or `null`
  if absent, so the filter never references a missing column.
- First filter pass kept 317k "actionable" rows — too loose: gene-level HPO/GO
  context alone was letting common intronic VUS through.
- **Refactor:** tightened the actionability gate — a variant now needs a
  *variant-specific* signal; only ClinVar P/LP bypasses the population-frequency
  ceiling. Result on Melinda: **481 → 378** kept (Tier1 7 / Tier2 93 / Tier3 278),
  which is a sensible actionable set. Tier 1 examples verified by hand
  (TMEM43 missense multi-predictor, POLG ClinVar P/LP, NIPSNAP2 SpliceAI-high splice).

---

## I3 — Renderer + orchestrator (`render_report.py`, `run_ontology_report.sh`)

- HTML groups by tier, shows AF/ClinVar/predictors + HPO & GO context + reason
  badges, with a live gene/reason filter and print-to-PDF.
- **Bug found & fixed:** TEST deletions carry multi-kb REF alleles, blowing the
  text report to 145 KB. Added `_fmt_allele()` truncation (`…(Nbp)`); text report
  dropped to ~4 KB. TSV keeps full alleles (it is data).
- Stage 6 emits a `--filtersql` of the actionable UIDs and calls native
  OpenCRAVAT `excel` + `vcf` reporters — verified both files are produced.

---

## I4 — Offline test suite (`tests/test_flow.py`)

- Builds a mock OC-style SQLite (11 hand-crafted variants) covering every path:
  ClinVar P/LP rare → T1, LoF-in-haploinsufficient → T1+PVS1, SpliceAI-high → T1,
  predictor-consensus → T1, VUS → T2, single-predictor+rare → T2, common-benign
  → dropped, non-panel → dropped, no-signal intronic → dropped, ClinVar-P-but-
  common → T3+flag, non-coding double-signal → T3.
- **Validation:** all 11 checks PASS.

---

## I5 — Full validation on TEST database

- TEST annotation: OpenCRAVAT converter skipped **152,834** symbolic-ALT (`<...>`)
  SV lines (`'_SV' object has no attribute 'sequence'`), annotating **69,801**
  SNV/indels + resolved deletions. Documented as a scope limitation.
- End-to-end run: 4,502 panel-gene variants scanned → **19 actionable**
  (Tier2 14 / Tier3 5), dominated by rare exon-loss deletions in cardiac panel
  genes (KCNB2, KCNIP4, LAMA2, TRPM3, …) plus non-coding regulatory hits.
- HTML rendered and visually verified in-browser.
- Native Excel (118 KB) + VCF (515 KB) exports confirmed.

---

## I6 — Refactor & re-validate

- Consolidated gene-level column handling (hpo/go/clingen) behind
  `GENE_TABLE_KEYS`; schema-aware `build_select()` substitutes `NULL` for absent
  columns so the same code runs on databases with different annotator sets
  (validated on both Melinda — no clingen/omim — and TEST — with clingen/omim).
- Re-ran `tests/test_flow.py` (green) and the TEST orchestrator (stable output).

## I7 — Generalization to any domain (de-cardio)

**Motivation:** the report must not be limited to cardiology or a fixed ontology
focus.

**Changes**
- Renamed project `cardio_ontology_report` → `ontology_report`.
- Config is now per-domain: `config/cardiology.yaml`, `config/hereditary_cancer.yaml`
  (second worked domain), `config/template.yaml` (blank scaffold). Added
  `report_title`, `domain_predictors`, `domain_evidence` sections.
- `ontology_filter.py`: removed hard-coded CardioBoost + ArrVars. They are now
  generic, **config-driven** mechanisms:
  - `domain_predictors` — arbitrary disease-tuned predictor columns + threshold →
    `PP3_<code>` + consensus contribution.
  - `domain_evidence` — arbitrary columns whose presence adds a reason and makes a
    variant actionable.
  Both resolve declared columns against the actual DB at runtime and ignore
  absent ones, so one config is portable across annotator sets.
- Renamed output column `cardio_tier` → `tier` (sqlite/json/renderer/tests).
- `schema_probe.py`: dropped cardio-specific fields (cardioboost/arrvars).
- `render_report.py`: title/header/footer now come from config `report_title`
  and `domain`.
- `run_ontology_report.sh`: added `-c CONFIG`; deliverables namespaced by domain.
- `tests/test_flow.py`: updated to `tier`; added a 12th check exercising the
  config-driven `domain_evidence` (ArrVars) path.

**Validation**
- Offline suite: **12/12 pass**.
- TEST DB, one database → two domains:
  - `cardiology` — panel 2,362 genes, **19 actionable** (active domain signals:
    CardioBoost, ArrVars).
  - `hereditary_cancer` — panel 2,303 genes (BRCA1/BRCA2/MLH1/MSH2/APC/TP53/PALB2/…
    present), **8 actionable**, cancer-specific reasons (`HPO_CARCINOMA`,
    `GO_DOUBLE-STRAND_BREAK`, `GO_REGULATION_OF_APOPTOTIC_PROCESS`).
  This is the key proof the engine is not cardiology-bound.

## I8 — Gene/genotype context + live-evidence autoimmunity report

**Motivation:** two requested field additions (NCBI Gene description per gene,
zygosity per variant) plus a creative new domain — autoimmunity — that exercises
"live" reporting (current study results, a visualization, loading extra modules).

**Changes**
- `schema_probe.py`: detect `sample`-table genotype, `vcfinfo__zygosity`/reads,
  and dbSNP rsID (`dbsnp__rsid`). Reports sample-column count.
- `ontology_filter.py`:
  - `zygosity_label()` normalizes het/hom/hemi (and 0/1, 1|1, … VCF spellings);
    `compute_vaf()` derives allele fraction from `alt_reads/tot_reads` when no
    explicit VAF. Both surfaced in `evidence`; rsID pulled through.
  - Multi-sample DBs: when the variant table lacks genotype, a `sample`-table
    uid→genotype map is built and merged per row.
  - `domain_evidence` gained `bypass_frequency`: a firing entry keeps a variant
    over the frequency ceiling (reason `RISK_ALLELE_COMMON` instead of
    `COMMON_AF_FLAG`). Needed because autoimmune risk alleles are common.
- `enrich_report.py` (new): cached, rate-limited, offline-safe HTTP-JSON client.
  - NCBI Gene: esearch→esummary → description/summary/map_location.
  - GWAS Catalog: per-rsID associations (trait, p-value reconstructed from
    mantissa/exponent so tiny p-values survive, OR/β, risk allele) + a follow to
    the study's PubMed citation; plus a cheap gene→rsID rollup. Cache stores both
    `ok` and `err` outcomes so failures don't retry within a run.
- `render_report.py`: cards/TSV/text now show zygosity, VAF, dbSNP, and (when
  enriched) the NCBI gene description. Degrades gracefully when absent.
- `render_autoimmune.py` (new): autoimmune renderer — inline SVG **trait-burden**
  chart (associations per trait, coloured by −log10 p) + per-variant GWAS
  study-evidence tables + NCBI gene summary. Reuses generic helpers/TSV/text.
- `config/autoimmunity.yaml` (new): autoimmune HPO + immune GO seeds; force-include
  HLA/PTPN22/CTLA4/IL23R/STAT4/TYK2/AIRE/FOXP3/…; `gwas_catalog` bypass evidence;
  `report:` block (renderer=autoimmune, enrichment genes+studies, visualization).
- `run_ontology_report.sh`: 7 stages (adds enrichment 5/7); reads `report.*` from
  config to pick renderer + enrichment; `-o` offline / `-E` no-enrich flags;
  added `dbsnp gwas_catalog vcfinfo` to the annotator set.
- Tests: `tests/test_autoimmune.py` (offline, deterministic, pre-seeded cache),
  `tests/demo_autoimmune.py` (buildable sample report without OpenCRAVAT).

**Validation**
- `tests/test_flow.py` — 12/12 (no regression).
- `tests/test_autoimmune.py` — all checks pass: schema detection, GWAS-common
  bypass, zygosity het/hom/hemi, VAF from depths, monogenic Tier1/PVS1, negatives
  dropped, offline enrichment injection (0 network calls), SVG + study tables.
- Online demo (`tests/demo_autoimmune.py`): 12/12 records enriched — real NCBI
  descriptions and current GWAS associations (PTPN22/rs2476601 → type 1 diabetes
  p≈2e-80, RA, hypothyroidism, …); 1 benign 404 for an rsID absent from the GWAS
  catalog, handled gracefully. Report HTML ~60 KB with a 14-bar trait chart and 8
  study tables. Offline re-run reproduces from cache with 0 calls.

## Open follow-ups (nice-to-have)
- HPO term-hierarchy expansion (descendant terms) instead of keyword/ID seeds.
- Optional SV-aware converter so large deletions get population AF.
- Per-gene collapsible summary section in the HTML.
- Map raw GWAS traits onto EFO parents to de-noise the trait chart (drop pure
  measurement traits like "platelet count").
- Weight the autoimmune tiering by GWAS effect size / p-value, not just the
  monogenic signals.
