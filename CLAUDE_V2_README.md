# Genomic Ontology Explorer — prototype

A frontend prototype of the Gene Inspector / GIP-style variant explorer,
built from the Facebook mockup screenshot, the Grok mockup spec, and the
dehrle0/genomics-ontology reference. **No build step** — unzip and open
`index.html` in a browser.

```
genomic-ontology-claude/
├── index.html          shell: top nav + 5 views
├── css/style.css        design tokens + all component styles
├── js/app.js             all interaction logic (vanilla JS, no framework)
├── data/mock-data.js     sample dataset (see "Data" below)
└── README.md             this file
```

## Run it

Just open `index.html` directly, or serve it locally if your browser
restricts local `<script src>` loading:

```bash
cd genomic-ontology-claude
python3 -m http.server 8000
# then open http://localhost:8000
```

## What's implemented

- **Ontology view** — left tree, selectable between **HPO**, **GO**, and
  **Organ/System**, each as Organ/category → term → gene. Tree/List
  layout toggle, "only systems with findings" scope toggle, and a
  search box that filters + auto-expands matches. Clicking a gene loads
  the right panel.
- **Gene detail (right panel)** — four tabs: **Overview** (KPI strip,
  gene summary, associated pathology), **Phenotypes** (HPO term cards),
  **Variants** (grouped into Potential Concerns / Protective
  Associations / Uncertain Findings / Not Categorized, colored exactly
  like the mockup screenshot, each variant expandable to show the
  "what studies say / condition / genotype relevance" rows), and
  **Publications** (curated cards with DOI links). "View in Genome
  Browser" opens a stub modal — the intended hook for an IGV.js /
  locus-view embed.
- **Genes view** — one row per gene (variant count, pathogenic/LP
  count, phased-het %, HPO term count), sortable, filterable by
  symbol, click-through to Ontology detail.
- **Variants view** — full cross-gene table with filters for gene,
  ClinVar class, phase, and zygosity; click a row to jump to that
  gene's Variants tab.
- **Analysis view** — KPI strip, Polygenic Risk Score cards (trait,
  percentile bar, HIGH/MODERATE/PROTECTIVE, PGS Catalog ID), and a
  Pharmacogenomics table (gene, diplotype, phenotype, drug, action
  tier, CPIC/DPWG text).
- **Reports view** — narrative summary, gene breakdown table, Print
  (browser print-to-PDF), and Download JSON (the same payload shape
  the UI consumes).

## Assumptions made (per your instruction to note all of them)

1. **No live data.** The 7.5 GB `.sqlite` output lives on your own
   server and wasn't uploaded; I have no network access from this
   session either. Everything here runs on `data/mock-data.js`, built
   to mirror your real fields (rsID, genotype, zygosity, phase, MAF,
   consequence, ClinVar, REVEL, HPO terms, PGS/PGx) so swapping in a
   real API is a data-layer change, not a UI rewrite (see "API seam"
   below). Genes VDR/MTHFR/AGXT2/CBS come straight from your
   screenshot; I added BRCA1/CFTR/LDLR/G6PD so the ontology trees and
   organ-system groupings had enough breadth to be worth navigating.
2. **"Verify with a different model, up to 3 times."** I can't invoke
   a second model instance inside this session, so I substituted three
   self-review passes instead of a second opinion — flagging this
   plainly rather than implying an independent review happened:
   - **Pass 1 — wiring**: every `getElementById` call cross-checked
     against the HTML (one flagged case, `genome-browser-btn`, is
     created dynamically inside the gene-detail render and is fine);
     confirmed every data field the UI reads exists in the mock
     dataset with matching names; verified `<div>` tag balance (50/50).
   - **Pass 2 — logic/consistency**: traced re-render and event-binding
     flows (variant row expand/collapse, cross-view navigation from
     Genes/Variants tables back into the Ontology detail panel) for
     stale-listener or state-mismatch bugs; checked CSS class name
     agreement between JS-generated classes (`cat-pill--*`,
     `badge--*`, `variant-section__head.*`, `action-tier.*`) and the
     stylesheet; the two flagged "duplicate selector" hits were both
     legitimate responsive/state overrides, not conflicts.
   - **Pass 3 — fidelity to spec**: re-read the Grok mockup notes and
     your screenshot against the built views to confirm each named
     panel (Ontology, Genes, Variants, Analysis, Reports; Overview /
     Phenotypes / Variants / Publications gene tabs) exists and that
     the category colors/labels match the screenshot's red
     (Potential Concerns) / green (Protective Associations) / blue
     (Uncertain Findings) / gray (Not Categorized) scheme.
   I did **not** run this in an actual browser (no network to fetch a
   headless-browser dependency in this sandbox) — I'd treat a real
   browser smoke-test as the next verification step before you rely
   on this beyond a prototype.
3. **Single-file-per-concern, no framework.** Given no backend and a
   need for a zero-install deliverable, I used plain HTML/CSS/JS
   rather than React — trivial to port later if you want the actual
   product on a component framework.
4. **Visual direction** — I matched the clinical/dense style of your
   real screenshot (white panels, colored severity badges, monospace
   for every genomic identifier) rather than a marketing-style
   redesign, since the goal was fidelity to existing mockups, not a
   new brand.
5. **Unused-but-present data** — `researchedVariants` and
   `evidenceLevel` are in the mock schema (they're in your real
   OpenCRAVAT-adjacent data) but not yet surfaced in the UI. Flagging
   rather than silently dropping them — natural next additions (see
   Roadmap).

## API seam (for wiring up the real backend)

`data/mock-data.js` is the entire seam. Replace the four top-level
consts (`GENES`, `ONTOLOGIES`, `PRS`/`PGX`, `REPORT`) with `fetch()`
calls to:

```
GET /api/genes                → GENES (list)
GET /api/genes/:symbol        → one GENES entry, full detail
GET /api/ontology/:type       → ONTOLOGIES[type]   (type = hpo | go | organ)
GET /api/analysis             → { prs: PRS, pgx: PGX }
GET /api/report               → REPORT
```

`js/app.js` never reaches into the data layer except through these
globals, so no other file needs to change.

## Revision 2 — what changed

Based on your feedback:

- **List is now the default** left-panel layout (Tree and a new **Graph**
  mode are still there as toggles). Graph draws an actual node/edge
  diagram (Group → Term, SVG lines and circles) with gene chips
  attached at each term row — a genuine third option, not a relabeled
  list.
- **Organ/System now has a real second layer** (System → Sub-system →
  Gene), matching HPO/GO's shape, so Tree and List actually look
  different for it too — previously Organ/System had no term layer and
  the two modes were identical.
- **Filled ontology gaps**: added an HPO "Immune system" group (VDR ↔
  multiple sclerosis), two more GO terms (BRCA1 double-strand-break
  repair, VDR vitamin-D response, LDLR lipoprotein clearance), and the
  Organ/System sub-system layer above. Every gene is now confirmed (see
  verification below) to appear in all three trees with no orphaned
  nodes.
- **Gene tab list**: added **Studies** between Variants and
  Publications — every individual study finding across a gene's
  variants, one card each, with evidence-strength dots and a
  click-through back to that variant.
- **Overview tab**: expanded the KPI row (added Protective, Uncertain,
  GO terms, Publications counts alongside the existing ones) and added
  a **Reference links** row — NCBI Gene, OMIM (gene + phenotype
  entries, separately), GeneCards, ClinVar gene search.
- **Variants tables** (both the per-gene tabs and the global Variants
  view) gained **CADD, SpliceAI, AlphaMissense, QUAL (call quality),
  and read support (alt/total reads, %)** columns, plus a quick
  NCBI Gene / OMIM link pair inside each expanded variant's detail row.
  Tables now scroll horizontally rather than compressing columns.
- **Faint colorful DNA watermark**: a generated double-helix (colored
  by base — A/T/G/C — using the app's own palette) sits fixed behind
  the page at 7% opacity. It's placed so it only shows through the
  app's transparent gaps (padding, empty states) and is always covered
  by the white cards and tables — it can't reduce data legibility by
  construction, not just by low opacity.
- **Accuracy fixes** (see the CHANGELOG comment at the bottom of
  `data/mock-data.js` for detail): AGXT2 was incorrectly linked to
  "Hyperoxaluria, primary, type III" — that phenotype is caused by
  HOGA1, not AGXT2; corrected. OMIM numbers were previously mixing
  gene-entry and phenotype-entry numbers under one ambiguous field for
  several genes (CFTR, LDLR, CBS, BRCA1, MTHFR, VDR) — split into
  `omimGene` / `omimPhenotype` and corrected each. Added real NCBI
  Entrez Gene IDs for the new reference-links row.

### Verification passes performed (4, per your instruction)

Same caveat as last time: I can't invoke a second model instance in
this session, so these are self-review passes against the spec and
against the data itself, not an independent second opinion — flagging
that plainly again rather than overstating it.

1. **Syntax + structure** — both JS files parse cleanly
   (`node -c`); `index.html`'s `<div>` tags balance (50 open / 50
   close).
2. **DOM wiring** — every `getElementById` call cross-checked against
   the HTML or confirmed as dynamically created before use
   (`genome-browser-btn`, `graph-host`); ontology/gene cross-reference
   integrity checked programmatically: every gene referenced in every
   tree exists in `GENES`, every gene is placed in all three trees, and
   no gene sits in a group without also sitting under one of that
   group's terms (which would silently vanish in Tree mode).
3. **New data fields** — every gene has the fields the Overview links
   row reads (`ncbiGeneId`, `omimGene`, all four `links`); every
   variant has the five new score fields, with sanity-checked ranges
   (CADD 0–60, AlphaMissense/SpliceAI 0–1, read counts where
   matching ≤ total).
4. **CSS/JS class agreement + table column integrity** — every new
   class string built in `app.js` (graph nodes/edges, study cards,
   score cells, reference links) has a matching rule in
   `style.css`; column counts were hand-traced against generated row
   cells for all three tables (gene-level variant table: 13 header
   cells vs. 13 generated cells, including the expandable row's
   `colspan="13"`; global Variants table: 14 vs. 14; Genes table: 7
   vs. 7).

**Not done**: an actual in-browser render. This sandbox has no
headless browser available and no network access to install one, so
the graph layout's pixel-level spacing (chip positions at `GENE_X`
increments) is verified by math/trace, not by looking at it. With the
gene counts currently in the mock data (max 2 genes per term) it fits
the fixed graph width comfortably; if you extend the dataset to genes
with many more per term, that width may need to become dynamic rather
than fixed — flagging now rather than waiting for it to surface as a
bug. Opening it in a real browser is the natural next check before you
lean on this further.

## Roadmap (flagged, not built)

- Real **IGV.js** (or exon/domain diagram) embed behind "View in
  Genome Browser," replacing the current stub modal.
- Surface `researchedVariants` and `evidenceLevel` in the Overview and
  Phenotypes tabs (evidence-strength badges on HPO cards).
- Predictor-score columns (AlphaMissense, CADD, SpliceAI, BayesDel,
  MetaRNN) as an expandable row under each variant, and an optional
  "predictor consensus" rollup in Analysis.
- Gene-level constraint (pLI/LOEUF) is already in the mock schema and
  shown in the gene header; wire to gnomAD constraint data for real
  samples.
