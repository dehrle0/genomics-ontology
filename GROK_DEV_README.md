# grok-dev Branch – Visual Ontology Explorer

This branch adds the **Visual Ontology Explorer** (Option 3) interactive report UI.

## What’s new

- `render_visual_ontology_explorer.py` – full self-contained HTML generator
  - Left panel: hierarchical HPO ontology tree (organ system → phenotype → gene)
  - Right panel: Gene Details with tabs **Overview · Variants · Phenotypes · Publications**
  - Explicit **maternal / paternal / undetermined** phasing badges
  - ~70 % of heterozygous variants in the demo data are phased
- `samples/demo_variant_report.json` – realistic sample that matches the `VariantReport` Pydantic / LinkML schema (fields already present in OpenCRAVAT pipeline output)
- Phasing field was already defined in `genomics_ontology_io/models.py` (`phasing: maternal | paternal | de_novo | undetermined`)

## Quick start (local test)

```bash
# From repo root
python render_visual_ontology_explorer.py --demo -o reports/visual_ontology_explorer.html

# Or feed a real report JSON produced by the existing pipeline
python render_visual_ontology_explorer.py -i path/to/your_variant_report.json -o reports/my_report.html

# Serve and open
python -m http.server 8080 --directory reports
# → http://localhost:8080/visual_ontology_explorer.html
```

## Data expectations (OpenCRAVAT / existing pipeline)

The renderer expects a JSON object that validates against `VariantReport`:

| Field | Source |
|-------|--------|
| `monogenic_findings[].gene_symbol` | OpenCRAVAT `hugo` / gene annotation |
| `monogenic_findings[].rsid` | dbSNP |
| `monogenic_findings[].chromosome`, `.position` | VCF |
| `monogenic_findings[].zygosity` | genotype call |
| `monogenic_findings[].phasing` | WhatsHap / SHAPEIT / long-range LD (already modelled) |
| `monogenic_findings[].clinvar_significance` | ClinVar |
| `monogenic_findings[].revel_score` | REVEL |
| `monogenic_findings[].impact_consequence` | VEP / OpenCRAVAT SO terms |
| `monogenic_findings[].associated_hpo_terms` | HPO gene/variant annotations |
| `polygenic_findings[]` | PGS Catalog roll-ups |
| `pharma_findings[]` | PharmCAT / CPIC |

## How to create the branch on GitHub

```bash
git clone https://github.com/dehrle0/genomics-ontology.git
cd genomics-ontology
git checkout dev
git checkout -b grok-dev

# Copy the new files from this package into the repo:
#   render_visual_ontology_explorer.py
#   samples/demo_variant_report.json
#   GROK_DEV_README.md
# (and any other files you want)

git add render_visual_ontology_explorer.py samples/ GROK_DEV_README.md
git commit -m "feat: Visual Ontology Explorer report UI (Option 3) with phasing"
git push -u origin grok-dev
```

Then open a PR against `dev` when ready.

## Notes

- The ontology tree is built dynamically from the HPO terms present on the variants; richer phenotype labels can later be injected from a local HPO OBO cache.
- Publications are currently stubbed for the demo genes (TTN, LMNA, MYH7). Production can pull from Open Targets / PubMed / EuropePMC caches already used by the enrichment step.
- The original `render_new_ontology_report.py` is left untouched so both UIs can coexist.
