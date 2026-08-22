# grok-dev – Visual Ontology Explorer (Option 3)

High-fidelity interactive HTML report UI for gene / variant / phenotype exploration.

## What's new in this upgrade

- Top navigation bar (Ontology · Genes · Variants · Analysis · Reports)
- Left hierarchical HPO ontology graph (organ system → phenotype → gene)
- Right Gene Details panel with tabs:
  - **Gene Overview** – counts + NCBI summary
  - **Phenotypes** – HPO term cards
  - **Variants** – full table (HGVS/rsID, ClinVar badge, REVEL bar, coordinate, zygosity, **Phase**, allele freq, last evaluated)
  - **Publications** – curated study cards with DOI links
- Explicit **Maternal / Paternal / Unknown** phase badges
- ~70 % of heterozygous variants phased in the demo data
- Works with or without `pydantic` installed (falls back to raw dict)

---

## Quick demo (no real data needed)

```bash
cd ~/My-Projects/genomics-ontology/genomics-ontology   # your repo root

# Use your existing venv (or create one)
source .venv/bin/activate          # if not already active
# pip install pydantic             # only needed for strict validation

python render_visual_ontology_explorer.py --demo -o reports/visual_ontology_explorer.html
python -m http.server 8081 --directory reports
# → http://localhost:8081/visual_ontology_explorer.html
```

---

## Using real data from the existing pipeline

### Input contract

The renderer expects a single JSON file that matches the `VariantReport` schema
(`genomic_ontology_schema.yaml` / `genomics_ontology_io.models.VariantReport`):

```json
{
  "patient_id": "…",
  "run_date": "…",
  "monogenic_findings": [
    {
      "gene_symbol": "TTN",
      "ncbi_description": "…",
      "rsid": "rs…",
      "chromosome": "chr2",
      "position": 179431234,
      "genotype": "C/T",
      "zygosity": "Heterozygous",
      "revel_score": 0.91,
      "impact_consequence": "Nonsense",
      "clinvar_significance": "Pathogenic",
      "phasing": "maternal",
      "associated_hpo_terms": ["HP:0001626", "HP:0001644"],
      "associated_mondo_terms": [],
      "gnomad_af": 1.2e-5,
      "last_evaluated": "2024-02-12"
    }
  ],
  "polygenic_findings": [],
  "pharma_findings": []
}
```

### End-to-end with OpenCRAVAT (already in this repo)

```bash
# 1. Annotate VCF / gVCF with OpenCRAVAT
#    Required annotators (already listed in run_ontology_master_pipeline.sh):
#      hpo  go  clinvar  clingen  omim  ncbigene  revel  alphamissense
#      bayesdel  metarnn  esm1b  varity  spliceai  cadd  dbsnp  vcfinfo
#      gwas_catalog  pharmgkb  …

./run_ontology_master_pipeline.sh \
    /path/to/sample.vcf.gz \
    /path/to/output_dir \
    SAMPLE_NAME

# Produces:
#   SAMPLE_NAME_master_actionable.json
#   SAMPLE_NAME.sqlite
#   SAMPLE_NAME_master_actionable.sqlite

# 2. (Optional) Enrich with NCBI Gene descriptions + GWAS studies
python lib/enrich_report.py \
    --genes --studies \
    -i /path/to/output_dir/SAMPLE_NAME_master_actionable.json \
    -o /path/to/output_dir/SAMPLE_NAME_enriched.json

# 3. Thin adapter: actionable JSON → VariantReport shape
python - <<'PY'
import json, re
from pathlib import Path

src = Path("/path/to/output_dir/SAMPLE_NAME_master_actionable.json")
raw = json.loads(src.read_text())

findings = []
for r in raw.get("records", raw if isinstance(raw, list) else []):
    hpo_raw = r.get("gene_hpo_id") or ""
    hpos = [h.strip() for h in re.split(r"[;,\s]+", hpo_raw) if h.strip().startswith("HP:")]
    findings.append({
        "gene_symbol": r.get("hugo") or r.get("gene_symbol") or "Unknown",
        "ncbi_description": r.get("ncbi_description") or r.get("gene_desc"),
        "rsid": r.get("rsid") or r.get("dbsnp"),
        "chromosome": r.get("chrom") or r.get("chromosome") or "",
        "position": int(r.get("pos") or r.get("position") or 0),
        "genotype": r.get("genotype") or r.get("alleles") or "N/A",
        "zygosity": r.get("zygosity") or ("Heterozygous" if "/" in str(r.get("genotype","")) else "Unknown"),
        "revel_score": r.get("revel"),
        "impact_consequence": r.get("so") or r.get("consequence") or r.get("impact_consequence") or "—",
        "clinvar_significance": r.get("clinvar_sig") or r.get("clinvar_significance"),
        "phasing": (r.get("phasing") or r.get("phase") or "undetermined").lower(),
        "associated_hpo_terms": hpos,
        "associated_mondo_terms": [],
        "gnomad_af": r.get("gnomad4_af") or r.get("gnomad_af"),
        "last_evaluated": r.get("clinvar_date") or "",
    })

out = {
    "patient_id": raw.get("patient") or "SAMPLE_NAME",
    "run_date": raw.get("run_date") or "",
    "monogenic_findings": findings,
    "polygenic_findings": [],
    "pharma_findings": [],
}
Path("/path/to/output_dir/SAMPLE_NAME_variant_report.json").write_text(json.dumps(out, indent=2))
print("Wrote VariantReport JSON")
PY

# 4. Render
python render_visual_ontology_explorer.py \
    -i /path/to/output_dir/SAMPLE_NAME_variant_report.json \
    -o reports/SAMPLE_NAME_visual_explorer.html
```

### OpenCRAVAT / SQLite requirements (already used by this repo)

| Resource | Typical path / note |
|----------|---------------------|
| OpenCRAVAT modules | `hpo`, `ncbigene`, `clinvar`, `revel`, `dbsnp`, `gwas_catalog`, … |
| Local NCBI Gene SQLite | `/data/opencravat/modules/annotators/ncbigene/data/ncbigene.sqlite` (used by `lib/enrich_report.py`) |
| Input | VCF / gVCF **or** existing OpenCRAVAT `.sqlite` |
| Phasing | Populate `phasing` field (WhatsHap, SHAPEIT, or long-range LD); values: `maternal`, `paternal`, `de_novo`, `undetermined` |

If you already have an annotated OpenCRAVAT SQLite or the `*_master_actionable.json` from a previous run, you only need steps 2–4 above.

---

## Commit the upgrade

```bash
cd ~/My-Projects/genomics-ontology/genomics-ontology
git checkout grok-dev
# copy upgraded render_visual_ontology_explorer.py + GROK_DEV_README.md into repo root
git add render_visual_ontology_explorer.py GROK_DEV_README.md samples/ reports/
git commit -m "feat: upgrade Visual Ontology Explorer to high-fidelity mockup UI"
git push origin grok-dev
```
