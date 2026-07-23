#!/usr/bin/env bash
# =============================================================================
# run_ontology_report.sh
# Generic ontology-driven (HPO + GO) actionable variant report. The phenotype
# DOMAIN is chosen with -c <config>; the engine itself is domain-agnostic.
#
# Usage:
#   ./run_ontology_report.sh [-c CONFIG] <INPUT_VCF|SQLITE> <OUTPUT_DIR> <PREFIX>
#
#   -c CONFIG   domain config yaml (default: config/cardiology.yaml)
#               e.g. -c config/hereditary_cancer.yaml
#
#   - If <INPUT> ends in .sqlite it is treated as an already-annotated OpenCRAVAT
#     database and Stage 1 (annotation) is skipped.
#   - Otherwise it is treated as a VCF and annotated first.
#
# Stages: 1 annotate -> 2 build panel -> 3 schema probe -> 4 actionable filter
#         -> 5 enrich (NCBI gene desc + live GWAS studies) -> 6 render
#         -> 7 native OpenCRAVAT Excel/VCF export
#
#   -c CONFIG   domain config (default config/cardiology.yaml)
#   -o          offline enrichment: use the cache only, never hit the network
#   -E          skip enrichment entirely
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

eval "$(micromamba shell hook --shell bash)"
micromamba activate cravat_env

CONFIG="$SCRIPT_DIR/config/cardiology.yaml"
OFFLINE=0        # -o : never hit the network during enrichment (cache only)
NO_ENRICH=0      # -E : skip the enrichment stage entirely
while getopts ":c:oE" opt; do
  case "$opt" in
    c) CONFIG="$OPTARG" ;;
    o) OFFLINE=1 ;;
    E) NO_ENRICH=1 ;;
    *) echo "Usage: $0 [-c CONFIG] [-o(ffline)] [-E no-enrich] <VCF|sqlite> <output_dir> <prefix>"; exit 1 ;;
  esac
done
shift $((OPTIND - 1))

ANNOTATORS="gnomad4 allofus250k spliceai cadd linsight clinvar clinvar_acmg clingen omim alphamissense revel cardioboost go gtex arrvars bayesdel esm1b varity_r metarnn ncer ccre_screen regulomedb pubmed hpo dbsnp gwas_catalog vcfinfo"

INPUT="${1:?Usage: run_ontology_report.sh [-c CONFIG] <VCF|sqlite> <output_dir> <prefix>}"
OUTDIR="${2:?output_dir required}"
PREFIX="${3:-Patient}"

[ -f "$CONFIG" ] || { echo "Config not found: $CONFIG"; exit 1; }

# Read domain + report settings (renderer, enrichment toggles) from the config.
read -r DOMAIN RENDERER ENRICH_GENES ENRICH_STUDIES < <(python3 - "$CONFIG" <<'PY'
import sys, yaml
cfg = yaml.safe_load(open(sys.argv[1])) or {}
rep = cfg.get("report", {}) or {}
enr = rep.get("enrichment", {}) or {}
print(cfg.get("domain", "domain"),
      rep.get("renderer", "generic"),
      "1" if enr.get("genes", True) else "0",
      "1" if enr.get("studies", False) else "0")
PY
)

mkdir -p "$OUTDIR"
RAW_DB="$OUTDIR/${PREFIX}.sqlite"
PANEL="$OUTDIR/${PREFIX}_${DOMAIN}_panel.json"
SCHEMA="$OUTDIR/${PREFIX}_schema.json"
ACT_DB="$OUTDIR/${PREFIX}_${DOMAIN}_actionable.sqlite"
ACT_JSON="$OUTDIR/${PREFIX}_${DOMAIN}_actionable.json"
ENRICH_CACHE="$OUTDIR/${PREFIX}_${DOMAIN}_enrich_cache.json"
HTML="$OUTDIR/${PREFIX}_${DOMAIN}_report.html"
TSV="$OUTDIR/${PREFIX}_${DOMAIN}_report.tsv"
TEXT="$OUTDIR/${PREFIX}_${DOMAIN}_report.txt"

echo "=================================================================="
echo "ONTOLOGY ACTIONABLE REPORT  |  domain=$DOMAIN  prefix=$PREFIX"
echo "config=$CONFIG"
echo "input=$INPUT"
echo "outdir=$OUTDIR"
echo "=================================================================="

# ---- Stage 1: annotation (skip if a .sqlite was passed) --------------------
case "$INPUT" in
  *.sqlite)
    echo "[1/7] Using existing annotated DB: $INPUT"
    RAW_DB="$INPUT"
    ;;
  *)
    echo "[1/7] Annotating VCF with OpenCRAVAT (this can take hours on WGS)..."
    oc run "$INPUT" -l hg38 -a $ANNOTATORS -d "$OUTDIR" --mp "$(nproc)" -n "$PREFIX"
    ;;
esac

echo "[2/7] Building ontology gene panel from HPO + GO..."
python3 lib/build_ontology_panel.py --config "$CONFIG" --out "$PANEL"

echo "[3/7] Probing database schema..."
python3 lib/schema_probe.py "$RAW_DB" --out "$SCHEMA"

echo "[4/7] Selecting actionable variants + tiering..."
python3 lib/ontology_filter.py \
  --raw-db "$RAW_DB" --panel "$PANEL" --schema "$SCHEMA" --config "$CONFIG" \
  --out-sqlite "$ACT_DB" --out-json "$ACT_JSON" --patient "$PREFIX"

# ---- Stage 5: enrichment (NCBI gene descriptions + live GWAS studies) ------
ENRICH_ARGS=""
[ "$ENRICH_GENES" = "1" ] && ENRICH_ARGS="$ENRICH_ARGS --genes"
[ "$ENRICH_STUDIES" = "1" ] && ENRICH_ARGS="$ENRICH_ARGS --studies"
[ "$OFFLINE" = "1" ] && ENRICH_ARGS="$ENRICH_ARGS --offline"
if [ "$NO_ENRICH" = "1" ] || [ -z "$ENRICH_ARGS" ]; then
  echo "[5/7] Enrichment skipped."
else
  STUDY_MSG=""; [ "$ENRICH_STUDIES" = "1" ] && STUDY_MSG=", live GWAS studies"
  echo "[5/7] Enriching (NCBI gene descriptions${STUDY_MSG})..."
  python3 lib/enrich_report.py --in-json "$ACT_JSON" --cache "$ENRICH_CACHE" $ENRICH_ARGS \
    || echo "  [warn] enrichment step reported an error (non-fatal, report still renders)"
fi

echo "[6/7] Rendering HTML / TSV / text deliverables (renderer: $RENDERER)..."
if [ "$RENDERER" = "autoimmune" ]; then
  python3 lib/render_autoimmune.py \
    --in-json "$ACT_JSON" --out-html "$HTML" --out-tsv "$TSV" --out-text "$TEXT"
else
  python3 lib/render_report.py \
    --in-json "$ACT_JSON" --out-html "$HTML" --out-tsv "$TSV" --out-text "$TEXT"
fi

echo "[7/7] Native OpenCRAVAT Excel + VCF export of actionable variants..."
FILTERSQL="$(python3 lib/make_filtersql.py --in-json "$ACT_JSON" --schema "$SCHEMA")"
if [ -n "$FILTERSQL" ]; then
  oc report "$RAW_DB" -t excel vcf --filtersql "$FILTERSQL" \
    -s "$OUTDIR/${PREFIX}_${DOMAIN}_actionable" --silent || \
    echo "  [warn] native oc report export skipped (non-fatal)"
else
  echo "  [info] no actionable variants; skipping native export"
fi

echo "------------------------------------------------------------------"
echo "DONE ($DOMAIN). Key deliverables:"
echo "  HTML  : $HTML"
echo "  TSV   : $TSV"
echo "  Text  : $TEXT"
echo "  Panel : $PANEL"
echo "  Excel : $OUTDIR/${PREFIX}_${DOMAIN}_actionable.xlsx  (native OpenCRAVAT)"
echo "  VCF   : $OUTDIR/${PREFIX}_${DOMAIN}_actionable.vcf   (native OpenCRAVAT)"
echo "------------------------------------------------------------------"
