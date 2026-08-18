#!/usr/bin/env bash
# =============================================================================
# run_ontology_master_pipeline.sh
# Complete pipeline: VCF -> OpenCRAVAT Annotation -> Multi-Domain HPO Level 1+2 Reports
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Ensure python environment has OpenCRAVAT
if command -v micromamba >/dev/null 2>&1; then
    eval "$(micromamba shell hook --shell bash)"
    micromamba activate cravat_env 2>/dev/null || true
fi

# Usage help
usage() {
    echo "Usage: $0 [options] <VCF|sqlite> <output_dir> <sample_name>"
    echo ""
    echo "Options:"
    echo "  -c CONFIG       Custom domain YAML config (default: config/ontology_domains.yaml)"
    echo "  -R RENDERER     Renderer style: master (default), autoimmune, browser, glass, dashboard"
    echo "  -o              Offline mode (use local ncbigene & cached GWAS evidence)"
    echo "  -h              Show this help message"
    echo ""
    exit 1
}

CONFIG="$SCRIPT_DIR/config/ontology_domains.yaml"
RENDERER="master"
OFFLINE=0

while getopts ":c:R:oh" opt; do
    case "$opt" in
        c) CONFIG="$OPTARG" ;;
        R) RENDERER="$OPTARG" ;;
        o) OFFLINE=1 ;;
        h) usage ;;
        *) usage ;;
    esac
done
shift $((OPTIND - 1))

if [ "$#" -lt 3 ]; then
    usage
fi

INPUT="$1"
OUTDIR="$2"
PREFIX="$3"

ANNOTATORS="hpo go clinvar clingen omim ncbigene revel alphamissense bayesdel metarnn esm1b varity spliceai cadd linsight ncer regulomedb ccre_screen gtex dbsnp vcfinfo gwas_catalog pharmgkb civic interpro"

mkdir -p "$OUTDIR"
RAW_DB="$OUTDIR/${PREFIX}.sqlite"
SCHEMA="$OUTDIR/${PREFIX}_schema.json"
PANEL="$OUTDIR/${PREFIX}_master_panel.json"
ACT_DB="$OUTDIR/${PREFIX}_master_actionable.sqlite"
ACT_JSON="$OUTDIR/${PREFIX}_master_actionable.json"
ENRICH_CACHE="$OUTDIR/${PREFIX}_enrich_cache.json"

HTML="$OUTDIR/${PREFIX}_master_ontology_report.html"
TSV="$OUTDIR/${PREFIX}_master_ontology_report.tsv"
TEXT="$OUTDIR/${PREFIX}_master_ontology_report.txt"

echo "=================================================================="
echo "ONTOLOGY MASTER PIPELINE  |  sample=$PREFIX  renderer=$RENDERER"
echo "input=$INPUT"
echo "outdir=$OUTDIR"
echo "=================================================================="

# -----------------------------------------------------------------------------
# Stage 1: Annotation (OpenCRAVAT)
# -----------------------------------------------------------------------------
case "$INPUT" in
    *.sqlite)
        echo "[1/6] Using pre-annotated SQLite database: $INPUT"
        RAW_DB="$INPUT"
        ;;
    *.vcf|*.vcf.gz|*.g.vcf|*.g.vcf.gz|*)
        echo "[1/6] Running OpenCRAVAT annotation on VCF/g.VCF input ($INPUT)..."
        oc run "$INPUT" -l hg38 -a $ANNOTATORS -d "$OUTDIR" --mp "$(nproc 2>/dev/null || echo 4)" -n "$PREFIX"
        RAW_DB="$OUTDIR/${PREFIX}.sqlite"
        ;;
esac

# -----------------------------------------------------------------------------
# Stage 2: Panel Construction & Schema Probe
# -----------------------------------------------------------------------------
echo "[2/6] Building HPO + GO multi-domain gene panel..."
python3 lib/build_ontology_panel.py --config "$CONFIG" --out "$PANEL"

echo "[3/6] Probing database schema & phasing columns..."
python3 lib/schema_probe.py "$RAW_DB" --out "$SCHEMA"

# -----------------------------------------------------------------------------
# Stage 3: Variant Selection & Tiering
# -----------------------------------------------------------------------------
echo "[4/6] Filtering & tiering actionable variants..."
python3 lib/ontology_filter.py \
    --raw-db "$RAW_DB" --panel "$PANEL" --schema "$SCHEMA" --config "$CONFIG" \
    --out-sqlite "$ACT_DB" --out-json "$ACT_JSON" --patient "$PREFIX"

# -----------------------------------------------------------------------------
# Stage 4: Local & Remote Enrichment
# -----------------------------------------------------------------------------
ENRICH_ARGS="--genes"
[ "$OFFLINE" -eq 1 ] && ENRICH_ARGS="$ENRICH_ARGS --offline"

echo "[5/6] Enriching gene descriptions & OMIM clinical synopses..."
python3 lib/enrich_report.py --in-json "$ACT_JSON" --cache "$ENRICH_CACHE" $ENRICH_ARGS || true

# -----------------------------------------------------------------------------
# Stage 5: Rendering Multi-Domain Master Deliverables
# -----------------------------------------------------------------------------
echo "[6/6] Rendering deliverables (Universal Master Hub + HPO Level 1 & 2 selector)..."
case "$RENDERER" in
    master|hub|portal)
        python3 lib/render_master_hub.py \
            --in-json "$ACT_JSON" --out-html "$HTML" --out-tsv "$TSV" --out-text "$TEXT" \
            --domain-config "$CONFIG"
        ;;
    autoimmune)
        python3 lib/render_autoimmune.py \
            --in-json "$ACT_JSON" --out-html "$HTML" --out-tsv "$TSV" --out-text "$TEXT"
        ;;
    browser)
        python3 lib/render_gene_browser.py \
            --in-json "$ACT_JSON" --out-html "$HTML" --out-tsv "$TSV" --out-text "$TEXT"
        ;;
    *)
        python3 lib/render_report.py \
            --in-json "$ACT_JSON" --out-html "$HTML" --out-tsv "$TSV" --out-text "$TEXT"
        ;;
esac

# -----------------------------------------------------------------------------
# Stage 6: Google Drive Delivery ("Ontology" folder)
# -----------------------------------------------------------------------------
echo "[7/7] Delivering reports to Google Drive ('Ontology')..."
python3 "$SCRIPT_DIR/cloud_delivery_service.py" "$HTML" "Ontology" || true
python3 "$SCRIPT_DIR/cloud_delivery_service.py" "$TSV" "Ontology" || true
python3 "$SCRIPT_DIR/cloud_delivery_service.py" "$TEXT" "Ontology" || true

echo "=================================================================="
echo "SUCCESS. Master Ontology Deliverables:"
echo "  Universal Master Portal HTML : $HTML"
echo "  Actionable Variant TSV       : $TSV"
echo "  Clinical Text Summary        : $TEXT"
echo "  Google Drive Location        : ~/Google Drive/My Drive/Ontology/"
echo "=================================================================="
