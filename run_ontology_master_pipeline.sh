#!/usr/bin/env bash
# =============================================================================
# run_ontology_master_pipeline.sh
# End-to-End Execution Wrapper for Genomic Ontology Reports
# Supports: .vcf, .vcf.gz, .sqlite, or OpenCRAVAT Job IDs
# Automatically outputs to Google Drive: ~/Google Drive/My Drive/Ontology/{Sample_ID}-{DD-MM-YYYY}/
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

usage() {
    echo "Usage: $0 [options] <sample_id> <input_source>"
    echo ""
    echo "Arguments:"
    echo "  sample_id      Patient/Sample identifier (e.g. DE_master, HG003)"
    echo "  input_source   Input: .vcf, .vcf.gz, .sqlite, or OC Job ID (e.g. 260706-105810)"
    echo ""
    echo "Options:"
    echo "  -c CONFIG      Domain YAML config (default: config/ontology_domains.yaml)"
    echo "  -g GDRIVE_DIR  Google Drive base directory (default: ~/Google Drive/My Drive/Ontology)"
    echo "  -l             Local only (skip Google Drive sync)"
    echo "  -h             Show this help message"
    echo ""
    echo "Examples:"
    echo "  $0 DE_master 260706-105810"
    echo "  $0 DE_master /data/opencravat/jobs/default/260706-105810/DE_master_phased_final.UCSC.vcf.gz.sqlite"
    echo "  $0 Patient101 sample.vcf.gz"
    echo ""
    exit 1
}

CONFIG="$SCRIPT_DIR/config/ontology_domains.yaml"
GDRIVE_DIR="$HOME/Google Drive/My Drive/Ontology"
LOCAL_ONLY=""

while getopts ":c:g:lh" opt; do
    case "$opt" in
        c) CONFIG="$OPTARG" ;;
        g) GDRIVE_DIR="$OPTARG" ;;
        l) LOCAL_ONLY="--local-only" ;;
        h) usage ;;
        *) usage ;;
    esac
done
shift $((OPTIND - 1))

if [ "$#" -lt 2 ]; then
    usage
fi

SAMPLE_ID="$1"
INPUT_SOURCE="$2"

python3 "$SCRIPT_DIR/run_ontology_pipeline.py" \
    --sample "$SAMPLE_ID" \
    --input "$INPUT_SOURCE" \
    --config "$CONFIG" \
    --gdrive-dir "$GDRIVE_DIR" \
    $LOCAL_ONLY
