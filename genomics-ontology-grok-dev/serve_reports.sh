#!/bin/bash
# serve_reports.sh — Local Web Server for Genomics Reports (PC, iPhone, iPad)

PORT=8080
IP=$(hostname -I | awk '{print $1}')
DIR="/data/Genomes/DE/Data/Final/2026-03-29"

echo "=================================================================="
echo "   GENOMICS REPORT MULTI-DEVICE WEB SERVER"
echo "=================================================================="
echo " Local Network IP : $IP"
echo " Server Port      : $PORT"
echo " Root Directory   : $DIR"
echo "=================================================================="
echo " To view on PC, iPhone, or iPad connected to Wi-Fi, open Safari/Chrome:"
echo ""
echo "   Master Ontology Hub : http://$IP:$PORT/ontology_reports/DE_master_hub_master_ontology_report.html"
echo "   Pharmacogenomics    : http://$IP:$PORT/pharma_reports/DE_patient_pharma_report.html"
echo "=================================================================="

python3 -m http.server $PORT --directory "$DIR"
