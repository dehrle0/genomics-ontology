import json
import datetime
import os

def export_ehr(json_path, out_path, subject_name='Daniel Ehrle', subject_id='Daniel_Ehrle'):
    with open(json_path, 'r', encoding='utf-8') as f:
        data = json.load(f)

    ehr_bundle = {
        'resourceType': 'Bundle',
        'id': f'genomics-ehr-{subject_id}',
        'type': 'collection',
        'timestamp': datetime.datetime.now().isoformat(),
        'meta': {
            'profile': ['http://hl7.org/fhir/uv/genomics-reporting/StructureDefinition/genomics-report'],
            'source': 'Genomic Ontology Reporting System v5.3 (AlphaGenome Enhanced)'
        },
        'patient': {
            'resourceType': 'Patient',
            'id': subject_id,
            'name': [{'use': 'official', 'text': subject_name}],
            'meta': {
                'source': 'Self-Reported & Client Phased Whole-Genome Sequencing'
            }
        },
        'reportMetadata': {
            'reportType': 'Clinical Genomics Decision Support & Self-Reported Evidence Summary',
            'intendedUse': 'Non-Diagnostic Exploratory Research & Clinical Consultation Only',
            'disclaimer': 'This record contains patient-initiated genomic analysis compiled using local open-weight foundation models. It is not an in vitro diagnostic test. Orthogonal clinical validation is recommended prior to therapeutic changes.',
            'sequencingModality': '40x Whole-Genome Sequencing (WGS, GRCh38.p14)',
            'generatedDate': datetime.datetime.now().strftime('%Y-%m-%d')
        },
        'actionableGenomicAlerts': [],
        'pharmacogenomicPrescribingAlerts': [],
        'organSystemRiskSummary': [],
        'polygenicRiskScoresSelfReported': []
    }

    # Extract Primary Findings
    for g in data.get('genes', []):
        sym = g.get('symbol')
        for v in g.get('variants', []):
            cv = str(v.get('clinvar'))
            cat = str(v.get('category'))
            if 'pathogenic' in cv.lower() and 'conflicting' not in cv.lower():
                ehr_bundle['actionableGenomicAlerts'].append({
                    'gene': sym,
                    'variantHgvsC': v.get('cchange'),
                    'variantHgvsP': v.get('achange'),
                    'coordinate': v.get('coordinate'),
                    'dbSnp': v.get('id'),
                    'zygosity': v.get('zygosity'),
                    'clinicalSignificance': cv,
                    'clinvarAccession': f'VCV{v.get("clinvarId")}' if v.get('clinvarId') else None,
                    'interpretation': 'Monogenic carrier state; asymptomatic under baseline conditions. Low-cost screening suggested.',
                    'alphagenomeAvi': v.get('aviPhred')
                })

    # Extract PGx Alerts
    pgx_interactions = data.get('pharmacogenomics', {}).get('interactions', [])
    for p in pgx_interactions[:15]:
        ehr_bundle['pharmacogenomicPrescribingAlerts'].append({
            'gene': p.get('gene'),
            'variant': p.get('variant'),
            'medication': p.get('drug') or p.get('medication'),
            'phenotype': p.get('phenotype'),
            'clinicalImplication': p.get('implication') or p.get('guideline') or p.get('action'),
            'cpicLevel': p.get('cpicLevel', 'Level A/B')
        })

    # Extract Organ System Risk
    for o in data.get('organRiskMatrix', []):
        ehr_bundle['organSystemRiskSummary'].append({
            'organSystem': o.get('system'),
            'riskTier': o.get('riskTier'),
            'polygenicPercentile': o.get('prsPercentile'),
            'pathogenicAlleleCount': o.get('pathogenicCount'),
            'primaryPathway': o.get('pathway'),
            'genesOfConcern': o.get('concernGenes')
        })

    # Extract top polygenic traits
    for prs in data.get('polygenicRisk', [])[:25]:
        ehr_bundle['polygenicRiskScoresSelfReported'].append({
            'trait': prs.get('trait'),
            'organSystem': prs.get('organSystem'),
            'percentile': prs.get('percentile'),
            'riskCategory': prs.get('category'),
            'pgsCatalogOrSource': prs.get('pgsId')
        })

    with open(out_path, 'w', encoding='utf-8') as f:
        json.dump(ehr_bundle, f, indent=2)
    print(f"[EHR Export] Written JSON bundle to {out_path} ({len(ehr_bundle['actionableGenomicAlerts'])} genomic alerts, {len(ehr_bundle['pharmacogenomicPrescribingAlerts'])} PGx alerts).")

if __name__ == '__main__':
    export_ehr('data/sanitized/proband_01_ontology_pharma_alphagenome.json', 'reports/alternatives/proband_01_ehr_import.json', 'PROBAND_01', 'PROBAND_WGS_40X')
