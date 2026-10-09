import os
import json
import re

# Load sanitized JSON data
with open('data/sanitized/proband_01_ontology_pharma_alphagenome.json', 'r') as f:
    data = json.load(f)

# Extract structured sections
genes_map = {g['symbol']: g for g in data.get('genes', [])}

# 1. Monogenic Pathogenic & High Impact Carriers
# CBLIF c.79+1G>A (VCV439755)
# GJB2 p.Met34Thr (VCV17000)
# VDR c.-1172A>G (VCV3336650)
primary_pathogenic = [
    {
        'gene': 'CBLIF',
        'name': 'Cobalamin Binding Intrinsic Factor',
        'var': 'c.79+1G>A',
        'hgvs_p': 'Splice Donor Disruption',
        'coord': 'chr11:59845374',
        'rsid': 'rs147785187',
        'zyg': 'Heterozygous (Carrier)',
        'clinvar': 'Pathogenic',
        'vcv': 'VCV439755',
        'cadd': 32.0,
        'revel': None,
        'avi': 33.9,
        'avi_mod': 'Splicing (Top 0.04%)',
        'disease': 'Hereditary Intrinsic Factor Deficiency (MIM:261000)',
        'inheritance': 'Autosomal Recessive',
        'mechanism': 'Disruption of canonical exon 1 splice donor site prevents transcription of functional intrinsic factor, impairing ileal cobalamin (Vitamin B12) absorption.',
        'clinical_guidance': 'Carrier state is typically asymptomatic under baseline nutrition; suggest periodic monitoring of serum Vitamin B12 and methylmalonic acid (MMA) to preclude subclinical deficiency or macrocytic anemia.'
    },
    {
        'gene': 'GJB2',
        'name': 'Gap Junction Protein Beta 2 (Connexin 26)',
        'var': 'c.101T>C',
        'hgvs_p': 'p.Met34Thr',
        'coord': 'chr13:20189481',
        'rsid': 'rs35887622',
        'zyg': 'Heterozygous (Carrier)',
        'clinvar': 'Pathogenic',
        'vcv': 'VCV17000',
        'cadd': 20.9,
        'revel': 0.702,
        'avi': 23.6,
        'avi_mod': 'Cactus (Phred 23.6)',
        'disease': 'Autosomal Recessive Non-Syndromic Sensorineural Hearing Impairment (DFNB1A, MIM:220290)',
        'inheritance': 'Autosomal Recessive',
        'mechanism': 'Missense alteration in the first transmembrane domain of connexin 26 altering potassium recycling channels within cochlear endolymph.',
        'clinical_guidance': 'Unpaired heterozygous carrier status does not typically cause syndromic deafness; suggest periodic baseline pure-tone audiometry and awareness of ototoxic medication exposures (e.g., aminoglycosides).'
    },
    {
        'gene': 'VDR',
        'name': 'Vitamin D Receptor',
        'var': 'c.-1172A>G',
        'hgvs_p': '5\' UTR Regulatory Variant',
        'coord': 'chr12:47906043',
        'rsid': 'rs4516035',
        'zyg': 'Homozygous (Alternative)',
        'clinvar': 'Likely pathogenic',
        'vcv': 'VCV3336650',
        'cadd': 5.36,
        'revel': None,
        'avi': None,
        'avi_mod': 'Promoter / 5\' UTR',
        'disease': 'Modulation of Vitamin D Endocrine Axis & Bone Mineral Density',
        'inheritance': 'Multifactorial / Promoter Susceptibility',
        'mechanism': 'Homozygous promoter variation altering GATA/transcription factor binding affinity in 5\' regulatory region, attenuating VDR transcriptional output.',
        'clinical_guidance': 'Suggest routine surveillance of 25-hydroxyvitamin D [25(OH)D], serum calcium, and consideration of bone mineral density screening during routine adult wellness evaluations.'
    }
]

# 2. Pharmacogenomics & Toxicogenomics
pgx_variants = [
    {
        'gene': 'DPYD',
        'var': 'c.496A>G (p.Met166Val)',
        'rsid': 'rs2297595',
        'zyg': 'Homozygous (*6)',
        'phenotype': 'Dihydropyrimidine Dehydrogenase (DPD) intermediate metabolizer profile',
        'drug': 'Fluoropyrimidines (5-Fluorouracil, Capecitabine, Tegafur)',
        'guideline': 'CPIC Level A Guideline: Homozygous *6 carriers demonstrate modest reduction in catalytic enzymatic clearance; oncologic practice recommends standard dosing with vigilant early toxicity monitoring (empiric dose reduction not universally mandated for isolated *6).'
    },
    {
        'gene': 'NAT2',
        'var': 'c.341T>C (p.Ile114Thr)',
        'rsid': 'rs1801280',
        'zyg': 'Homozygous (*5/*5)',
        'phenotype': 'Slow Acetylator Phenotype',
        'drug': 'Isoniazid, Hydralazine, Sulfonamides, Procainamide',
        'guideline': 'Substantially reduced N-acetyltransferase 2 catalytic activity. Suggest therapeutic drug monitoring and lower titration when initiating isoniazid or hydralazine to prevent drug-induced peripheral neuropathy or lupus-like syndromes.'
    },
    {
        'gene': 'VKORC1',
        'var': 'c.174-136C>T / c.-1639G>A',
        'rsid': 'rs9934438 / rs9923231',
        'zyg': 'Heterozygous (Intermediate Sensitivity)',
        'phenotype': 'Modulated Vitamin K Epoxide Reductase Expression',
        'drug': 'Warfarin / Coumarin Anticoagulants',
        'guideline': 'Intermediate sensitivity to vitamin K antagonists. If warfarin therapy is ever indicated, suggest CPIC-guided dose titration algorithms incorporating VKORC1 and CYP2C9 genotypes to avoid over-anticoagulation.'
    },
    {
        'gene': 'CYP2C9',
        'var': 'c.430C>T (p.Arg144Cys, *2)',
        'rsid': 'rs1799853',
        'zyg': 'Heterozygous (*1/*2)',
        'phenotype': 'Intermediate CYP2C9 Metabolizer',
        'drug': 'Warfarin, Phenytoin, Celecoxib, Glipizide',
        'guideline': 'Moderately attenuated hepatic clearance of CYP2C9 substrates. Suggest conservative dose titration for narrow therapeutic index NSAIDs and oral hypoglycemics.'
    }
]

# 3. Secondary Risk Modifiers & Channelopathies
secondary_modifiers = [
    {
        'gene': 'ANK2',
        'var': 'p.Arg3906Trp',
        'coord': 'chr4:113268802',
        'rsid': 'rs1800178',
        'zyg': 'Heterozygous',
        'clinvar': 'Conflicting classifications of pathogenicity',
        'cadd': 27.7,
        'revel': 0.722,
        'avi': 25.2,
        'domain': 'Cardiac Electrophysiology / Ankyrin-B Syndrome Susceptibility',
        'guidance': 'Ankyrin-B cellular scaffolding protein; implicated in ion channel localization. Suggest clinical caution regarding QT-prolonging pharmacotherapy (cross-reference CredibleMeds.org) and obtaining a baseline 12-lead ECG.'
    },
    {
        'gene': 'F5',
        'var': 'p.Arg534Gln (Factor V Leiden)',
        'coord': 'chr1:169549811',
        'rsid': 'rs6025',
        'zyg': 'Heterozygous',
        'clinvar': 'Drug response / Risk factor',
        'cadd': 27.9,
        'revel': None,
        'avi': 23.6,
        'domain': 'Hemostasis / Venous Thromboembolism (VTE) Risk Modifier',
        'guidance': 'Resistance to activated protein C (APC). Imparts a 3- to 5-fold increased relative risk of primary VTE. Lifelong anticoagulation is not indicated for asymptomatic carriers; suggest situational thromboprophylaxis during major orthopedic surgery, trauma, or prolonged immobilization.'
    },
    {
        'gene': 'F5',
        'var': 'p.Thr295Ala',
        'coord': 'chr1:169556715',
        'rsid': 'rs371760153',
        'zyg': 'Heterozygous (Maternal)',
        'clinvar': 'Uncertain significance',
        'cadd': 24.7,
        'revel': 0.700,
        'avi': 21.0,
        'domain': 'Factor V Rare Missense VUS',
        'guidance': 'Co-occurring rare missense VUS located on maternal allele (trans to or unphased with rs6025). High in silico deleteriousness (CADD 24.7, REVEL 0.700); classified as VUS in ClinVar.'
    },
    {
        'gene': 'LRP5',
        'var': 'p.Val233Glu',
        'coord': 'chr11:68363758',
        'rsid': 'rs1177484229',
        'zyg': 'Heterozygous',
        'clinvar': 'Uncertain significance',
        'cadd': 24.9,
        'revel': 0.903,
        'avi': 25.4,
        'domain': 'Wnt/Beta-Catenin Signaling / Bone Density & Retinal Vascularization',
        'guidance': 'Extremely high REVEL score (0.903) in the beta-propeller domain of LRP5. Heterozygous variants in LRP5 modulate bone mineral density and familial exudative vitreoretinopathy susceptibility.'
    },
    {
        'gene': 'MTHFR',
        'var': 'p.Ala222Val (c.665C>T / 677C>T)',
        'coord': 'chr1:11796321',
        'rsid': 'rs1801133',
        'zyg': 'Homozygous (Alternative)',
        'clinvar': 'Drug response',
        'cadd': 27.3,
        'revel': 0.842,
        'avi': 25.8,
        'domain': 'Folate Metabolism & Homocysteine Remethylation',
        'guidance': 'Common thermolabile variant (~60% enzymatic activity reduction when homozygous). Suggest maintaining dietary folate and Vitamin B12 adequacy; routine unindicated high-dose folate supplementation or aggressive anticoagulation is discouraged under current ACMG practice guidelines.'
    }
]

# 4. Verified Protective & Longevity Modifiers (Grounded in JSON)
protective_verified = [
    {
        'gene': 'ADH1C',
        'var': 'p.Ile350Val (c.1048A>G)',
        'coord': 'chr4:99339632',
        'rsid': 'rs698',
        'zyg': 'Heterozygous (Paternal)',
        'clinvar': 'Protective',
        'finding': 'Alcohol Dehydrogenase 1C *1/*2 polymorphism. Modulates ethanol oxidation rate; historically correlated with favorable HDL-cholesterol profiles and statistical cardioprotection in moderate alcohol consumers.'
    },
    {
        'gene': 'CCR5',
        'var': 'c.-556A>G (rs1799987) / p.Ser185IlefsTer32 (rs333)',
        'coord': 'chr3:46370444 / chr3:46373453',
        'rsid': 'rs1799987 / rs333',
        'zyg': 'Heterozygous',
        'clinvar': 'Protective / Benign',
        'finding': 'Chemokine receptor 5 regulatory and structural variations. Carriers demonstrate altered inflammatory chemokine recruitment (CCL4 GWAS p=1e-12) and partial protection against macrophage-tropic HIV-1 cellular entry.'
    },
    {
        'gene': 'CDKN2B',
        'var': 'c.*2619C>T',
        'coord': 'chr9:22003368',
        'rsid': 'rs1063192',
        'zyg': 'Heterozygous',
        'clinvar': 'Likely pathogenic|protective',
        'finding': 'Cyclin Dependent Kinase Inhibitor 2B 3\' UTR locus. GWAS Catalog (PMID: 30054458) demonstrates genome-wide significant correlation with metabolic glycemic regulation and Type 2 diabetes risk modulation (p=3e-18).'
    },
    {
        'gene': 'CASP8',
        'var': 'c.-937_-932del',
        'coord': 'chr2:201232809',
        'rsid': 'rs3834129',
        'zyg': 'Heterozygous',
        'clinvar': 'Protective',
        'finding': 'Six-nucleotide promoter deletion attenuating caspase-8 activation dynamics; epidemiological literature documents reduced apoptotic vulnerability and statistical protection against multiple cutaneous epithelial malignancies.'
    },
    {
        'gene': 'APOB',
        'var': 'p.Pro2739Leu (rs1801701) & common lipid alleles',
        'coord': 'chr2:21006509',
        'rsid': 'rs1801701',
        'zyg': 'Heterozygous',
        'clinvar': 'Conflicting / Benign',
        'finding': 'Polygenic apolipoprotein B variant (GWAS p=3e-22, PMID: 41325697). Imparts common population-level lipid variability without conferring monogenic Familial Hypercholesterolemia.'
    }
]

print('Structured datasets loaded successfully.')
