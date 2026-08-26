/**
 * MOCK DATA — Genomic Ontology Explorer
 * ---------------------------------------------------------------------
 * Stands in for the real backend (OpenCRAVAT sqlite output +
 * enrich_report.py caches). Field names track the OpenCRAVAT annotators
 * actually run on the job (see JOB_META.annotators below) so this is a
 * direct swap point, not a rewrite, when a real API is available:
 *
 *   GET /api/genes                -> GENES (summary fields only)
 *   GET /api/genes/:symbol        -> one GENES entry, full detail
 *   GET /api/ontology/:type       -> ONTOLOGIES[type]
 *   GET /api/analysis             -> PRS + PGX
 *   GET /api/report               -> REPORT
 *
 * Job metadata: HG003_GRCh38_1_22_v4.2.1_benchmark (GIAB benchmark,
 * public), OpenCRAVAT 3.1.1, run 2026-07-15. NOTE: this session has no
 * access to that job's actual .sqlite output (private server, no
 * upload, no network) — every score/finding below is illustrative,
 * shaped to match real annotator output fields.
 *
 * Reference-link identifiers (NCBI Gene ID, OMIM gene/phenotype
 * numbers) are real public identifiers for these genes, corrected in
 * this revision — see CHANGELOG at the bottom of this file.
 * ---------------------------------------------------------------------
 */

const JOB_META = {
  sample: "HG003_GRCh38_1_22_v4.2.1_benchmark",
  opencravatVersion: "3.1.1",
  submitted: "2026-07-15 23:09:35",
  uniqueVariants: 4045427,
  annotators: [
    "allofus250k", "alphamissense", "arrvars", "bayesdel", "cadd",
    "cardioboost", "ccre_screen", "clinvar", "clinvar_acmg", "dbsnp",
    "esm1b", "gerp", "gnomad4", "go", "gtex", "gwas_catalog", "hpo",
    "linsight", "litvar", "metarnn", "ncbigene", "ncer", "omim",
    "phastcons", "phylop", "pubmed", "regulomedb", "revel", "spliceai",
    "varity_r"
  ]
};

function refLinks(symbol, ncbiGeneId, omimGene) {
  return {
    ncbiGene: "https://www.ncbi.nlm.nih.gov/gene/" + ncbiGeneId,
    omim: "https://omim.org/entry/" + omimGene,
    genecards: "https://www.genecards.org/cgi-bin/carddisp.pl?gene=" + symbol,
    clinvarGene: "https://www.ncbi.nlm.nih.gov/clinvar/?term=" + symbol + "%5Bgene%5D"
  };
}

// ---------------------------------------------------------------------
// GENES
// ---------------------------------------------------------------------
const GENES = [
  {
    symbol: "VDR",
    name: "Vitamin D Receptor",
    chromosome: "12q13.11",
    organSystem: "Endocrine",
    ncbiGeneId: "7421",
    omimGene: "601769",
    omimPhenotype: "277440",
    links: refLinks("VDR", "7421", "601769"),
    summary: "The Vitamin D receptor (VDR) gene plays a crucial role in the body's use of vitamin D, which is vital for bone health and other cellular functions. The protein coded by the VDR gene acts as an intracellular hormone receptor, binding to the active form of vitamin D, 1,25-dihydroxyvitamin D3 (calcitriol), and facilitating its biological effects. When the VDR gene is dysfunctional, it can lead to a condition known as Vitamin D-dependent rickets type 2A (VDDR2A).",
    associatedPathology: [
      { name: "Rickets, vitamin D-resistant, type IIA", inheritance: "Autosomal recessive", omim: "277440" }
    ],
    pli: 0.02,
    loeuf: 1.14,
    variantsDetected: 161,
    researchedVariants: 34,
    hpoTermCount: 7,
    goTermCount: 2,
    variants: [
      {
        id: "rs987849", genotype: "G/A", zygosity: "Heterozygous", phase: "Unknown",
        maf: 0.60, coordinate: "chr12:47844974", consequence: ["Intron", "NMD transcript variant"],
        category: "concern", clinvar: "Conflicting", revel: 0.21, cadd: 12.3, spliceai: 0.02, alphamissense: null,
        qual: 712.4, reads: { matching: 41, total: 63 }, lastEvaluated: "2025-02-11",
        studies: [
          { finding: "Associated with a 1.53-fold increased risk of multiple sclerosis", condition: "Multiple sclerosis", genotypeRelevance: "Your G/A fully matches the studied rs987849 risk allele", evidenceLevel: 2, source: "Huang J, et al. J Neuroimmunol 2024" },
          { finding: "Increased risk of multiple sclerosis in an independent replication cohort", condition: "Multiple sclerosis", genotypeRelevance: "Your G/A fully matches the studied rs987849 risk allele", evidenceLevel: 2, source: "GWAS Catalog replication set, 2023" }
        ]
      },
      {
        id: "rs886441", genotype: "G/A", zygosity: "Heterozygous", phase: "Paternal",
        maf: 0.77, coordinate: "chr12:47845221", consequence: ["Intron", "NMD transcript variant"],
        category: "protective", clinvar: "Benign", revel: 0.05, cadd: 3.1, spliceai: 0.00, alphamissense: null,
        qual: 891.0, reads: { matching: 38, total: 55 }, lastEvaluated: "2024-11-02",
        studies: [
          { finding: "G allele associated with reduced non-Hodgkin lymphoma risk", condition: "Non-Hodgkin lymphoma", genotypeRelevance: "Your G/A carries the protective G allele", evidenceLevel: 3, source: "Lindqvist A, et al. Cancer Epidemiol 2023" }
        ]
      },
      {
        id: "rs2544043", genotype: "C/C", zygosity: "Homozygous", phase: "N/A",
        maf: 0.95, coordinate: "chr12:47846010", consequence: ["3' UTR"],
        category: "uncertain", clinvar: "Uncertain significance", revel: 0.03, cadd: 5.6, spliceai: 0.01, alphamissense: null,
        qual: 655.8, reads: { matching: 60, total: 61 }, lastEvaluated: "2023-08-19",
        studies: [
          { finding: "No significant differences in genotype or allele frequencies for VDR rs2544043 (G>C) were observed between individuals with diabetic nephropathy and controls", condition: "Diabetic nephropathy", genotypeRelevance: "Your C/C matches the studied genotype; no association found", evidenceLevel: 2, source: "Renal Genomics Consortium, 2022" }
        ]
      },
      { id: "rs3847987", genotype: "C/A", zygosity: "Heterozygous", phase: "Maternal", maf: 0.13, coordinate: "chr12:47847112", consequence: ["3' UTR"], category: "uncategorized", clinvar: "Not reviewed", revel: 0.04, cadd: 4.0, spliceai: 0.00, alphamissense: null, qual: 601.2, reads: { matching: 22, total: 47 }, lastEvaluated: null, studies: [] },
      { id: "rs11168293", genotype: "G/T", zygosity: "Heterozygous", phase: "Unknown", maf: 0.28, coordinate: "chr12:47848330", consequence: ["5' UTR"], category: "uncategorized", clinvar: "Not reviewed", revel: 0.06, cadd: 6.2, spliceai: 0.01, alphamissense: null, qual: 733.9, reads: { matching: 29, total: 52 }, lastEvaluated: null, studies: [] },
      { id: "rs9729", genotype: "G/T", zygosity: "Heterozygous", phase: "Paternal", maf: 0.52, coordinate: "chr12:47849501", consequence: ["3' UTR"], category: "uncategorized", clinvar: "Not reviewed", revel: 0.02, cadd: 2.8, spliceai: 0.00, alphamissense: null, qual: 812.5, reads: { matching: 33, total: 58 }, lastEvaluated: null, studies: [] },
      { id: "rs2228570", genotype: "G/G", zygosity: "Homozygous", phase: "N/A", maf: 0.63, coordinate: "chr12:47844622", consequence: ["Missense", "5' UTR"], category: "uncategorized", clinvar: "Not reviewed", revel: 0.31, cadd: 14.7, spliceai: 0.03, alphamissense: 0.28, qual: 940.1, reads: { matching: 57, total: 59 }, lastEvaluated: null, studies: [] },
      { id: "rs739837", genotype: "G/T", zygosity: "Heterozygous", phase: "Maternal", maf: 0.53, coordinate: "chr12:47850980", consequence: ["3' UTR"], category: "uncategorized", clinvar: "Not reviewed", revel: 0.04, cadd: 3.4, spliceai: 0.00, alphamissense: null, qual: 688.0, reads: { matching: 25, total: 49 }, lastEvaluated: null, studies: [] }
    ],
    hpoTerms: [
      { id: "HP:0002748", label: "Rickets", evidence: "Gene-level (OMIM 277440)" },
      { id: "HP:0004349", label: "Reduced bone mineral density", evidence: "Gene-level" },
      { id: "HP:0000007", label: "Autosomal recessive inheritance", evidence: "Curated" },
      { id: "HP:0002960", label: "Autoimmunity", evidence: "Variant-level (rs987849, MS association)" }
    ],
    publications: [
      { title: "VDR gene polymorphisms and multiple sclerosis susceptibility: a meta-analysis", authors: "Huang J, et al.", journal: "J Neuroimmunol", year: 2024, doi: "10.1016/j.jneuroim.2024.578112", tags: ["Large Cohort", "ClinVar"], finding: "Pooled analysis across 9 cohorts found a modest but consistent MS risk increase for the rs987849 A allele." },
      { title: "Vitamin D receptor variation and lymphoma risk in a Scandinavian registry", authors: "Lindqvist A, et al.", journal: "Cancer Epidemiol", year: 2023, doi: "10.1016/j.canep.2023.102289", tags: ["Phase-aware"], finding: "Registry-linked analysis suggested a protective association for the VDR G allele at rs886441." }
    ]
  },
  {
    symbol: "MTHFR",
    name: "Methylenetetrahydrofolate Reductase",
    chromosome: "1p36.22",
    organSystem: "Metabolic",
    ncbiGeneId: "4524",
    omimGene: "607093",
    omimPhenotype: "236250",
    links: refLinks("MTHFR", "4524", "607093"),
    summary: "MTHFR encodes the enzyme that converts 5,10-methylenetetrahydrofolate to 5-methyltetrahydrofolate, the primary circulating form of folate used in homocysteine remethylation. Reduced-function variants are common and associated with mild-to-moderate hyperhomocysteinemia.",
    associatedPathology: [
      { name: "Homocystinuria due to MTHFR deficiency", inheritance: "Autosomal recessive", omim: "236250" }
    ],
    pli: 0.01,
    loeuf: 0.98,
    variantsDetected: 89,
    researchedVariants: 21,
    hpoTermCount: 5,
    goTermCount: 1,
    variants: [
      {
        id: "rs1801133", genotype: "A/G", zygosity: "Heterozygous", phase: "Maternal",
        maf: 0.34, coordinate: "chr1:11796321", consequence: ["Missense"],
        category: "concern", clinvar: "Risk factor", revel: 0.44, cadd: 19.8, spliceai: 0.01, alphamissense: 0.51,
        qual: 878.3, reads: { matching: 34, total: 60 }, lastEvaluated: "2025-05-30",
        studies: [
          { finding: "C677T variant associated with mildly elevated homocysteine and modestly increased cardiovascular risk", condition: "Cardiovascular disease", genotypeRelevance: "Your A/G is heterozygous for the reduced-function 677T allele", evidenceLevel: 3, source: "Chen R, et al. Atherosclerosis 2024" }
        ]
      },
      {
        id: "rs1801131", genotype: "T/T", zygosity: "Homozygous", phase: "N/A",
        maf: 0.31, coordinate: "chr1:11794419", consequence: ["Missense"],
        category: "uncertain", clinvar: "Uncertain significance", revel: 0.18, cadd: 11.2, spliceai: 0.00, alphamissense: 0.22,
        qual: 903.7, reads: { matching: 58, total: 61 }, lastEvaluated: "2024-01-14",
        studies: [
          { finding: "A1298C variant shows inconsistent association with neural tube defect risk across studies", condition: "Neural tube defects", genotypeRelevance: "Your T/T is homozygous wild-type at this position", evidenceLevel: 2, source: "Folate Genetics Working Group, 2021" }
        ]
      }
    ],
    hpoTerms: [
      { id: "HP:0003166", label: "Hyperhomocysteinemia", evidence: "Gene-level" },
      { id: "HP:0001677", label: "Coronary artery atherosclerosis", evidence: "Curated" }
    ],
    publications: [
      { title: "MTHFR C677T polymorphism and cardiovascular risk: updated meta-analysis", authors: "Chen R, et al.", journal: "Atherosclerosis", year: 2024, doi: "10.1016/j.atherosclerosis.2024.117201", tags: ["Large Cohort"], finding: "Confirmed a small but significant homocysteine elevation with the 677T allele, effect size unchanged from prior estimates." }
    ]
  },
  {
    symbol: "AGXT2",
    name: "Alanine--Glyoxylate Aminotransferase 2",
    chromosome: "5p13.2",
    organSystem: "Metabolic",
    ncbiGeneId: "64902",
    omimGene: "614696",
    omimPhenotype: null,
    links: refLinks("AGXT2", "64902", "614696"),
    summary: "AGXT2 encodes a mitochondrial enzyme involved in the metabolism of asymmetric dimethylarginine (ADMA) and beta-aminoisobutyrate. Variants affect circulating ADMA levels, a marker linked to endothelial function, and have been studied as modifiers of cardiovascular and renal traits rather than as a classic single-gene disease.",
    associatedPathology: [],
    pli: 0.00,
    loeuf: 1.32,
    variantsDetected: 47,
    researchedVariants: 9,
    hpoTermCount: 2,
    goTermCount: 1,
    variants: [
      { id: "rs37370", genotype: "A/G", zygosity: "Heterozygous", phase: "Paternal", maf: 0.41, coordinate: "chr5:35035678", consequence: ["Missense"], category: "uncertain", clinvar: "Uncertain significance", revel: 0.09, cadd: 9.4, spliceai: 0.00, alphamissense: 0.11, qual: 744.6, reads: { matching: 27, total: 51 }, lastEvaluated: "2022-10-03", studies: [
        { finding: "Missense variant associated with modestly elevated circulating ADMA levels", condition: "Endothelial dysfunction marker", genotypeRelevance: "Your A/G is heterozygous for the ADMA-elevating allele", evidenceLevel: 2, source: "Cardiometabolic Traits Consortium, 2022" }
      ] },
      { id: "rs37369", genotype: "G/G", zygosity: "Homozygous", phase: "N/A", maf: 0.55, coordinate: "chr5:35036201", consequence: ["Intron"], category: "uncategorized", clinvar: "Not reviewed", revel: 0.01, cadd: 1.9, spliceai: 0.00, alphamissense: null, qual: 812.0, reads: { matching: 49, total: 50 }, lastEvaluated: null, studies: [] }
    ],
    hpoTerms: [
      { id: "HP:0003259", label: "Elevated ADMA / renal function marker", evidence: "Variant-level" }
    ],
    publications: []
  },
  {
    symbol: "CBS",
    name: "Cystathionine Beta-Synthase",
    chromosome: "21q22.3",
    organSystem: "Metabolic",
    ncbiGeneId: "875",
    omimGene: "613381",
    omimPhenotype: "236200",
    links: refLinks("CBS", "875", "613381"),
    summary: "CBS catalyzes the first committed step of the transsulfuration pathway, converting homocysteine to cystathionine. Loss-of-function variants cause classical homocystinuria; common variants have smaller effects on homocysteine levels.",
    associatedPathology: [
      { name: "Homocystinuria, CBS-related", inheritance: "Autosomal recessive", omim: "236200" }
    ],
    pli: 0.34,
    loeuf: 0.71,
    variantsDetected: 122,
    researchedVariants: 18,
    hpoTermCount: 7,
    goTermCount: 1,
    variants: [
      {
        id: "rs5742905", genotype: "C/T", zygosity: "Heterozygous", phase: "Unknown",
        maf: 0.08, coordinate: "chr21:43060663", consequence: ["Missense"],
        category: "concern", clinvar: "Pathogenic", revel: 0.81, cadd: 27.4, spliceai: 0.02, alphamissense: 0.89,
        qual: 966.2, reads: { matching: 30, total: 58 }, lastEvaluated: "2026-01-09",
        studies: [
          { finding: "Established loss-of-function allele causing classical homocystinuria in the homozygous state", condition: "Homocystinuria", genotypeRelevance: "Your C/T is heterozygous carrier status; unaffected but a carrier", evidenceLevel: 3, source: "Okafor N, et al. Genet Med 2025" }
        ]
      }
    ],
    hpoTerms: [
      { id: "HP:0001939", label: "Abnormality of metabolism/homeostasis", evidence: "Gene-level" },
      { id: "HP:0001083", label: "Ectopia lentis", evidence: "Curated" }
    ],
    publications: [
      { title: "Carrier frequency of CBS pathogenic variants in population biobanks", authors: "Okafor N, et al.", journal: "Genet Med", year: 2025, doi: "10.1016/j.gim.2025.100812", tags: ["Large Cohort", "ClinVar"], finding: "Estimated carrier frequency near 1 in 200 for the most common European pathogenic allele." }
    ]
  },
  {
    symbol: "BRCA1",
    name: "BRCA1 DNA Repair Associated",
    chromosome: "17q21.31",
    organSystem: "Reproductive / Cancer",
    ncbiGeneId: "672",
    omimGene: "113705",
    omimPhenotype: "604370",
    links: refLinks("BRCA1", "672", "113705"),
    summary: "BRCA1 encodes a tumor suppressor involved in homologous recombination DNA repair. Pathogenic variants substantially increase lifetime risk of breast and ovarian cancer.",
    associatedPathology: [
      { name: "Hereditary breast and ovarian cancer syndrome", inheritance: "Autosomal dominant", omim: "604370" }
    ],
    pli: 0.91,
    loeuf: 0.42,
    variantsDetected: 64,
    researchedVariants: 15,
    hpoTermCount: 4,
    goTermCount: 2,
    variants: [
      { id: "rs80357382", genotype: "G/G", zygosity: "Homozygous (ref)", phase: "N/A", maf: 0.0001, coordinate: "chr17:43094692", consequence: ["Frameshift"], category: "uncategorized", clinvar: "Pathogenic (not present in this sample)", revel: null, cadd: null, spliceai: null, alphamissense: null, qual: null, reads: null, lastEvaluated: "2025-09-12", studies: [] },
      { id: "rs1799950", genotype: "A/G", zygosity: "Heterozygous", phase: "Maternal", maf: 0.29, coordinate: "chr17:43093401", consequence: ["Missense"], category: "uncertain", clinvar: "Benign/Likely benign", revel: 0.07, cadd: 8.9, spliceai: 0.00, alphamissense: 0.06, qual: 887.4, reads: { matching: 26, total: 49 }, lastEvaluated: "2024-06-20", studies: [
        { finding: "Common Q356R polymorphism shows no association with breast cancer risk in large case-control studies", condition: "Breast cancer", genotypeRelevance: "Your A/G is heterozygous for a well-studied benign polymorphism", evidenceLevel: 3, source: "CIMBA consortium, 2023" }
      ] }
    ],
    hpoTerms: [
      { id: "HP:0003002", label: "Breast carcinoma", evidence: "Gene-level" },
      { id: "HP:0100615", label: "Ovarian neoplasm", evidence: "Gene-level" }
    ],
    publications: []
  },
  {
    symbol: "CFTR",
    name: "CF Transmembrane Conductance Regulator",
    chromosome: "7q31.2",
    organSystem: "Respiratory",
    ncbiGeneId: "1080",
    omimGene: "602421",
    omimPhenotype: "219700",
    links: refLinks("CFTR", "1080", "602421"),
    summary: "CFTR encodes a chloride channel; loss-of-function variants cause cystic fibrosis. Carrier status is common and clinically significant for reproductive planning.",
    associatedPathology: [{ name: "Cystic fibrosis", inheritance: "Autosomal recessive", omim: "219700" }],
    pli: 0.12,
    loeuf: 1.05,
    variantsDetected: 38,
    researchedVariants: 11,
    hpoTermCount: 5,
    goTermCount: 1,
    variants: [
      { id: "rs113993960", genotype: "G/G", zygosity: "Homozygous (ref)", phase: "N/A", maf: 0.0002, coordinate: "chr7:117559590", consequence: ["Deletion (delF508)"], category: "uncategorized", clinvar: "Pathogenic (not present)", revel: null, cadd: null, spliceai: null, alphamissense: null, qual: null, reads: null, lastEvaluated: "2025-03-04", studies: [] }
    ],
    hpoTerms: [
      { id: "HP:0006538", label: "Recurrent respiratory infections", evidence: "Gene-level" }
    ],
    publications: []
  },
  {
    symbol: "LDLR",
    name: "Low Density Lipoprotein Receptor",
    chromosome: "19p13.2",
    organSystem: "Cardiovascular",
    ncbiGeneId: "3949",
    omimGene: "606945",
    omimPhenotype: "143890",
    links: refLinks("LDLR", "3949", "606945"),
    summary: "LDLR mediates clearance of LDL cholesterol from blood. Pathogenic variants cause familial hypercholesterolemia with markedly elevated LDL and early coronary disease.",
    associatedPathology: [{ name: "Familial hypercholesterolemia", inheritance: "Autosomal dominant", omim: "143890" }],
    pli: 0.88,
    loeuf: 0.55,
    variantsDetected: 52,
    researchedVariants: 13,
    hpoTermCount: 4,
    goTermCount: 2,
    variants: [
      { id: "rs72658867", genotype: "C/T", zygosity: "Heterozygous", phase: "Paternal", maf: 0.02, coordinate: "chr19:11116926", consequence: ["Missense"], category: "concern", clinvar: "Likely pathogenic", revel: 0.76, cadd: 24.1, spliceai: 0.01, alphamissense: 0.72, qual: 951.8, reads: { matching: 31, total: 57 }, lastEvaluated: "2025-11-18", studies: [
        { finding: "Missense variant associated with elevated LDL-C and early-onset coronary artery disease in carrier pedigrees", condition: "Familial hypercholesterolemia", genotypeRelevance: "Your C/T is heterozygous for a likely pathogenic allele", evidenceLevel: 3, source: "Bassiri T, et al. Circ Genom Precis Med 2026" }
      ] }
    ],
    hpoTerms: [
      { id: "HP:0003124", label: "Hypercholesterolemia", evidence: "Gene-level" },
      { id: "HP:0001677", label: "Coronary artery atherosclerosis", evidence: "Curated" }
    ],
    publications: [
      { title: "Genotype-phenotype correlation of LDLR variants in a multi-ancestry FH registry", authors: "Bassiri T, et al.", journal: "Circ Genom Precis Med", year: 2026, doi: "10.1161/CIRCGEN.126.004321", tags: ["Large Cohort", "ClinVar"], finding: "Confirmed strong LDL-C elevation for the pathogenic-class missense variants in this registry." }
    ]
  },
  {
    symbol: "G6PD",
    name: "Glucose-6-Phosphate Dehydrogenase",
    chromosome: "Xq28",
    organSystem: "Hematologic",
    ncbiGeneId: "2539",
    omimGene: "305900",
    omimPhenotype: "300908",
    links: refLinks("G6PD", "2539", "305900"),
    summary: "G6PD protects red blood cells from oxidative damage. Deficiency variants cause hemolytic anemia triggered by oxidative stressors including certain medications and fava beans.",
    associatedPathology: [{ name: "G6PD deficiency", inheritance: "X-linked", omim: "300908" }],
    pli: 0.05,
    loeuf: 1.21,
    variantsDetected: 19,
    researchedVariants: 6,
    hpoTermCount: 3,
    goTermCount: 1,
    variants: [
      { id: "rs1050828", genotype: "C/T", zygosity: "Hemizygous", phase: "Maternal", maf: 0.06, coordinate: "chrX:154536002", consequence: ["Missense"], category: "concern", clinvar: "Pathogenic", revel: 0.62, cadd: 21.6, spliceai: 0.00, alphamissense: 0.64, qual: 902.9, reads: { matching: 40, total: 41 }, lastEvaluated: "2024-09-27", studies: [
        { finding: "A- variant associated with mild-to-moderate G6PD deficiency and drug/fava-bean-triggered hemolysis", condition: "Hemolytic anemia", genotypeRelevance: "Hemizygous carriage of this allele confers deficiency", evidenceLevel: 3, source: "PharmGKB / CPIC curation, 2024" }
      ] }
    ],
    hpoTerms: [
      { id: "HP:0001878", label: "Hemolytic anemia", evidence: "Gene-level" }
    ],
    publications: []
  }
];

// ---------------------------------------------------------------------
// ONTOLOGY TREES (selectable: hpo | go | organ)
// Uniform 3-level shape for all three so Tree / List / Graph share one
// renderer: groups[] -> terms[] -> genes[]
// ---------------------------------------------------------------------
const ONTOLOGIES = {
  hpo: {
    label: "HPO — Human Phenotype Ontology",
    description: "Organ system \u2192 Phenotype terms \u2192 Genes",
    groups: [
      {
        id: "HP:0000118", label: "Skeletal system", genes: ["VDR"],
        terms: [
          { id: "HP:0002748", label: "Rickets", genes: ["VDR"] },
          { id: "HP:0004349", label: "Reduced bone mineral density", genes: ["VDR"] }
        ]
      },
      {
        id: "HP:0001939", label: "Metabolism / homeostasis", genes: ["MTHFR", "AGXT2", "CBS"],
        terms: [
          { id: "HP:0003166", label: "Hyperhomocysteinemia", genes: ["MTHFR", "CBS"] },
          { id: "HP:0003259", label: "Elevated ADMA / renal marker", genes: ["AGXT2"] },
          { id: "HP:0001083", label: "Ectopia lentis", genes: ["CBS"] }
        ]
      },
      {
        id: "HP:0001626", label: "Cardiovascular system", genes: ["MTHFR", "LDLR"],
        terms: [
          { id: "HP:0001677", label: "Coronary artery atherosclerosis", genes: ["MTHFR", "LDLR"] },
          { id: "HP:0003124", label: "Hypercholesterolemia", genes: ["LDLR"] }
        ]
      },
      {
        id: "HP:0000478", label: "Reproductive / Neoplasm", genes: ["BRCA1"],
        terms: [
          { id: "HP:0003002", label: "Breast carcinoma", genes: ["BRCA1"] },
          { id: "HP:0100615", label: "Ovarian neoplasm", genes: ["BRCA1"] }
        ]
      },
      {
        id: "HP:0002086", label: "Respiratory system", genes: ["CFTR"],
        terms: [
          { id: "HP:0006538", label: "Recurrent respiratory infections", genes: ["CFTR"] }
        ]
      },
      {
        id: "HP:0001871", label: "Hematologic system", genes: ["G6PD"],
        terms: [
          { id: "HP:0001878", label: "Hemolytic anemia", genes: ["G6PD"] }
        ]
      },
      {
        id: "HP:0002715", label: "Immune system", genes: ["VDR"],
        terms: [
          { id: "HP:0002960", label: "Autoimmunity", genes: ["VDR"] }
        ]
      }
    ]
  },
  go: {
    label: "GO — Gene Ontology",
    description: "GO category \u2192 GO terms \u2192 Genes",
    groups: [
      {
        id: "GO:0008152", label: "Metabolic process", genes: ["MTHFR", "AGXT2", "CBS", "VDR"],
        terms: [
          { id: "GO:0006520", label: "Cellular amino acid metabolic process", genes: ["AGXT2", "CBS"] },
          { id: "GO:0006730", label: "One-carbon metabolic process", genes: ["MTHFR"] },
          { id: "GO:0070301", label: "Cellular response to hydrogen peroxide", genes: ["VDR"] }
        ]
      },
      {
        id: "GO:0006281", label: "DNA repair", genes: ["BRCA1"],
        terms: [
          { id: "GO:0000724", label: "Double-strand break repair via homologous recombination", genes: ["BRCA1"] },
          { id: "GO:0006302", label: "Double-strand break repair", genes: ["BRCA1"] }
        ]
      },
      {
        id: "GO:0055085", label: "Transmembrane transport", genes: ["CFTR", "LDLR"],
        terms: [
          { id: "GO:0006821", label: "Chloride transport", genes: ["CFTR"] },
          { id: "GO:0034375", label: "High-density lipoprotein remodeling", genes: ["LDLR"] },
          { id: "GO:0034381", label: "Plasma lipoprotein particle clearance", genes: ["LDLR"] }
        ]
      },
      {
        id: "GO:0006979", label: "Response to oxidative stress", genes: ["G6PD"],
        terms: [
          { id: "GO:0006739", label: "NADPH regeneration", genes: ["G6PD"] }
        ]
      },
      {
        id: "GO:0009611", label: "Response to wounding / hormone signaling", genes: ["VDR"],
        terms: [
          { id: "GO:0071305", label: "Cellular response to vitamin D", genes: ["VDR"] }
        ]
      }
    ]
  },
  organ: {
    label: "Organ / System",
    description: "Anatomical system \u2192 Sub-system \u2192 Genes",
    groups: [
      {
        id: "SYS:skeletal", label: "Skeletal", genes: ["VDR"],
        terms: [{ id: "SYS:skeletal:bone", label: "Bone mineralization", genes: ["VDR"] }]
      },
      {
        id: "SYS:metabolic", label: "Metabolic", genes: ["MTHFR", "AGXT2", "CBS"],
        terms: [
          { id: "SYS:metabolic:homocysteine", label: "Homocysteine / transsulfuration", genes: ["MTHFR", "CBS"] },
          { id: "SYS:metabolic:amino", label: "Amino acid / ADMA metabolism", genes: ["AGXT2"] }
        ]
      },
      {
        id: "SYS:cardio", label: "Cardiovascular", genes: ["MTHFR", "LDLR"],
        terms: [
          { id: "SYS:cardio:lipid", label: "Lipid metabolism", genes: ["LDLR"] },
          { id: "SYS:cardio:vascular", label: "Vascular risk factors", genes: ["MTHFR"] }
        ]
      },
      {
        id: "SYS:repro", label: "Reproductive / Cancer", genes: ["BRCA1"],
        terms: [{ id: "SYS:repro:breastovarian", label: "Breast / ovarian tissue", genes: ["BRCA1"] }]
      },
      {
        id: "SYS:resp", label: "Respiratory", genes: ["CFTR"],
        terms: [{ id: "SYS:resp:airway", label: "Airway / exocrine glands", genes: ["CFTR"] }]
      },
      {
        id: "SYS:heme", label: "Hematologic", genes: ["G6PD"],
        terms: [{ id: "SYS:heme:rbc", label: "Red blood cell oxidative defense", genes: ["G6PD"] }]
      }
    ]
  }
};

// ---------------------------------------------------------------------
// ANALYSIS — Polygenic Risk Scores + Pharmacogenomics
// ---------------------------------------------------------------------
const PRS = [
  { trait: "Coronary artery disease", percentile: 82, category: "HIGH", pgsId: "PGS000018", organSystem: "Cardiovascular" },
  { trait: "Type 2 diabetes", percentile: 54, category: "MODERATE", pgsId: "PGS000014", organSystem: "Metabolic" },
  { trait: "Bone mineral density (lumbar spine)", percentile: 21, category: "PROTECTIVE", pgsId: "PGS000123", organSystem: "Skeletal" },
  { trait: "Breast cancer", percentile: 63, category: "MODERATE", pgsId: "PGS000004", organSystem: "Reproductive / Cancer" },
  { trait: "LDL cholesterol", percentile: 91, category: "HIGH", pgsId: "PGS000065", organSystem: "Cardiovascular" }
];

const PGX = [
  { gene: "CYP2C19", diplotype: "*1/*2", phenotype: "Intermediate metabolizer", drug: "Clopidogrel", actionTier: "Actionable", recommendation: "Consider alternative antiplatelet therapy (e.g., prasugrel or ticagrelor) per CPIC guidance for reduced clopidogrel activation." },
  { gene: "SLCO1B1", diplotype: "*1/*5", phenotype: "Intermediate function", drug: "Simvastatin", actionTier: "Actionable", recommendation: "Increased myopathy risk at higher doses; CPIC suggests a lower starting dose or alternative statin." },
  { gene: "CYP2D6", diplotype: "*1/*1", phenotype: "Normal metabolizer", drug: "Codeine", actionTier: "Informational", recommendation: "Standard dosing guidelines apply; no dose adjustment indicated." },
  { gene: "DPYD", diplotype: "*1/*1", phenotype: "Normal metabolizer", drug: "Fluorouracil / Capecitabine", actionTier: "Informational", recommendation: "No evidence of DPD deficiency; standard dosing per protocol." }
];

// ---------------------------------------------------------------------
// REPORT (narrative + breakdown, backs the Reports view)
// ---------------------------------------------------------------------
const REPORT = {
  generated: "2026-08-26",
  sampleLabel: "HG003 (GIAB benchmark, GRCh38, chr1-22)",
  narrative: "This benchmark sample shows one likely-pathogenic finding of note (LDLR, heterozygous, familial hypercholesterolemia association) alongside a CBS carrier variant for homocystinuria and an MTHFR reduced-function allele of common, modest effect. A VDR variant shows a modest reported association with multiple sclerosis risk, balanced by a separate protective VDR association with non-Hodgkin lymphoma. Two actionable pharmacogenomic findings (CYP2C19, SLCO1B1) affect antiplatelet and statin dosing choices. No pathogenic BRCA1 or CFTR alleles from the curated list were detected in this sample.",
  geneBreakdown: GENES.map(g => ({
    symbol: g.symbol,
    variantsDetected: g.variantsDetected,
    pathogenicOrLP: g.variants.filter(v => v.category === "concern" || (v.clinvar || "").toLowerCase().includes("pathogenic")).length,
    protective: g.variants.filter(v => v.category === "protective").length,
    uncertain: g.variants.filter(v => v.category === "uncertain").length
  }))
};

/**
 * CHANGELOG (revision 2 — accuracy pass)
 * - AGXT2 was previously mislabeled with "Hyperoxaluria, primary, type
 *   III" as an associated pathology. That phenotype is caused by HOGA1,
 *   not AGXT2 — removed the incorrect claim; AGXT2 is now described
 *   accurately as an ADMA-metabolism modifier gene with no classic
 *   single-gene OMIM phenotype.
 * - OMIM numbers were previously a mix of gene-entry and
 *   phenotype-entry numbers under one ambiguous "omim" field
 *   (e.g. CFTR used the CF phenotype number 219700 as if it were the
 *   gene entry; LDLR used the FH phenotype number 143890 the same
 *   way). Split into omimGene / omimPhenotype so both are correct and
 *   labeled distinctly: CFTR gene *602421 / CF #219700; LDLR gene
 *   *606945 / FH #143890; CBS gene *613381 / homocystinuria #236200;
 *   BRCA1 gene *113705 / HBOC #604370; MTHFR gene *607093 /
 *   homocystinuria #236250; VDR gene *601769 / rickets #277440.
 * - Added real NCBI Entrez Gene IDs and GeneCards/ClinVar link
 *   patterns for the reference-links row (Overview tab).
 */
