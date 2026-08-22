# PLAN.md: Unified Genomic Ontology & Pharmacogenomics Reporting

## 1. Executive Summary & Core Objectives
The objective is to update `/home/daniel-ehrle/My-Projects/genomics/ontology_report` into a world-class, production-grade genomic and pharmacogenomic reporting engine. 
The system enforces a **domain-agnostic, ontology-driven architecture** where gene panels and actionability criteria are dynamically derived from Human Phenotype Ontology (HPO), MONDO, and EFO mappings.

### Key Architectural Pillars
- **Strict Clinical Segmentation**: High-penetrance monogenic drivers are classified using the standard 3-Tier clinical framework. Polygenic risk factors are **rolled up** using an additive/multiplicative scoring model into trait-level percentiles rather than treated as separate clinical SNPs.
- **Hierarchical HPO Partitioning**: Grouping and filtering are organized by Level 1 (organ systems) and Level 2 (Morphology vs. Physiology) clinical design patterns.
- **Bi-Directional Pharmacogenomics**: Diplotype and star-allele calling (via PharmCAT) mapped to CPIC and DPWG guidelines.
- **FAIR CURIE Standardization**: Every node validated against strict OBO Foundry schema rules using LinkML and Pydantic.

---

## 2. Technical Design & Data Schema
### 2.1 LinkML Schema Specification (`genomic_ontology_schema.yaml`)
We use LinkML to govern data exchange formats. The schema models:
1. `VariantRecord`: Represents raw VCF annotations (rsID, Genotype, MAF, REVEL, Chrom, Pos, Ref, Alt, Phasing).
2. `GeneRecord`: Encapsulates gene metadata (NCBI description, associated HPO terms).
3. `PolygenicRollup`: Represents aggregated risk scores (EFO trait ID, computed score, population percentile, risk categorization).
4. `PharmaRecommendation`: Represents drug-gene interactions (Gene, Diplotype, Phenotype, Drug, Recommendation, Source).

### 2.2 Google Drive Integration via MCP Server
An MCP (Model Context Protocol) Server for Google Drive enables automated, secure, and bidirectional read/write access to the personal Google Drive folder `"Ontology"` (and `"Pharmacogenomics"`).
- **Assumptions**: 
  - OAuth2 Client credentials are stored locally.
  - The MCP server exposes tools: `create_file`, `upload_file`, `list_folder_contents`, and `update_file_content`.
  - Directory IDs for "Ontology" and "Pharmacogenomics" are resolved dynamically.

---

## 3. Multi-Agent Design Review & Collaboration Log

### Agent 1 (Primary Architect - Genomics Specialist)
- **Proposed Design**: Integrate monogenic and polygenic results in a single, high-contrast HTML table. List each polygenic SNP in the table next to monogenic variants.
- **Logic**: Maximizes visibility of all annotated variants in one view.

### Agent 2 (Clinical Genetics Peer Reviewer)
- **Review Critique**: *Disagree with unified table layout.* Polygenic variants (e.g., GWAS SNPs with low odds ratios) should not be mixed with high-penetrance Mendelian variants (e.g., pathogenic *BRCA1* or *POLG* frameshifts). Mixing them creates clinical confusion. Polygenic data is only meaningful when **rolled up** into a multi-locus trait percentile (e.g., 94th percentile risk for CAD).
- **Resolution**: Separate the monogenic and polygenic report components. Monogenic variants are displayed under gene headings with their ACMG/ClinVar/REVEL Tiers. Polygenic risk alleles are mathematically aggregated into an EFO-mapped risk profile shown as horizontal bar charts under the correct HPO Level 1/Level 2 physiological nodes.

### Agent 3 (Technical Systems / Validation Engineer)
- **Review Critique**: Ensure that the phasing data (WhatsHap / SHAPEIT5 maternal vs. paternal haplotype blocks) and short-read quality scores are explicitly rendered. Add strict LinkML-to-Pydantic validation to the pipeline before rendering.
- **Resolution**: Incorporated `models.py` validation inside the build pipeline. Enforced rendering of `maternal`, `paternal`, or `undetermined` phasing states.

### Consensus Agreement (95% Confidence)
The agents have reached **100% agreement** on the split/modular architecture:
- Monogenic reporting uses ClinVar/REVEL variant-level classification.
- Polygenic reporting implements an **additive log-odds algorithm** to roll up SNPs into a trait percentile.
- Data structures are validated using Pydantic classes generated directly from the LinkML schema.
