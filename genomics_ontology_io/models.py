from typing import List, Optional
from pydantic import BaseModel, Field, field_validator
import re

class MonogenicFinding(BaseModel):
    gene_symbol: str = Field(..., description="Canonical gene symbol")
    ncbi_description: Optional[str] = Field(None, description="NCBI Gene description")
    rsid: Optional[str] = Field(None, description="dbSNP rsID")
    chromosome: str = Field(..., description="Chromosome name")
    position: int = Field(..., description="1-based chromosome position")
    genotype: str = Field(..., description="Patient genotype calls (e.g. A/G, -/CAGT)")
    zygosity: str = Field(..., description="Homozygous, Heterozygous, or Hemizygous")
    revel_score: Optional[float] = Field(None, description="REVEL pathogenicity score")
    impact_consequence: str = Field(..., description="Variant sequence consequence (e.g. Missense, Frameshift)")
    clinvar_significance: Optional[str] = Field(None, description="ClinVar clinical significance category")
    phasing: str = Field("undetermined", description="Phased haplotype: maternal, paternal, de_novo, or undetermined")
    associated_hpo_terms: List[str] = Field(default_factory=list, description="OBO-compliant HPO terms (e.g. HP:0001626)")
    associated_mondo_terms: List[str] = Field(default_factory=list, description="OBO-compliant MONDO disease terms")

    @field_validator("associated_hpo_terms")
    @classmethod
    def validate_hpo_terms(cls, values: List[str]) -> List[str]:
        for val in values:
            if not re.match(r"^HP:\d{7}$", val):
                raise ValueError(f"Invalid HPO CURIE format: '{val}'. Must match 'HP:\\d{{7}}'.")
        return values

    @field_validator("associated_mondo_terms")
    @classmethod
    def validate_mondo_terms(cls, values: List[str]) -> List[str]:
        for val in values:
            if not re.match(r"^MONDO:\d{7}$", val):
                raise ValueError(f"Invalid MONDO CURIE format: '{val}'. Must match 'MONDO:\\d{{7}}'.")
        return values


class PolygenicRollup(BaseModel):
    efo_trait_id: str = Field(..., description="Experimental Factor Ontology CURIE (e.g. EFO:0000400)")
    trait_name: str = Field(..., description="Clinical trait name")
    pgs_catalog_id: Optional[str] = Field(None, description="PGS Catalog unique identifier (e.g. PGS000018)")
    computed_score: float = Field(..., description="Raw computed polygenic score")
    percentile: float = Field(..., description="Ancestry-normalized population percentile (0-100)")
    risk_category: str = Field(..., description="Triage classification: HIGH, MODERATE, or PROTECTIVE")
    hpo_level1_system: str = Field(..., description="Level 1 HPO organ system CURIE (e.g. HP:0001626)")
    hpo_level2_subcategory: str = Field(..., description="Level 2 HPO subcategory: Morphology or Physiology")

    @field_validator("efo_trait_id")
    @classmethod
    def validate_efo_id(cls, val: str) -> str:
        if not re.match(r"^EFO:\d{7}$", val):
            raise ValueError(f"Invalid EFO CURIE format: '{val}'. Must match 'EFO:\\d{{7}}'.")
        return val


class PharmaRecommendation(BaseModel):
    gene: str = Field(..., description="PGx pharmacogene name")
    diplotype: str = Field(..., description="Assigned star-allele diplotype (e.g. *1/*4)")
    phenotype: Optional[str] = Field(None, description="Predicted metabolizer phenotype")
    affected_drug: str = Field(..., description="Drug affected by genotype")
    clinical_recommendation: str = Field(..., description="Clinical actionability summary")
    action_tier: str = Field(..., description="STANDARD, FAVOUR, CAUTION, AVOID, DOSE_UP, DOSE_DOWN, MONITOR")
    guideline_source: str = Field("CPIC", description="PGx guideline source (CPIC, DPWG, FDA)")


class VariantReport(BaseModel):
    patient_id: str = Field(..., description="Unique patient identifier")
    run_date: str = Field(..., description="Run timestamp")
    monogenic_findings: List[MonogenicFinding] = Field(default_factory=list)
    polygenic_findings: List[PolygenicRollup] = Field(default_factory=list)
    pharma_findings: List[PharmaRecommendation] = Field(default_factory=list)
