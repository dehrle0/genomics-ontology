#!/usr/bin/env python3
"""
genomics_utils.py
Centralized Genomics & Annotation Utilities for the Ontology Report Pipeline.

Provides shared, robust functions for:
  - String and numeric normalization
  - Variant coordinate serialization and parsing (chr:pos:ref>alt)
  - Zygosity and VAF normalization
  - Quality score (Phred / percentile) formatting
  - Clinical annotation cleaning (ClinVar disease strings, OMIM IDs)
"""

from typing import Any, Optional, Tuple, List

def safe_str(val: Any, default: str = "") -> str:
    """Safely convert a value to string, handling None and missing values."""
    if val is None or str(val).strip().lower() in ("none", "null", "nan", "-"):
        return default
    return str(val).strip()

def safe_num(val: Any, default: Optional[float] = None) -> Optional[float]:
    """Safely convert a value to float, handling None and invalid strings."""
    if val is None or val == "":
        return default
    try:
        f = float(val)
        return default if (f != f) else f  # check for NaN
    except (TypeError, ValueError):
        return default

def format_variant_key(chrom: Any, pos: Any, ref: Any, alt: Any) -> str:
    """Format standard 1-based variant key: chr:pos:ref>alt."""
    c = str(chrom).strip()
    if not c.startswith("chr") and c.upper() not in ("MT", "M"):
        c = f"chr{c}"
    elif c.upper() in ("MT", "M"):
        c = "chrM"
    return f"{c}:{pos}:{ref}>{alt}"

def parse_variant_key(key_str: str) -> Optional[Tuple[str, int, str, str]]:
    """Parse standard 1-based variant key chr:pos:ref>alt into (chrom, pos, ref, alt)."""
    try:
        if ":" not in key_str or ">" not in key_str:
            return None
        chrom, rest = key_str.split(":", 1)
        pos_str, alleles = rest.split(":", 1)
        ref, alt = alleles.split(">", 1)
        return chrom.strip(), int(pos_str.strip()), ref.strip().upper(), alt.strip().upper()
    except (ValueError, IndexError):
        return None

def normalize_zygosity(raw: Any) -> Optional[str]:
    """Normalize varied zygosity representations into a clean label."""
    if raw is None:
        return None
    s = str(raw).strip().lower()
    if not s or s in ("-", "na", "none", "unknown", "."):
        return None
    if s in ("het", "heterozygous", "0/1", "1/0", "0|1", "1|0"):
        return "Heterozygous"
    if s in ("hom", "homozygous", "1/1", "1|1"):
        return "Homozygous"
    if s in ("hemi", "hemizygous", "1", "1/.", "./1"):
        return "Hemizygous"
    if s in ("ref", "0/0", "0|0", "homref"):
        return "Reference"
    return str(raw).strip()

def compute_vaf(vaf: Any, alt_reads: Any, tot_reads: Any) -> Optional[float]:
    """Prefer explicit VAF; otherwise derive from alt_reads / tot_reads."""
    v = safe_num(vaf)
    if v is not None:
        return v
    a = safe_num(alt_reads)
    t = safe_num(tot_reads)
    if a is not None and t not in (None, 0.0):
        return round(a / t, 4)
    return None

def format_phred(score: Any, prefix: str = "Q", precision: int = 1) -> str:
    """Format Phred score string (e.g. Q32.0) or return em-dash if missing."""
    n = safe_num(score)
    if n is None:
        return "—"
    return f"{prefix}{n:.{precision}f}"

def format_percentile(pct: Any, precision: int = 2) -> str:
    """Format genome-wide percentile string (e.g. Top 0.04%)."""
    n = safe_num(pct)
    if n is None:
        return "—"
    return f"Top {n:.{precision}f}%"

def clean_clinvar_disease(disease_str: Any, max_terms: int = 3) -> List[str]:
    """Parse pipe-separated ClinVar diseases into cleaned, unique terms."""
    raw = safe_str(disease_str)
    if not raw:
        return []
    terms = []
    seen = set()
    for t in raw.split("|"):
        cl = t.strip()
        if not cl or cl.lower() in ("not specified", "not provided", "-"):
            continue
        if cl.lower() not in seen:
            seen.add(cl.lower())
            terms.append(cl)
    return terms[:max_terms]

def format_omim_ids(omim_raw: Any) -> str:
    """Normalize OMIM IDs into standard formatted string (e.g. OMIM: 123456; 654321)."""
    raw = safe_str(omim_raw)
    if not raw:
        return ""
    parts = [p.strip().replace("OMIM:", "").replace("MIM:", "").strip() for p in raw.replace(";", ",").split(",")]
    valid = [p for p in parts if p.isdigit()]
    if not valid:
        return ""
    return "OMIM: " + "; ".join(valid)

def get_clinvar_url(clinvar_id: Any) -> Optional[str]:
    """Generate NCBI ClinVar URL for a given variation ID."""
    cid = safe_str(clinvar_id).replace("VCV", "").replace("vcv", "").strip()
    if cid and cid.isdigit():
        return f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{cid}/"
    return None

def get_dbsnp_url(rsid: Any) -> Optional[str]:
    """Generate NCBI dbSNP URL for a given rsID."""
    r = safe_str(rsid).lower().strip()
    if r.startswith("rs") and r[2:].isdigit():
        return f"https://www.ncbi.nlm.nih.gov/snp/{r}"
    return None

def get_omim_url(omim_id: Any) -> Optional[str]:
    """Generate OMIM entry URL for the primary OMIM accession."""
    raw = safe_str(omim_id)
    if not raw:
        return None
    parts = [p.strip().replace("OMIM:", "").replace("MIM:", "").strip() for p in raw.replace(";", ",").split(",")]
    digits = [p for p in parts if p.isdigit()]
    if digits:
        return f"https://www.omim.org/entry/{digits[0]}"
    return None

def get_clingen_url(hugo: Any) -> Optional[str]:
    """Generate ClinGen gene validity URL for a HUGO gene symbol."""
    h = safe_str(hugo).strip().upper()
    if h:
        return f"https://search.clinicalgenome.org/kb/genes/{h}"
    return None

def get_alphagenome_url(chrom: Any, pos: Any, ref: Any, alt: Any) -> Optional[str]:
    """Generate DeepMind AlphaGenome Atlas variant exploration URL."""
    c = safe_str(chrom).strip()
    if not c.startswith("chr") and c:
        c = f"chr{c}"
    p = safe_str(pos).strip()
    r = safe_str(ref).strip().upper()
    a = safe_str(alt).strip().upper()
    if c and p and r and a and r != "-" and a != "-":
        return f"https://alphagenome.deepmind.com/variant/{c}:{p}:{r}>{a}"
    elif c:
        return f"https://alphagenome.deepmind.com/locus/{c}:{p}-{p}"
    return None

def get_pubmed_url(pmid: Any) -> Optional[str]:
    """Generate NCBI PubMed URL for a PMID."""
    p = safe_str(pmid).replace("PMID:", "").replace("PMID", "").strip()
    if p and p.isdigit():
        return f"https://pubmed.ncbi.nlm.nih.gov/{p}/"
    return None

def get_gnomad_url(chrom: Any, pos: Any, ref: Any, alt: Any) -> Optional[str]:
    """Generate Broad Institute gnomAD v4 URL for variant coordinates."""
    c = safe_str(chrom).strip().replace("chr", "")
    p = safe_str(pos).strip()
    r = safe_str(ref).strip().upper()
    a = safe_str(alt).strip().upper()
    if c and p and r and a:
        return f"https://gnomad.broadinstitute.org/variant/{c}-{p}-{r}-{a}?dataset=gnomad_r4"
    return None

