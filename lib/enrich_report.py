#!/usr/bin/env python3
"""
enrich_report.py
Enrich an actionable-variant JSON (from ontology_filter.py) with information
pulled live from public genomics web services, then write the enriched JSON back
for the renderers to consume.

Two independent enrichments (either or both, chosen by flags):

  --genes    NCBI Gene description + summary + cytogenetic location per gene
             (NCBI E-utilities). This is the "NCBI Gene description" field
             requested for every report (cardiology, hereditary cancer, ...).

  --studies  Current GWAS study associations for each variant / gene from the
             EBI GWAS Catalog REST API — trait, p-value, odds ratio, risk
             allele, and the underlying PubMed publication. This powers the
             autoimmunity report's "current info on study results" feature.

Design goals:
  * Offline-safe. --offline (or any network failure) degrades gracefully: the
    report is still produced, just without the enriched fields.
  * Cached. Every remote answer is memoised to a JSON cache file so re-runs and
    repeated genes/variants cost nothing and stay reproducible.
  * Polite. Requests are rate-limited and carry a descriptive User-Agent.

The engine stays domain-agnostic: which enrichments run is decided by the
orchestrator / config, not hard-coded per disease.
"""
import argparse
import json
import os
import sys
import time
import urllib.parse
import urllib.request

USER_AGENT = "ontology_report/1.0 (+research; genomics domain pipeline)"
NCBI_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
GWAS_API = "https://www.ebi.ac.uk/gwas/rest/api"

# Cache schema version — bump if the stored record shape changes.
CACHE_VERSION = 1


class Fetcher:
    """Small cached, rate-limited, offline-aware HTTP-JSON client."""

    def __init__(self, cache_path, offline=False, min_interval=0.34, timeout=15,
                 ncbi_api_key=None):
        self.cache_path = cache_path
        self.offline = offline
        self.min_interval = min_interval
        self.timeout = timeout
        self.ncbi_api_key = ncbi_api_key or os.environ.get("NCBI_API_KEY")
        self._last = 0.0
        self.calls = 0
        self.errors = 0
        self.cache = {"_version": CACHE_VERSION}
        if cache_path and os.path.exists(cache_path):
            try:
                loaded = json.load(open(cache_path))
                if loaded.get("_version") == CACHE_VERSION:
                    self.cache = loaded
            except (OSError, ValueError):
                self.cache = {"_version": CACHE_VERSION}

    def save(self):
        if not self.cache_path:
            return
        tmp = self.cache_path + ".tmp"
        with open(tmp, "w") as f:
            json.dump(self.cache, f)
        os.replace(tmp, self.cache_path)

    def _throttle(self):
        dt = time.time() - self._last
        if dt < self.min_interval:
            time.sleep(self.min_interval - dt)
        self._last = time.time()

    def get_json(self, url, cache_key):
        """Return parsed JSON for url, using the cache. Cache stores {"ok":..}
        or {"err":..}; both short-circuit future calls within a run/re-run."""
        if cache_key in self.cache:
            entry = self.cache[cache_key]
            return entry.get("ok")
        if self.offline:
            return None
        self._throttle()
        try:
            req = urllib.request.Request(
                url, headers={"User-Agent": USER_AGENT, "Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=self.timeout) as r:
                data = json.loads(r.read())
            self.calls += 1
            self.cache[cache_key] = {"ok": data, "ts": int(time.time())}
            return data
        except Exception as e:  # network, HTTP, JSON — all non-fatal
            self.errors += 1
            self.cache[cache_key] = {"err": str(e)[:200], "ts": int(time.time())}
            return None

    def _ncbi_key(self, url):
        if self.ncbi_api_key:
            sep = "&" if "?" in url else "?"
            return f"{url}{sep}api_key={self.ncbi_api_key}"
        return url


# --------------------------------------------------------------------------- #
# NCBI Gene description
# --------------------------------------------------------------------------- #
def gene_info(fetcher, symbol):
    """Return {ncbi_gene_id, description, summary, map_location, aliases} for a
    HUGO symbol, or None. Two-step E-utilities: esearch -> esummary."""
    if not symbol:
        return None
    ck = f"gene:{symbol}"
    if ck in fetcher.cache and "ok" in fetcher.cache[ck]:
        return fetcher.cache[ck]["ok"]
    if fetcher.offline and ck not in fetcher.cache:
        return None

    term = urllib.parse.quote(f"{symbol}[sym] AND Homo sapiens[orgn]")
    search_url = fetcher._ncbi_key(
        f"{NCBI_EUTILS}/esearch.fcgi?db=gene&term={term}&retmode=json&retmax=1")
    sd = fetcher.get_json(search_url, f"gene_search:{symbol}")
    ids = (sd or {}).get("esearchresult", {}).get("idlist", []) if sd else []
    if not ids:
        fetcher.cache[ck] = {"ok": None, "ts": int(time.time())}
        return None
    gid = ids[0]
    sum_url = fetcher._ncbi_key(
        f"{NCBI_EUTILS}/esummary.fcgi?db=gene&id={gid}&retmode=json")
    su = fetcher.get_json(sum_url, f"gene_summary:{gid}")
    g = (su or {}).get("result", {}).get(gid) if su else None
    if not g:
        fetcher.cache[ck] = {"ok": None, "ts": int(time.time())}
        return None
    info = {
        "ncbi_gene_id": gid,
        "description": g.get("description") or g.get("nomenclaturename"),
        "summary": (g.get("summary") or "").strip() or None,
        "map_location": g.get("maplocation") or None,
        "aliases": g.get("otheraliases") or None,
    }
    fetcher.cache[ck] = {"ok": info, "ts": int(time.time())}
    return info


# --------------------------------------------------------------------------- #
# GWAS Catalog study evidence
# --------------------------------------------------------------------------- #
def _assoc_pvalue(a):
    """GWAS Catalog reports very small p-values as mantissa x 10^exponent; the
    plain `pvalue` field underflows to 0.0. Reconstruct from the parts so tiny
    p-values (e.g. 3e-198) survive."""
    m = a.get("pvalueMantissa")
    e = a.get("pvalueExponent")
    if m not in (None, "") and e not in (None, ""):
        try:
            return f"{float(m)}e{int(e)}"
        except (TypeError, ValueError):
            pass
    return a.get("pvalue")


def _summarize_assoc(a):
    traits = [t.get("trait") for t in (a.get("efoTraits") or []) if t.get("trait")]
    href = a.get("_links", {}).get("study", {}).get("href")
    return {
        "traits": traits,
        "pvalue": _assoc_pvalue(a),
        "or_beta": a.get("orPerCopyNum") or a.get("betaNum"),
        "risk_allele": _risk_allele(a),
        "risk_freq": a.get("riskFrequency"),
        "_study_href": href,
        "pubmed": None,
    }


def _risk_allele(a):
    for locus in (a.get("loci") or []):
        for ra in (locus.get("strongestRiskAlleles") or []):
            if ra.get("riskAlleleName"):
                return ra["riskAlleleName"]
    return None


def _study_pub(fetcher, href, max_hops=1):
    """Follow an association's study link to grab the PubMed citation."""
    if not href or fetcher.offline:
        return None
    ck = f"study:{href}"
    s = fetcher.get_json(href, ck)
    pi = (s or {}).get("publicationInfo") if s else None
    if not pi:
        return None
    return {
        "pubmed_id": pi.get("pubmedId"),
        "title": pi.get("title"),
        "journal": pi.get("publication"),
        "date": pi.get("publicationDate"),
        "author": (pi.get("author") or {}).get("fullname"),
    }


def snp_studies(fetcher, rsid, max_studies=6, with_pub=True):
    """Return current GWAS associations for a dbSNP rsID, trait-ranked by
    p-value (smallest first)."""
    if not rsid or not str(rsid).startswith("rs"):
        return None
    ck = f"gwas_snp:{rsid}"
    if ck in fetcher.cache and "ok" in fetcher.cache[ck]:
        return fetcher.cache[ck]["ok"]
    if fetcher.offline and ck not in fetcher.cache:
        return None
    url = (f"{GWAS_API}/singleNucleotidePolymorphisms/{rsid}/associations"
           f"?projection=associationBySnp")
    d = fetcher.get_json(url, f"gwas_snp_raw:{rsid}")
    assoc = (d or {}).get("_embedded", {}).get("associations", []) if d else []
    out = []
    for a in assoc:
        s = _summarize_assoc(a)
        if s["traits"]:
            out.append(s)

    def _pv(x):
        try:
            return float(x["pvalue"])
        except (TypeError, ValueError, KeyError):
            return 1.0
    out.sort(key=_pv)
    out = out[:max_studies]
    if with_pub:
        for s in out:
            s["pubmed"] = _study_pub(fetcher, s.pop("_study_href", None))
    else:
        for s in out:
            s.pop("_study_href", None)
    result = {"rsid": rsid, "n_associations": len(assoc), "top": out}
    fetcher.cache[ck] = {"ok": result, "ts": int(time.time())}
    return result


def gene_trait_summary(fetcher, symbol, max_traits=8):
    """Aggregate distinct EFO traits reported for a gene across the GWAS Catalog
    (gene-level view, cheap: one request, no per-study follow)."""
    if not symbol:
        return None
    ck = f"gwas_gene:{symbol}"
    if ck in fetcher.cache and "ok" in fetcher.cache[ck]:
        return fetcher.cache[ck]["ok"]
    if fetcher.offline and ck not in fetcher.cache:
        return None
    url = (f"{GWAS_API}/singleNucleotidePolymorphisms/search/findByGene"
           f"?geneName={urllib.parse.quote(symbol)}")
    d = fetcher.get_json(url, f"gwas_gene_raw:{symbol}")
    snps = (d or {}).get("_embedded", {}).get(
        "singleNucleotidePolymorphisms", []) if d else []
    rsids = sorted({s.get("rsId") for s in snps if s.get("rsId")})
    result = {"gene": symbol, "n_gwas_snps": len(rsids), "rsids": rsids[:25]}
    fetcher.cache[ck] = {"ok": result, "ts": int(time.time())}
    return result


# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #
def enrich(data, cache_path, do_genes=True, do_studies=False, offline=False,
           max_variant_studies=6):
    fetcher = Fetcher(cache_path, offline=offline)
    records = data.get("records", [])

    gene_cache = {}          # symbol -> gene_info (per-run memo, dedup work)
    genes_enriched = 0
    studies_enriched = 0

    for rec in records:
        sym = rec.get("hugo")
        if do_genes and sym:
            if sym not in gene_cache:
                gene_cache[sym] = gene_info(fetcher, sym)
            gi = gene_cache[sym]
            if gi:
                rec["gene_info"] = gi
                genes_enriched += 1
        if do_studies:
            se = {}
            snp = snp_studies(fetcher, rec.get("rsid"),
                              max_studies=max_variant_studies)
            if se_ok(snp):
                se["snp"] = snp
            gt = gene_trait_summary(fetcher, sym)
            if gt and gt.get("n_gwas_snps"):
                se["gene"] = gt
            if se:
                rec["study_evidence"] = se
                studies_enriched += 1

    fetcher.save()
    meta = {
        "offline": offline,
        "genes_requested": do_genes,
        "studies_requested": do_studies,
        "records_with_gene_info": genes_enriched,
        "records_with_study_evidence": studies_enriched,
        "remote_calls": fetcher.calls,
        "remote_errors": fetcher.errors,
        "cache": cache_path,
    }
    data["enrichment"] = meta
    return data, meta


def se_ok(snp):
    return bool(snp and snp.get("top"))


def main():
    ap = argparse.ArgumentParser(description="Enrich actionable report JSON with "
                                             "NCBI gene descriptions + GWAS studies")
    ap.add_argument("--in-json", required=True)
    ap.add_argument("--out-json", help="defaults to overwriting --in-json")
    ap.add_argument("--cache", help="JSON cache file (defaults next to in-json)")
    ap.add_argument("--genes", action="store_true", help="fetch NCBI gene descriptions")
    ap.add_argument("--studies", action="store_true", help="fetch GWAS study evidence")
    ap.add_argument("--offline", action="store_true",
                    help="never hit the network; use cache only")
    ap.add_argument("--max-variant-studies", type=int, default=6)
    args = ap.parse_args()

    data = json.load(open(args.in_json))
    out_json = args.out_json or args.in_json
    cache = args.cache or (os.path.splitext(args.in_json)[0] + "_enrich_cache.json")
    do_genes = args.genes or not (args.genes or args.studies)  # default: genes on

    data, meta = enrich(data, cache, do_genes=do_genes, do_studies=args.studies,
                        offline=args.offline,
                        max_variant_studies=args.max_variant_studies)
    with open(out_json, "w") as f:
        json.dump(data, f, indent=2)

    print(f"[enrich] offline={meta['offline']} calls={meta['remote_calls']} "
          f"errors={meta['remote_errors']}")
    print(f"[enrich] gene descriptions on {meta['records_with_gene_info']} records; "
          f"study evidence on {meta['records_with_study_evidence']} records")
    print(f"[enrich] wrote {out_json} (cache {cache})")


if __name__ == "__main__":
    main()
