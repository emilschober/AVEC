import os
import sys
import traceback
from flask import Flask, request, jsonify, send_from_directory, abort
import time
import re
import requests
import threading
from collections import OrderedDict
from concurrent.futures import ThreadPoolExecutor
from decimal import Decimal, getcontext
from typing import Dict, Any, Optional, Tuple, List
from Bio.Seq import Seq
import pandas as pd
from html import escape
from urllib.parse import urlencode, urljoin, quote

# ===== CONFIGURATION =====
_BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DOCS_DIR = os.path.abspath(os.path.join(_BASE_DIR, '../docs'))
ENSEMBL_REST = "https://rest.ensembl.org"
# The new Ensembl platform moved from beta.ensembl.org to www.ensembl.org (the old host now
# 301-redirects). Keep these on the live host: a redirected POST is downgraded to a GET by
# requests, which the GraphQL endpoint answers with its HTML explorer page instead of JSON.
ENSEMBL_GRAPHQL = "https://www.ensembl.org/data/graphql"
ENSEMBL_REFGET = "https://www.ensembl.org/data/refget"
GRCH38_ASSEMBLY_ACCESSION = "GCA_000001405.29"
HEADERS = {"Content-Type": "application/json", "Accept": "application/json"}

NCBI_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
# Optional NCBI E-utilities API key (raises NCBI's rate limit from 3 to 10 requests/s) and contact
# e-mail, which NCBI asks registered tools to send; both are read from the environment.
NCBI_API_KEY = #insert API-Key Here to accelerate the ensembl requests
NCBI_EMAIL = os.environ.get("NCBI_EMAIL") or None
NCBI_TOOL = "AVEC"
# ClinVar esummary returns at most 500 records per request; larger regions are fetched page by
# page, a few pages at a time (see _fetch_clinvar_region).
CLINVAR_PAGE_SIZE = 500
_NCBI_CONCURRENCY = threading.Semaphore(3)
# Exon variant lookups are widened by this many intronic bases on each side so the exon's own
# canonical splice sites (+1/+2 donor, -1/-2 acceptor) are included in the splice-variant counts.
SPLICE_FLANK_BP = 2
LITVAR2_SENSOR = "https://www.ncbi.nlm.nih.gov/research/litvar2-api/sensor"
LITVAR2_DOCSUM = "https://www.ncbi.nlm.nih.gov/research/litvar2/docsum"
MIN_REVIEW_STARS = 2
# Longest accepted variant query. Real HGVS input is far shorter; the cap also keeps
# parse_hgvs_query's regexes (quadratic on long input without a colon) cheap.
MAX_QUERY_LENGTH = 200
# Only these file types are served from DOCS_DIR under /static/ (logo, figures); everything else
# there, notably docs/data with the underlying databases, stays private.
STATIC_ASSET_EXTENSIONS = {'.png', '.jpg', '.jpeg', '.svg', '.ico', '.webp'}
# Figures of the N1C VARIANT guidelines page (docs/guidelines.html), served under /guidelines/images/
GUIDELINES_IMAGES_DIR = os.path.join(DOCS_DIR, 'images-N1CVARIANTGuidelines2.0')

TREXOME_API_BASE = "https://tandem-repeat-domain-database.onrender.com"
TREXOME_FRONTEND_BASE = "https://trexomedb.rnatherapy.nl"

REVIEW_STATUS_TO_STARS = {
    "no assertion provided": 0,
    "no assertion criteria provided": 0,
    "no classification provided": 0,
    "no classifications from unflagged records": 0,
    "criteria provided, single submitter": 1,
    "criteria provided, conflicting classifications": 1,
    "criteria provided, conflicting interpretations of pathogenicity": 1,
    "criteria provided, multiple submitters, no conflicts": 2,
    "reviewed by expert panel": 3,
    "practice guideline": 4,
}
N1C_API_URL = "https://gene-registry.onrender.com/api/data?table=N1C_projects" 
N1C_API_ASSESSED_URL = "https://gene-registry.onrender.com/api/data?table=assessed_variants"
# The N1C tables are re-fetched in the background once a day, and every 10 minutes while
# either one is unavailable (see _refresh_n1c_if_stale).
N1C_REFRESH_SECONDS = 24 * 3600
N1C_RETRY_SECONDS = 10 * 60
N1C_UNAVAILABLE_NOTE = "N1C registry currently unavailable; existing-ASO matches may be missing."
_N1C_REFRESH_LOCK = threading.Lock()
_n1c_last_attempt: float = 0.0
# Global DataFrames to be loaded at startup
clingen_df: Optional[pd.DataFrame] = None
splicevar_df: Optional[pd.DataFrame] = None
n1c_variants_df: Optional[pd.DataFrame] = None 
n1c_assessed_df: Optional[pd.DataFrame] = None 
sscvdb_df: Optional[pd.DataFrame] = None 
rcnv_df: Optional[pd.DataFrame] = None

# Supplementary N1C table with gene-level features (uORF, NAT, PE, NMD-inducing exons, etc.)
n1c_supp_df: Optional[pd.DataFrame] = None

# Columns surfaced from the N1C supplementary table for WT-upregulation assessment
N1C_SUPP_DISPLAY_COLUMNS = [
    'HUGO Antisense', 'uORF', 'NAT', 'PE',
    'A3_NMD', 'A5_NMD', 'AI_NMD', 'EI_NMD', 'ES_NMD', 'NMD_in',
    'exon_id', 'UKBB_variants_at_exon_splice_sites',
    'high_confidence_uORF_initiates_in_exon',
    'high_confidence_close_uORF_present_after_skipping',
    'high_confidence_ouORF_present_after_skipping',
    'high_confidence_close_uORF_removed_by_skipping',
    'high_confidence_ouORF_removed_by_skipping',
]

# --- Data Loading ---
def load_databases():
    """
    Loads all necessary data files and fetches N1C registry data
    into global pandas DataFrames.
    """
    DATA_DIR = os.path.join(DOCS_DIR, 'data')

    clingen_path = os.path.join(DATA_DIR, 'Clingen-Curation-Activity-Summary-Report-2026-09-21.csv')
    splicevar_path = os.path.join(DATA_DIR, 'SpliceVarDbHgvs.xlsx')
    sscvdb_path = os.path.join(DATA_DIR, 'SSCVDB_hgvs.xlsx')
    n1c_supp_path = os.path.join(DATA_DIR, 'Supp File for Combined TANGO Approaches (May2026).xlsx')

    global clingen_df, splicevar_df, sscvdb_df, n1c_supp_df, rcnv_df
    try:
        clingen_df = pd.read_csv(clingen_path).set_index('gene_symbol')

        # Load SpliceVarDB from Excel
        splicevar_df = pd.read_excel(splicevar_path)
        # Sanitize SpliceVarDB data (critical for lookups)
        splicevar_df.columns = splicevar_df.columns.str.strip()
        for col in ['hgvs', 'gene', 'hg38', 'hgvs_genomic']:
            if col in splicevar_df.columns:
                splicevar_df[col] = splicevar_df[col].astype(str).str.strip()
        # Upper-cased once here so lookups can compare directly against _build_hg38_key_from_vep
        if 'hg38' in splicevar_df.columns:
            splicevar_df['hg38'] = splicevar_df['hg38'].str.upper()

        # Load SSCVDB from Excel
        try:
            sscvdb_df = pd.read_excel(sscvdb_path)
            sscvdb_df.columns = sscvdb_df.columns.str.strip()
            if 'Variant_ID' in sscvdb_df.columns:
                sscvdb_df['Variant_ID'] = sscvdb_df['Variant_ID'].astype(str).str.strip().str.upper()
        except Exception as e:
            print(f"Warning: Could not load SSCVDB.xlsx: {e}")

        # Load N1C supplementary table (uORF / NAT / PE / NMD-inducing exons per gene).
        # 'Supp Table' has a paper-grouping header row above the real column names, hence header=1.
        try:
            n1c_supp_df = pd.read_excel(n1c_supp_path, sheet_name='Supp Table', header=1)
            n1c_supp_df.columns = n1c_supp_df.columns.str.strip()
            if 'Gene' in n1c_supp_df.columns:
                n1c_supp_df['Gene'] = n1c_supp_df['Gene'].astype(str).str.strip()

        except Exception as e:
            print(f"Warning: Could not load N1C supplementary table: {e}")

        _load_n1c_tables()

        rcnv_path = os.path.join(DATA_DIR, 'rCNV.gene_scores.tsv')
        if os.path.exists(rcnv_path):
            rcnv_df = pd.read_csv(rcnv_path, sep='\t')
            # Standardize gene column for matching
            if 'gene' in rcnv_df.columns:
                rcnv_df['gene'] = rcnv_df['gene'].astype(str).str.strip().str.upper()
        else:
            print(f"Warning: rCNV file not found at {rcnv_path}")

    except Exception as e:
        traceback.print_exc()
        print(f"An error occurred during database loading: {e}")
        sys.exit(1)

def _fetch_n1c_variants() -> pd.DataFrame:
    """Fetches the N1C projects registry table."""
    response = requests.get(N1C_API_URL, timeout=120)
    response.raise_for_status()
    df = pd.DataFrame(response.json())
    # Sanitize the N1C columns we will search on
    if 'Gene' in df.columns:
        df['Gene'] = df['Gene'].astype(str).str.strip()
    if 'Coding DNA change (c.)' in df.columns:
        df['Coding DNA change (c.)'] = df['Coding DNA change (c.)'].astype(str).str.strip()
    return df

def _fetch_n1c_assessed() -> pd.DataFrame:
    """Fetches the N1C assessed variants (curated) table."""
    response = requests.get(N1C_API_ASSESSED_URL, timeout=120)
    response.raise_for_status()
    df = pd.DataFrame(response.json())
    # Normalize every string-named column to stripped strings; fillna first so missing
    # values become '' rather than the literal string 'nan'
    for col in list(df.columns):
        if isinstance(col, str):
            df[col] = df[col].fillna('').astype(str).str.strip()
    return df

def _load_n1c_tables() -> None:
    """(Re-)fetches both N1C registry tables. The registry is an external service (cold starts,
    outages), so a failed fetch only disables the N1C checks instead of taking the whole app
    down (downstream lookups treat a None DataFrame as "no match"), and a failed refresh keeps
    the previously loaded table rather than dropping it."""
    global n1c_variants_df, n1c_assessed_df, _n1c_last_attempt
    _n1c_last_attempt = time.monotonic()
    try:
        n1c_variants_df = _fetch_n1c_variants()
    except Exception as e:
        print(f"Warning: Could not load N1C registry: {e}")
    try:
        n1c_assessed_df = _fetch_n1c_assessed()
    except Exception as e:
        print(f"Warning: Could not load N1C assessed variants: {e}")

def _n1c_unavailable() -> bool:
    return n1c_variants_df is None or n1c_assessed_df is None

def _refresh_n1c_if_stale() -> None:
    """Starts a background re-fetch of the N1C tables when they are older than
    N1C_REFRESH_SECONDS, or N1C_RETRY_SECONDS after the last attempt while either is
    unavailable. Never blocks the calling request, which keeps using the current tables."""
    max_age = N1C_RETRY_SECONDS if _n1c_unavailable() else N1C_REFRESH_SECONDS
    if time.monotonic() - _n1c_last_attempt < max_age:
        return
    if not _N1C_REFRESH_LOCK.acquire(blocking=False):
        return  # a refresh is already running
    def run():
        try:
            _load_n1c_tables()
        finally:
            _N1C_REFRESH_LOCK.release()
    threading.Thread(target=run, daemon=True).start()

_COMPLEMENT_TABLE = str.maketrans('ACGTN', 'TGCAN')

def _reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT_TABLE)[::-1]

def _left_align_indel(client: 'EnsemblClient', chrom: str, pos: int, ref: str, alt: str, max_shift: int = 50):
    """Shifts a single-base-anchor indel (ref or alt of length 1) to its leftmost
    equivalent position within a repeat/homopolymer, matching the left-aligned
    VCF convention used by varaico/ClinVar (e.g. pos=X ref=C alt=CC -> pos=X-1 ref=T alt=TC)."""
    if len(ref) == len(alt):
        return pos, ref, alt
    if len(ref) < len(alt):
        anchor, indel_seq, is_insertion = ref, alt[1:], True
    else:
        anchor, indel_seq, is_insertion = alt, ref[1:], False
    if len(anchor) != 1:
        return pos, ref, alt
    window_start = max(1, pos - max_shift)
    if window_start >= pos:
        return pos, ref, alt
    window_seq = client.get_sequence_region(chrom, window_start, pos - 1)
    if not window_seq or len(window_seq) != pos - window_start:
        return pos, ref, alt
    window_seq = window_seq.upper()
    i = len(window_seq) - 1
    while i >= 0 and indel_seq[-1] == anchor:
        anchor, indel_seq = window_seq[i], anchor + indel_seq[:-1]
        pos -= 1
        i -= 1
    if is_insertion:
        return pos, anchor, anchor + indel_seq
    return pos, anchor + indel_seq, anchor

def _vep_to_vcf(vep_entry: Dict[str, Any], client: Optional['EnsemblClient'] = None, anchor_indels: bool = True) -> Optional[Tuple[str, int, str, str]]:
    """Converts a VEP entry to a forward-strand, left-aligned VCF-style (chrom, pos, ref, alt),
    with chrom unprefixed (e.g. '17'). This is the shared form behind the SpliceVarDB/SSCVDB
    lookup keys, the varaico link and the genomic notation.
    VEP's start/seq_region_name are always forward-strand, but for HGVS input allele_string is
    reported on the transcript's strand (BRCA1 c.5074G>A comes back as 'G/A', strand -1, while
    the genome and SpliceVarDB have C>T), so strand -1 alleles are reverse-complemented.
    VEP also reports indels without an anchor base ('A/-', '-/A'); with anchor_indels the
    preceding reference base is fetched (e.g. ref=TA, alt=T) and the indel left-aligned.
    Returns None for missing/ambiguous alleles, or for indels when anchor_indels is False."""
    chrom = str(vep_entry.get('seq_region_name') or '').strip()
    if chrom.lower().startswith('chr'):
        chrom = chrom[3:]
    pos = vep_entry.get('start')
    allele = vep_entry.get('allele_string') or ''
    if not chrom or not pos or '/' not in allele:
        return None
    pos = int(pos)
    ref, alt = allele.split('/', 1)
    if vep_entry.get('strand') == -1:
        ref = ref if ref == '-' else _reverse_complement(ref)
        alt = alt if alt == '-' else _reverse_complement(alt)
    if ref == '-' or alt == '-':
        if not anchor_indels:
            return None
        def fetch():
            seq_client = client or EnsemblClient()
            anchor_pos = pos - 1
            anchor_base = seq_client.get_sequence_region(chrom, anchor_pos, anchor_pos)
            anchor_base = anchor_base.upper() if anchor_base else None
            if not anchor_base or not re.fullmatch(r'[ACGT]', anchor_base):
                return None
            a_ref = anchor_base if ref == '-' else anchor_base + ref
            a_alt = anchor_base if alt == '-' else anchor_base + alt
            return _left_align_indel(seq_client, chrom, anchor_pos, a_ref, a_alt)
        # Several callers convert the same variant within one assessment; cache so the
        # reference-sequence fetches for anchoring/left-alignment only happen once.
        anchored = _ENSEMBL_LOOKUP_CACHE.get_or_set(('vep_to_vcf', chrom, pos, ref, alt), fetch)
        if not anchored:
            return None
        pos, ref, alt = anchored
    if not re.fullmatch(r'[ACGTN]+', ref) or not re.fullmatch(r'[ACGTN]+', alt):
        return None
    return chrom, pos, ref, alt

def _build_hg38_key_from_vep(vep_entry: Dict[str, Any], client: Optional['EnsemblClient'] = None) -> Optional[str]:
    """SpliceVarDB 'hg38' lookup key: <chrom>-<pos>-<ref>-<alt> (e.g. '17-43045803-C-T')."""
    vcf = _vep_to_vcf(vep_entry, client)
    return "{}-{}-{}-{}".format(*vcf).upper() if vcf else None

def _format_sscvdb_variant_id_from_vep(vep_entry: Dict[str, Any], client: Optional['EnsemblClient'] = None) -> Optional[str]:
    """SSCVDB 'Variant_ID' lookup key: chr<chrom>-<pos>-<ref>-<alt> (e.g. 'chr1-100007090-T-A')."""
    vcf = _vep_to_vcf(vep_entry, client)
    return "CHR{}-{}-{}-{}".format(*vcf).upper() if vcf else None

def _build_varaico_url(gene_symbol: Optional[str], vep_entry: Dict[str, Any], client: Optional['EnsemblClient'] = None) -> Optional[str]:
    """Builds a varaico.com paper-search URL for this variant (GRCh38/hg38), using the
    forward-strand, anchored VCF form from _vep_to_vcf. Returns None if the gene symbol,
    coordinates or alleles are unavailable."""
    try:
        if not gene_symbol:
            return None
        vcf = _vep_to_vcf(vep_entry, client)
        if not vcf:
            return None
        chrom, pos, ref, alt = vcf
        params = {"gene": gene_symbol, "chrom": f"chr{chrom}", "pos": pos, "ref": ref, "alt": alt, "assembly": "hg38"}
        return f"https://varaico.com/papers?{urlencode(params)}"
    except Exception:
        return None

def _extract_rsid(vep_entry: Dict[str, Any]) -> Optional[str]:
    """Returns the dbSNP rsID VEP reports as co-located with this variant, or None.
    VEP allele-matches its co-located variants, but the same position can also carry HGMD,
    COSMIC and somatic entries, so only germline `rs` identifiers are considered. dbSNP
    anchors indels one base to the left of VEP's start (e.g. VEP 117559592-117559594 CTT/-
    vs rs113993960 at 117559591-117559594), so coordinates are not required to match: an
    exact positional match wins, otherwise the first germline rsID reported is used."""
    start, end = vep_entry.get('start'), vep_entry.get('end')
    fallback = None
    for colocated in vep_entry.get('colocated_variants') or []:
        rsid = str(colocated.get('id') or '').strip()
        if not re.fullmatch(r'rs\d+', rsid) or colocated.get('somatic'):
            continue
        if colocated.get('start') == start and colocated.get('end') == end:
            return rsid
        if fallback is None:
            fallback = rsid
    return fallback

def _litvar_docsum_url(rsid: str) -> str:
    """The LitVar2 website URL for an rsID, built locally instead of read off the sensor
    response (it uses the same construction), so the link survives a slow, rate-limited or
    unreachable LitVar2 the same way the varaico link does."""
    return f"{LITVAR2_DOCSUM}?{urlencode({'variant': f'litvar@{rsid}##', 'query': rsid})}"

def _litvar_sensor(rsid: str) -> Optional[Dict[str, Any]]:
    """Asks LitVar2's sensor endpoint for the publication count of an rsID.
    Returns None when LitVar2 knows no publications for the rsID (the endpoint answers 400
    'Variant not found') or the service is unreachable. Deliberately kept on a short timeout
    with no retries: the count only decorates an otherwise locally built link, so a slow
    LitVar2 must neither stall an assessment nor cost us the link itself. Cached via the
    TTL cache, which (unlike lru_cache) doesn't pin a failed lookup's None."""
    def fetch():
        _t0 = time.perf_counter()
        try:
            resp = requests.get(f"{LITVAR2_SENSOR}/{rsid}", headers=HEADERS, timeout=10)
            if resp.status_code != 200:
                return None
            data = resp.json()
            return data if isinstance(data, dict) and data.get('link') else None
        except (requests.RequestException, ValueError):
            return None
        finally:
            _log_timing("LitVar2: sensor", _t0)
    return _ENSEMBL_LOOKUP_CACHE.get_or_set(('litvar_sensor', rsid), fetch)

def _build_litvar_info(vep_entry: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    """Builds the LitVar2 link block for this variant (link, rsID, publication count).
    Returns None only when the variant has no dbSNP rsID; the sensor merely contributes the
    publication count, so an unavailable LitVar2 no longer suppresses the link entirely."""
    rsid = _extract_rsid(vep_entry)
    if not rsid:
        return None
    sensor = _litvar_sensor(rsid) or {}
    return {
        "rsid": rsid,
        "url": sensor.get("link") or _litvar_docsum_url(rsid),
        "pmids_count": sensor.get("pmids_count"),
    }

_AA_3TO1 = {
    'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C', 'Gln': 'Q', 'Glu': 'E',
    'Gly': 'G', 'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F',
    'Pro': 'P', 'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V', 'Ter': '*',
    'Sec': 'U', 'Pyl': 'O',
}

def _build_variant_notations(target_consequence: Dict[str, Any], protein_effect: Optional[str], vep_entry: Dict[str, Any]) -> List[str]:
    """
    Builds every HGVS-style notation this variant is known by (cDNA, single-letter
    protein, genomic), regardless of which one the user originally typed, so
    variant-specific literature searches (Pathomechanism/Splicing) can match papers
    that cite the variant in any convention.
    """
    notations = []

    hgvsc = target_consequence.get('hgvsc') or ''
    cdna_match = re.search(r'(c\.[^\s\)\]]+)', hgvsc, re.IGNORECASE)
    if cdna_match:
        notations.append(cdna_match.group(1))

    protein_match = re.search(r'p\.\(?([A-Za-z]{3})(\d+)([A-Za-z]{3})\)?', protein_effect or '')
    if protein_match:
        aa_from = _AA_3TO1.get(protein_match.group(1).capitalize())
        aa_to = _AA_3TO1.get(protein_match.group(3).capitalize())
        if aa_from and aa_to:
            pos_num = protein_match.group(2)
            notations.append(f"p.{aa_from}{pos_num}{aa_to}")
            notations.append(f"{aa_from}{pos_num}{aa_to}")

    vcf = _vep_to_vcf(vep_entry, anchor_indels=False)
    if vcf and len(vcf[2]) == 1 and len(vcf[3]) == 1:
        _, pos, ref, alt = vcf
        notations.append(f"g.{pos}{ref}>{alt}")

    return notations

# --- Refget sequence access (Ensembl beta platform) ---
# GraphQL resolves a sequence's refget checksum (chromosome regions, transcript CDS, etc.);
# refget then serves the actual bases for that checksum. Neither service supports looking
# a sequence up by chromosome name or transcript ID directly, so GraphQL is the resolver step.

def _graphql_query(query: str, variables: Optional[Dict[str, Any]] = None, max_retries: int = 3, label: str = "Ensembl GraphQL") -> Optional[Dict[str, Any]]:
    """POSTs a query to Ensembl's GraphQL service and returns the 'data' payload, or None on failure."""
    url = ENSEMBL_GRAPHQL
    backoff = 1.0
    _t0 = time.perf_counter()
    try:
        for attempt in range(max_retries):
            try:
                # Redirects are followed by hand: requests would turn a redirected POST into a body-less
                # GET, which this endpoint answers with its GraphiQL HTML page (HTTP 200) rather than JSON.
                resp = requests.post(url, json={'query': query, 'variables': variables or {}}, timeout=120, allow_redirects=False)
                if resp.is_redirect and resp.headers.get('Location'):
                    url = urljoin(url, resp.headers['Location'])
                    continue
                if resp.status_code == 200:
                    if 'json' not in resp.headers.get('Content-Type', ''):
                        print(f"[WARN] {label}: non-JSON response ({resp.headers.get('Content-Type')}) from {url}", flush=True)
                        return None
                    return resp.json().get('data')
                if resp.status_code not in (429, 500, 502, 503, 504) or attempt + 1 == max_retries:
                    return None
            except (requests.RequestException, ValueError):
                if attempt + 1 == max_retries:
                    return None
            time.sleep(backoff)
            backoff *= 2
        return None
    finally:
        _log_timing(label, _t0)

# These resolvers cache through _ENSEMBL_LOOKUP_CACHE (defined further down; resolved at call
# time) rather than functools.lru_cache, because lru_cache would memoize a None from one failed
# lookup - and a pinned-None genome id silently disables every CDS/region fetch until restart.

def _grch38_genome_id() -> Optional[str]:
    """Resolves and caches the genome_id of the current GRCh38 release on Ensembl's new platform."""
    def fetch():
        data = _graphql_query(
            'query($acc: String) { genomes(by_keyword: { assembly_accession_id: $acc }) { genome_id release_number } }',
            {'acc': GRCH38_ASSEMBLY_ACCESSION},
            label="Ensembl GraphQL: genome id"
        )
        genomes = data.get('genomes') if data else None
        if not genomes:
            return None
        return max(genomes, key=lambda g: g['release_number'])['genome_id']
    return _ENSEMBL_LOOKUP_CACHE.get_or_set(('grch38_genome_id',), fetch)

def _region_checksum(chrom: str) -> Optional[str]:
    """Resolves and caches a chromosome/region's refget checksum via GraphQL."""
    def fetch():
        genome_id = _grch38_genome_id()
        if not genome_id:
            return None
        data = _graphql_query(
            'query($gid: String!, $name: String!) { region(by_name: { genome_id: $gid, name: $name }) { sequence { checksum } } }',
            {'gid': genome_id, 'name': chrom},
            label="Ensembl GraphQL: region checksum"
        )
        region = data.get('region') if data else None
        return region['sequence']['checksum'] if region else None
    return _ENSEMBL_LOOKUP_CACHE.get_or_set(('region_checksum', chrom), fetch)

def _cds_checksum(transcript_id: str) -> Optional[str]:
    """Resolves and caches a transcript's CDS refget checksum via GraphQL."""
    def fetch():
        genome_id = _grch38_genome_id()
        if not genome_id:
            return None
        data = _graphql_query(
            'query($gid: String!, $sid: String!) { transcript(by_id: { genome_id: $gid, stable_id: $sid }) '
            '{ product_generating_contexts { cds { sequence { checksum } } } } }',
            {'gid': genome_id, 'sid': transcript_id},
            label="Ensembl GraphQL: CDS checksum"
        )
        transcript = data.get('transcript') if data else None
        if not transcript:
            return None
        for ctx in transcript.get('product_generating_contexts', []):
            cds = ctx.get('cds')
            checksum = cds.get('sequence', {}).get('checksum') if cds else None
            if checksum:
                return checksum
        return None
    return _ENSEMBL_LOOKUP_CACHE.get_or_set(('cds_checksum', transcript_id), fetch)

def _refget_sequence(checksum: str, start: Optional[int] = None, end: Optional[int] = None, max_retries: int = 3) -> Optional[str]:
    """Fetches a sequence (or 0-based half-open subsequence) from Ensembl refget by checksum."""
    if not checksum:
        return None
    params = {}
    if start is not None: params['start'] = start
    if end is not None: params['end'] = end
    backoff = 1.0
    _t0 = time.perf_counter()
    try:
        for attempt in range(max_retries):
            try:
                resp = requests.get(
                    f"{ENSEMBL_REFGET}/sequence/{checksum}", headers={'Accept': 'text/plain'}, params=params, timeout=120
                )
                if resp.status_code == 200:
                    return resp.text
                if resp.status_code not in (429, 500, 502, 503, 504) or attempt + 1 == max_retries:
                    return None
            except requests.RequestException:
                if attempt + 1 == max_retries:
                    return None
            time.sleep(backoff)
            backoff *= 2
        return None
    finally:
        _log_timing("Ensembl refget: sequence fetch", _t0)

class _TTLCache:
    """Thread-safe, size-bounded, time-expiring cache, shared by every EnsemblClient
    instance (i.e. across requests) so that re-assessing a variant after only a MOI/MoA/
    guideline override changed doesn't re-fetch Ensembl reference data (VEP consequences,
    transcript/exon structure, protein domains, overlapping genes, ...) that didn't change.

    In-memory only (nothing is ever written to disk) and never holds user-entered content
    (overrides, commentary, the raw query string) - only Ensembl's own public reference
    annotations, keyed by variant/transcript/protein/gene identifier. Entries expire after
    `ttl_seconds` and the oldest-used entry is evicted once `maxsize` is exceeded, so nothing
    accumulates in the process indefinitely.
    """
    def __init__(self, maxsize: int = 2000, ttl_seconds: float = 24 * 3600):
        self._maxsize = maxsize
        self._ttl = ttl_seconds
        self._lock = threading.Lock()
        self._data: "OrderedDict[Any, Tuple[float, Any]]" = OrderedDict()

    def get_or_set(self, key: Any, fetch):
        now = time.monotonic()
        with self._lock:
            entry = self._data.get(key)
            if entry is not None:
                ts, value = entry
                if now - ts < self._ttl:
                    self._data.move_to_end(key)
                    return value
                del self._data[key]
        # Fetched outside the lock so one slow network call can't block unrelated cache
        # lookups; a concurrent duplicate fetch of the same not-yet-cached key is possible
        # but harmless (last write wins, no data is corrupted).
        value = fetch()
        # A None result means the lookup failed (network blip, 4xx/5xx, service moved); don't
        # memoize it, or one failure would be served as the answer for the rest of the TTL.
        if value is None:
            return None
        with self._lock:
            self._data[key] = (time.monotonic(), value)
            self._data.move_to_end(key)
            while len(self._data) > self._maxsize:
                self._data.popitem(last=False)
        return value

# Shared across requests (see _TTLCache docstring above) - the sole cache backing
# EnsemblClient._cached, so identical Ensembl lookups made by different requests (e.g. a
# reassessment triggered by a MOI/pathomechanism override) are served from memory instead
# of hitting Ensembl again.
_ENSEMBL_LOOKUP_CACHE = _TTLCache(maxsize=2000, ttl_seconds=24 * 3600)
# TReXome only annotates notes, and its host (onrender) can be slow or asleep. Its lookups
# cache "no data" and failures alike (as empty results), for one hour, so an outage costs
# one slow request per gene per hour instead of stalling every assessment on retries.
_TREXOME_CACHE = _TTLCache(maxsize=500, ttl_seconds=3600)

# Shared across requests: pure connection-pooling/rate-limiting infrastructure. Holds
# open sockets and an in-flight-request counter only.
_ENSEMBL_SESSION = requests.Session()
_ENSEMBL_SESSION.headers.update(HEADERS)
# Bounds how many Ensembl REST calls this process has in flight at once, so parallelizing
# independent lookups within one assessment can't turn into a burst that trips Ensembl's
# own rate limiting (which the existing 429 backoff below would then have to absorb).
_ENSEMBL_CONCURRENCY = threading.Semaphore(6)
# Shared pool for parallelizing independent leaf-level I/O within a single assessment
# (e.g. per-exon variant/ClinVar lookups, per-N1C-row VEP lookups). Tasks submitted here
# must not themselves submit to _IO_EXECUTOR and block on the result - see _PAIR_EXECUTOR.
_IO_EXECUTOR = ThreadPoolExecutor(max_workers=16)
# assess_consecutive_exons evaluates its (up to 2) candidate exon pairs concurrently, and
# each pair's evaluation itself submits leaf tasks to _IO_EXECUTOR and waits on them. Using
# a separate pool for this outer level keeps the dependency one-directional (_PAIR_EXECUTOR
# waits on _IO_EXECUTOR, never the reverse), which avoids the classic thread-pool deadlock
# where every worker in a shared pool ends up blocked waiting for a free worker in that same
# pool.
_PAIR_EXECUTOR = ThreadPoolExecutor(max_workers=4)
# TReXome lookups are started at the beginning of an assessment and run alongside it (see
# _assess_for_transcript). Its host (onrender) can take tens of seconds to wake up, so it gets
# its own pool: a slow TReXome must not occupy the _IO_EXECUTOR workers other lookups need.
_TREXOME_EXECUTOR = ThreadPoolExecutor(max_workers=16)
# Fetches the pages of a large ClinVar region lookup side by side. The lookup itself may run on
# _IO_EXECUTOR, so the pages get their own pool (they submit nothing further); NCBI's rate limit
# is respected by _NCBI_CONCURRENCY inside _ncbi_request.
_NCBI_EXECUTOR = ThreadPoolExecutor(max_workers=3)

def _log_timing(label: str, start: float) -> None:
    """Prints how long a step took. `label` must identify the kind of call only
    (endpoint/purpose), never per-request content (gene symbol, HGVS, IDs) - this is
    stdout-only diagnostic output, not a persisted log, but keep it that way regardless."""
    print(f"[TIMING] {label}: {(time.perf_counter() - start) * 1000:.1f} ms", flush=True)

class EnsemblClient:
    def __init__(self, base_url=ENSEMBL_REST, delay=0.1):
        self.base_url = base_url.rstrip('/')
        self.session = _ENSEMBL_SESSION
        self.delay = delay

    def _cached(self, key, fetch):
        # Backed by the shared, TTL-bounded _ENSEMBL_LOOKUP_CACHE (see its docstring) -
        # a fresh EnsemblClient is created per request, but the underlying cache persists
        # across requests, so an override-only reassessment reuses previously-fetched
        # Ensembl reference data instead of re-fetching it.
        return _ENSEMBL_LOOKUP_CACHE.get_or_set(key, fetch)

    def _get(self, path, params=None, max_retries=5, label=None):
        url = f"{self.base_url}{path}"
        backoff = 1.0
        _t0 = time.perf_counter()
        try:
            for attempt in range(max_retries):
                time.sleep(self.delay)
                try:
                    with _ENSEMBL_CONCURRENCY:
                        resp = self.session.get(url, params=params, timeout=120)
                    if resp.status_code == 200:
                        try: return resp.json()
                        except ValueError: return resp.text
                    elif resp.status_code in (429, 503):
                        wait = float(resp.headers.get('Retry-After', backoff))
                        time.sleep(wait); backoff *= 2
                    elif 500 <= resp.status_code < 600:
                        time.sleep(backoff); backoff *= 2
                    else:
                        if 400 <= resp.status_code < 500: return None
                except requests.RequestException:
                    if attempt + 1 == max_retries: raise
                    time.sleep(backoff); backoff *= 2
            return None
        finally:
            _log_timing(label or "Ensembl REST call", _t0)

    # Stays on REST: GraphQL's CDS.start/end snap to the spanning exon's outer boundary rather than the
    # true start/stop codon position (verified against BRCA1/ENST00000357654 by checking the genomic bases
    # at each boundary for the expected reverse-complemented start codon). Using them would silently shift
    # extract_exons_from_transcript's per-exon CDS-overlap math, which feeds exon-skipping frame/stop-codon checks.
    def lookup_id_expand(self, identifier):
        return self._cached(('lookup_id_expand', identifier), lambda: self._get(f"/lookup/id/{identifier}", params={'expand': '1'}, label="Ensembl: transcript lookup (expand)"))
    # Stays on REST: Ensembl's beta GraphQL schema has no variant/consequence/VEP type at all, so there is
    # nothing to migrate this to.
    def vep_hgvs(self, hgvs_string):
        hgvs_string = hgvs_string.strip()
        # Percent-encode user input so '/', '?', '#' etc. can't alter the Ensembl request path
        encoded = quote(hgvs_string, safe=":>+-_.()*")
        return self._cached(('vep_hgvs', hgvs_string), lambda: self._get(
            f"/vep/human/hgvs/{encoded}/", params={'hgvs': 1, 'variant_class': 1, 'mane': 1, 'canonical': 1, 'xref_refseq': 1},
            label="Ensembl: VEP lookup"
        ))
    def get_protein_meta(self, protein_id):
        """Fetches /lookup/id metadata (used for protein length) for a translation ID."""
        return self._cached(('protein_meta', protein_id), lambda: self._get(f"/lookup/id/{protein_id}", label="Ensembl: protein metadata"))
    def get_cds_sequence(self, transcript_id):
        """Fetches a transcript's CDS sequence via GraphQL (checksum) + refget (bases)."""
        def fetch():
            checksum = _cds_checksum(transcript_id)
            return _refget_sequence(checksum) if checksum else None
        return self._cached(('cds_sequence', transcript_id), fetch)
    def get_sequence_region(self, chrom, start, end):
        """Fetches a 1-based inclusive genomic region via GraphQL (checksum) + refget (bases).
        Refget uses 0-based half-open coordinates, so the start is shifted by one; the end is
        numerically identical (1-based inclusive end == 0-based exclusive end)."""
        checksum = _region_checksum(chrom)
        return _refget_sequence(checksum, start - 1, end) if checksum else None
    def get_domains(self, protein_id):
        return self._cached(('domains', protein_id), lambda: self._fetch_domains(protein_id))
    def _fetch_domains(self, protein_id):
        all_features = self._get(f"/overlap/translation/{protein_id}", params={"feature": "protein_feature"}, label="Ensembl: protein domains")
        if not all_features or not isinstance(all_features, list): return []
        allowed_sources = {src.lower() for src in (
            'CDD','Pfam','SMART','PROSITE profiles','DisProt','PROSITE patterns','PRINTS','ProDom', 'mobidb-lite', 'RepeatsDB'
        )}
        def _domain_url(feat: Dict[str, Any]) -> Optional[str]:
            src = str(feat.get('type', '')).strip().lower()
            identifier = str(feat.get('id') or feat.get('hit_id') or '').strip()
            interpro_id = str(feat.get('interpro') or '').strip()
            ident_upper, interpro_upper = identifier.upper(), interpro_id.upper()
            if src == 'pfam' and identifier:
                return f"https://www.ebi.ac.uk/interpro/entry/pfam/{ident_upper}/"
            if src == 'smart' and identifier:
                return f"https://www.ebi.ac.uk/interpro/entry/smart/{ident_upper}/"
            if src.startswith('prosite') and identifier:
                return f"https://prosite.expasy.org/{ident_upper}"
            if src == 'prints' and identifier:
                return f"https://www.ebi.ac.uk/interpro/entry/prints/{ident_upper}/"
            if src == 'cdd' and identifier:
                return f"https://www.ncbi.nlm.nih.gov/Structure/cdd/cddsrv.cgi?uid={identifier}"
            if src == 'superfamily' and identifier:
                return f"https://www.ebi.ac.uk/interpro/entry/superfamily/{ident_upper}/"
            if src == 'disprot' and identifier:
                return f"https://disprot.org/{identifier}"
            if interpro_upper.startswith('IPR'):
                return f"https://www.ebi.ac.uk/interpro/entry/InterPro/{interpro_upper}/"
            return None
        def _display_name(feat: Dict[str, Any]) -> str:
            desc = feat.get('description') or feat.get('id') or feat.get('hit_id') or 'Domain'
            src = feat.get('type') or feat.get('source')
            return f"{desc} ({src})" if src else str(desc)
        def _is_domain_feature(feature):
            ftype = str(feature.get('type', '')).strip()
            desc = str(feature.get('description', '') or '').lower()
            return (
                ftype.lower() in allowed_sources
                or any(k in desc for k in ('disorder','repeat','coiled-coil','coiled coil','low complexity'))
            )
        preliminary_domains = []
        for f in all_features:
            if not _is_domain_feature(f): continue
            feat = dict(f)
            feat['source'] = feat.get('type')
            feat['url'] = _domain_url(feat)
            feat['display_name'] = _display_name(feat)
            preliminary_domains.append(feat)
        unique_domains = {}
        for feat in preliminary_domains:
            src_norm = str(feat.get('source') or feat.get('type') or '').lower()
            ftype_norm = str(feat.get('type') or '').lower()
            desc_norm = str(feat.get('description') or '').lower()
            id_norm = str(feat.get('id') or feat.get('hit_id') or '').lower()
            logic_norm = str(feat.get('logic_name') or '').lower()
            is_mobidb = any('mobidb' in val for val in (src_norm, ftype_norm, desc_norm, id_norm, logic_norm))
            is_disorder_like = any('disorder' in val for val in (ftype_norm, desc_norm, logic_norm, id_norm))
            is_mobidb_disorder = is_mobidb and is_disorder_like
            if is_mobidb_disorder:
                feat['source'] = 'MobiDBLite'
                feat['display_name'] = 'MobiDBLite disorder prediction'
                feat['id'] = feat.get('id') or feat.get('hit_id') or 'MobiDBLite'
            key = None
            if is_mobidb_disorder:
                key = ('mobidb-lite', 'disorder_prediction')
            else:
                key = feat.get('interpro') or (feat.get('type'), feat.get('start'), feat.get('end'), feat.get('description') or feat.get('id'))
            if is_mobidb_disorder and key in unique_domains:
                existing = unique_domains[key]
                merged = dict(existing)
                merged['start'] = min(existing.get('start', float('inf')), feat.get('start', float('inf')))
                merged['end'] = max(existing.get('end', float('-inf')), feat.get('end', float('-inf')))
                merged['description'] = existing.get('description') or feat.get('description') or "MobiDBLite disorder prediction"
                merged['display_name'] = existing.get('display_name') or "MobiDBLite disorder prediction"
                merged['id'] = existing.get('id') or feat.get('id') or feat.get('hit_id')
                unique_domains[key] = merged
                continue
            if key not in unique_domains: unique_domains[key] = feat
        return list(unique_domains.values())
    # Stays on REST: Ensembl's beta GraphQL schema has no variant/variation type at all, so ClinVar/dbSNP
    # overlap (used for star-rating filtering) has nothing to migrate to.
    def overlap_region_variation(self, chrom, start, end):
        def fetch():
            data = self._get(f"/overlap/region/human/{chrom}:{start}-{end}", params={'feature': 'variation'}, label="Ensembl: region variants")
            return data if isinstance(data, list) else []
        return self._cached(('overlap_region_variation', chrom, start, end), fetch)
    # Stays on REST: GraphQL's overlap_region(by_slice) returned a different gene set than this endpoint
    # for the same region in a live comparison (e.g. it dropped SLX4IP from JAG1's overlap set) with no
    # clear explanation, so it isn't a verified drop-in replacement for the WT-upregulation/NAT search.
    def get_overlapping_genes(self, gene_id):
        """Fetches all genes that overlap with a given Ensembl Gene ID."""
        def fetch():
            data = self._get(f"/overlap/id/{gene_id}", params={"feature": "gene"}, label="Ensembl: overlapping genes")
            return data if isinstance(data, list) else []
        return self._cached(('overlapping_genes', gene_id), fetch)
    def lookup_symbol(self, symbol):
        """Fetches gene data for a given symbol via GraphQL (REST-shaped result; no exon/CDS expansion,
        so the unverified GraphQL CDS-boundary issue affecting lookup_symbol_expand doesn't apply here)."""
        def fetch():
            genome_id = _grch38_genome_id()
            if not genome_id:
                return None
            data = _graphql_query(
                'query($gid: String!, $sym: String!) { genes(by_symbol: { genome_id: $gid, symbol: $sym }) '
                '{ unversioned_stable_id symbol metadata { biotype { value } } '
                'slice { region { name } location { start end } strand { value } } } }',
                {'gid': genome_id, 'sym': symbol}
            )
            genes = data.get('genes') if data else None
            if not genes:
                return None
            gene = genes[0]
            slice_ = gene.get('slice') or {}
            location = slice_.get('location') or {}
            return {
                'id': gene.get('unversioned_stable_id'),
                'display_name': gene.get('symbol'),
                'biotype': ((gene.get('metadata') or {}).get('biotype') or {}).get('value'),
                'seq_region_name': (slice_.get('region') or {}).get('name'),
                'start': location.get('start'),
                'end': location.get('end'),
                'strand': (slice_.get('strand') or {}).get('value'),
            }
        return self._cached(('lookup_symbol', symbol), fetch)

    # Stays on REST: same exon/CDS-boundary imprecision as lookup_id_expand (see comment above it) applies
    # here too, since this also expands to Transcript/Exon/Translation for extract_exons_from_transcript.
    def lookup_symbol_expand(self, symbol):
        """Fetches expanded gene data for a given symbol."""
        def fetch():
            data = self._get(f"/lookup/symbol/human/{quote(symbol, safe='-_.')}", params={'expand': '1'}, label="Ensembl: gene lookup (expand)")
            return data if isinstance(data, dict) else None
        return self._cached(('lookup_symbol_expand', symbol), fetch)
    
# --- Helper & Parsing Functions ---

def _evaluate_splice_variant_position(
    variant_hgvs: str,
    vep_data: Dict[str, Any],
    details: Dict[str, Any],
    exonic_pathogenic_user_input: Optional[str] = None,
    moa_user_input: Optional[str] = None,
    overrides: Optional[Dict[str, bool]] = None,
    target_consequence: Optional[Dict[str, Any]] = None
) -> Optional[Dict[str, Any]]:
    """
    Assesses a validated splice-altering variant based on its genomic position.
    This logic is shared by database-found variants and user-validated variants.
    """
    overrides = overrides or {}
    if overrides.get('ASO already exists', False):
        return {
            "classification": "Eligible",
            "reason": "User confirms an ASO for this variant already exists, therefore it is eligible for splice switching.",
            "checks": {
                "ASO already exists": {"passed": True, "original_passed": False, "overridden": True}
            }
        }

    client = EnsemblClient()
    # Find the most relevant consequence to get gene and transcript, unless one was already resolved by the caller
    all_consequences = vep_data.get('transcript_consequences', [])
    if target_consequence is None:
        target_consequence = choose_best_consequence(all_consequences)
    if not target_consequence:
        return {"classification": "Unable to Assess", "reason": "Could not determine a target transcript from VEP response."}

    gene_symbol = target_consequence.get('gene_symbol')
    transcript_id = target_consequence.get('transcript_id')
    if not gene_symbol or not transcript_id:
        return {"classification": "Unable to Assess", "reason": "Could not determine gene symbol or transcript ID."}

    core_hgvs_match = re.search(r'(c\..*)', variant_hgvs, re.IGNORECASE)
    if not core_hgvs_match:
        # Fall back to the cDNA HGVS from the VEP transcript consequence (handles genomic-format inputs)
        hgvsc = target_consequence.get('hgvsc', '')
        core_hgvs_match = re.search(r'(c\..*)', hgvsc, re.IGNORECASE) if hgvsc else None
    if not core_hgvs_match:
        return {"classification": "Unable to Assess", "reason": "Could not parse HGVS notation for splice-site position assessment.", "details": details}
    core_canonical_hgvs = core_hgvs_match.group(1).lower()

    result = {"details": details}
    consequence_terms = set(target_consequence.get('consequence_terms', []))

    is_intronic_by_consequence = 'intron_variant' in consequence_terms or 'splice_acceptor_variant' in consequence_terms or 'splice_donor_variant' in consequence_terms
    # An intronic offset is a sign directly after a position number (c.123+5, c.124-2, c.-50+3);
    # a leading minus alone (c.-15A>G) is a 5'UTR position, not an intronic offset.
    offset_match = re.search(r'\d([+-])(\d+)', core_canonical_hgvs)
    is_intronic = is_intronic_by_consequence or offset_match is not None

    if is_intronic:
        if not offset_match:
            return {"classification": "Unable to Assess", "reason": "Could not determine distance to canonical splice site from HGVS notation.", "details": details}
        offset_sign, dist = offset_match.group(1), int(offset_match.group(2))
        result['details']['Distance to Canonical Splice Site'] = f"{dist} bp"

        bp_check_name = "Variant does not weaken branchpoint"

        if offset_sign == '+': # Downstream (e.g., c.123+1G>A)
            if dist <= 5:
                result.update({"classification": "Not Eligible", "reason": "Variant is a splice-altering variant located too close to the canonical splice site (<=+5bp). If there is functional evidence variant results in leaky exon skipping, consider exon inclusion."})
            else:
                original_passed = dist > 50
                has_override = bp_check_name in overrides
                final_passed = overrides.get(bp_check_name, original_passed)
                bp_checks = {bp_check_name: {"passed": final_passed, "original_passed": original_passed, "overridden": has_override}}
                if final_passed:
                    result.update({"classification": "Likely Eligible", "reason": "Variant is a splice-altering variant in a favorable deep-intronic position (>+50bp).", "checks": bp_checks})
                else:
                    result.update({"classification": "Unlikely Eligible", "reason": "Variant is a splice-altering variant located near the canonical splice site (+6-+50bp).", "checks": bp_checks})

        else: # Upstream (e.g., c.124-2A>G)
            if dist <= 5:
                result.update({"classification": "Not Eligible", "reason": "Variant is a splice-altering variant located too close to the canonical splice site (>=-5bp). If there is functional evidence variant results in leaky exon skipping, consider exon inclusion."})
            else:
                original_passed = dist > 100
                has_override = bp_check_name in overrides
                final_passed = overrides.get(bp_check_name, original_passed)
                bp_checks = {bp_check_name: {"passed": final_passed, "original_passed": original_passed, "overridden": has_override}}
                if final_passed:
                    result.update({"classification": "Likely Eligible", "reason": "Variant is a splice-altering variant in a favorable deep-intronic position (<-100bp).", "checks": bp_checks})
                else:
                    result.update({"classification": "Unlikely Eligible", "reason": "Variant is a splice-altering variant located near the canonical splice site (-6 to -100bp).", "checks": bp_checks})
        
    elif any(c in consequence_terms for c in ['missense_variant', 'synonymous_variant', 'stop_gained', 'frameshift_variant', 'inframe_deletion', 'inframe_insertion']):
        # This is an exonic variant.
        if exonic_pathogenic_user_input is None:
            # Prompt the user for input on pathogenicity
            return {
                "classification": "Awaiting User Input",
                "reason": "This is an exonic splice-altering variant. To assess for splice correction, please confirm if the variant is pathogenic *without* its splicing effect.",
                "user_validation_prompt_exonic_pathogenicity": True,
                "moa_user_input": moa_user_input
            }
        elif exonic_pathogenic_user_input == 'yes':
            result.update({
                "classification": "Manual Assessment Required",
                "reason": "User confirmed the variant is pathogenic independent of its splicing effect. Manual lookup should be performed. Please refer to the about/methods splice correction page for alternative strategies."
            })
        else: # User answered 'no'
            # Assess based on distance to exon boundaries.
            gene_data = client.lookup_symbol_expand(gene_symbol)
            if not gene_data or 'Transcript' not in gene_data:
                 return {"classification": "Unable to Assess", "reason": f"Could not fetch expanded gene data for {gene_symbol}."}

            transcript_object = next((t for t in gene_data['Transcript'] if t.get('id') == transcript_id), None)
            if not transcript_object:
                return {"classification": "Unable to Assess", "reason": f"Could not find transcript {transcript_id} in expanded gene data."}
            
            all_exons = extract_exons_from_transcript(transcript_object)
            
            variant_pos = vep_data.get('start')

            if not variant_pos:
                return {"classification": "Unable to Assess", "reason": "Could not determine variant position for distance calculation."}

            all_distances = []
            for exon in all_exons:
                if exon.get('start') and exon.get('end'):
                    all_distances.append(abs(variant_pos - exon['start']))
                    all_distances.append(abs(variant_pos - exon['end']))
            
            if not all_distances:
                return {"classification": "Unable to Assess", "reason": "Could not calculate distances to any exon boundaries."}

            dist = min(all_distances)
            
            result['details']['Distance to Canonical Splice Site'] = f"{dist} bp"

            if dist <= 5:
                result.update({"classification": "Not Eligible", "reason": "User confirmed benign protein effect, but variant is within 5bp of an exon boundary. This may be eligible for exon inclusion if leaky exon skipping is observed."})
            elif 6 <= dist <= 15:
                result.update({"classification": "Unlikely Eligible", "reason": "User confirmed benign protein effect. Variant is between 6-15bp of an exon boundary."})
            else: # dist > 15
                result.update({"classification": "Likely Eligible", "reason": "User confirmed benign protein effect. Variant is >15bp from an exon boundary, making it a potential candidate for splice correction."})

    else:
        result.update({"classification": "Unable to Assess", "reason": "The variant's consequence type is neither intronic nor a recognized exonic change, so its position relative to the splice site could not be assessed."})

    return result

def _split_transcript_and_gene(identifier: str) -> Tuple[str, Optional[str]]:
    """
    Splits HGVS identifiers of the form 'NM_12345.6(GENE)' into the bare transcript ID
    and the parenthesized gene symbol. Returns (identifier, None) when there's no suffix.
    """
    match = re.match(r'^([A-Za-z0-9_.]+)\(([A-Za-z0-9\-]+)\)$', identifier)
    if match:
        return match.group(1), match.group(2)
    return identifier, None

def parse_hgvs_query(query: str) -> Tuple[Optional[str], Optional[str], Optional[str]]:
    """
    Parses a query string into a VEP-compatible HGVS string, an optional gene symbol,
    and an optional transcript ID explicitly typed by the user.
    """
    query = query.strip()

    # Pattern 1: Handles formats like "NM_12345.6:c.123A>G", "NM_12345.6(GENE):c.123A>G", or "GENE:c.123A>G"
    match_colon = re.search(r'([^:]+):([cgnmp]\..*)', query, re.IGNORECASE)
    if match_colon:
        identifier = match_colon.group(1).strip()
        variant = match_colon.group(2).strip()
        hgvs_string = f"{identifier}:{variant}"
        # If the identifier is a transcript, we don't have a gene symbol from the query,
        # unless one was given in parentheses (e.g. "NM_12345.6(GENE)")
        if identifier.startswith("NM_") or identifier.startswith("ENST"):
            transcript_id, gene_from_parens = _split_transcript_and_gene(identifier)
            return hgvs_string, gene_from_parens, transcript_id
        else:
            return hgvs_string, identifier, None

    # Pattern 2: Handles formats like "GENE c.123A>G" or "NM_12345.6 c.123A>G"
    match_space = re.search(r'([A-Z0-9\-_()]+)\s+([cgnmp]\..*)', query, re.IGNORECASE)
    if match_space:
        identifier = match_space.group(1).strip()
        variant = match_space.group(2).strip()
        if identifier.startswith("NM_") or identifier.startswith("ENST"):
            transcript_id, gene_from_parens = _split_transcript_and_gene(identifier)
            return f"{identifier}:{variant}", gene_from_parens, transcript_id
        return f"{identifier}:{variant}", identifier, None

    return None, None, None

def classify_variant_clinsig(clinsig_field):
    if clinsig_field is None: return 'other'
    vals = [v.lower() for v in (clinsig_field if isinstance(clinsig_field, list) else [clinsig_field]) if isinstance(v, str)]
    # 'conflicting interpretations/classifications of pathogenicity' contains the word
    # 'pathogenic' but is not a pathogenic assertion
    vals = [v for v in vals if 'conflicting' not in v]
    if any('pathogenic' in v for v in vals): return 'pathogenic'
    if any('uncertain' in v for v in vals): return 'VUS'
    if any('benign' in v for v in vals): return 'benign'
    return 'other'

def _build_variant_link(variant: Dict[str, Any]) -> Optional[str]:
    """
    Builds a best-effort external link for a variant returned by Ensembl overlap endpoints.
    Prefers ClinVar links when available, otherwise falls back to dbSNP or Ensembl Variation explorer.
    """
    variant_id = str(variant.get('id') or variant.get('variation_name') or '').strip()
    if not variant_id:
        return None

    clinvar_id = variant.get('clinvar_variation_id') or variant.get('clinvar_id')
    if clinvar_id:
        cv_id = str(clinvar_id).strip()
        if cv_id:
            return f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{cv_id}/"

    source = str(variant.get('source') or '').lower()
    if source == 'clinvar' and variant_id.isdigit():
        return f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{variant_id}/"

    if variant_id.lower().startswith('rs'):
        return f"https://www.ncbi.nlm.nih.gov/snp/{variant_id}"

    return f"https://www.ensembl.org/Homo_sapiens/Variation/Explore?v={variant_id}"

def _ncbi_request(endpoint, params, method="GET", max_retries=3):
    p = dict(params)
    p["retmode"] = "json"
    p["tool"] = NCBI_TOOL
    if NCBI_EMAIL:
        p["email"] = NCBI_EMAIL
    if NCBI_API_KEY:
        p["api_key"] = NCBI_API_KEY
    url = f"{NCBI_EUTILS}/{endpoint}"
    backoff = 0.5
    _t0 = time.perf_counter()
    try:
        for _ in range(max_retries):
            try:
                with _NCBI_CONCURRENCY:
                    if method == "POST":
                        r = requests.post(url, data=p, timeout=120)
                    else:
                        r = requests.get(url, params=p, timeout=120)
                if r.status_code == 200:
                    return r.json()
                if r.status_code == 429:
                    time.sleep(backoff)
                    backoff *= 2
                    continue
                return None
            except requests.RequestException:
                time.sleep(backoff)
                backoff *= 2
        return None
    finally:
        _log_timing(f"NCBI ClinVar: {endpoint}", _t0)

def fetch_clinvar_region_cached(chrom, start, end):
    """
    Returns dict keyed by (chrom, pos) → list of {stars, ref, alt, variation_id}.
    Returns None on API failure so callers can skip the filter rather than silently
    treating every variant as unreviewed (which would bias toward eligibility).
    Returns {} when the API succeeds but finds no ClinVar variants in the region.
    Cached per region in the shared TTL cache, which doesn't memoize a failure's None -
    lru_cache would pin a transient NCBI error as "filter off" until restart.
    """
    return _ENSEMBL_LOOKUP_CACHE.get_or_set(('clinvar_region', chrom, start, end), lambda: _fetch_clinvar_region(chrom, start, end))

def _fetch_clinvar_region(chrom, start, end):
    # esearch only counts the region's records and leaves their IDs on NCBI's history server;
    # the summaries are then paged from there, so a large exon (BRCA1 exon 10 has >6,000 ClinVar
    # records) is not silently cut off after the first page.
    search_resp = _ncbi_request("esearch.fcgi", {
        "db": "clinvar",
        "term": f"{chrom}[chr] AND {start}:{end}[chrpos38]",
        "retmax": 0,
        "usehistory": "y",
    })
    if search_resp is None:
        return None
    result = search_resp.get("esearchresult", {})
    try:
        count = int(result.get("count", 0))
    except (TypeError, ValueError):
        return None
    if count == 0:
        return {}
    webenv, query_key = result.get("webenv"), result.get("querykey")
    if not webenv or not query_key:
        return None

    def fetch_page(retstart):
        return _ncbi_request("esummary.fcgi", {
            "db": "clinvar", "WebEnv": webenv, "query_key": query_key,
            "retstart": retstart, "retmax": CLINVAR_PAGE_SIZE,
        })
    region = {}
    for page in _NCBI_EXECUTOR.map(fetch_page, range(0, count, CLINVAR_PAGE_SIZE)):
        # A partial result would leave the missing records' variants unfiltered, so a failed page
        # fails the whole lookup (the caller then reports the star filter as not applied)
        if page is None:
            return None
        _add_clinvar_summaries(region, chrom, page)
    return region

def _add_clinvar_summaries(region: Dict[Tuple[str, int], List[Dict[str, Any]]], chrom: str, summary_resp: Dict[str, Any]) -> None:
    """Adds the records of one ClinVar esummary response to region, keyed by (chrom, position).
    A record is keyed under both its GRCh38 location start and its SPDI position (0-based + 1)
    when the two differ: for indels in repeats they can, and Ensembl's variant start matches
    one or the other."""
    for variation_id, rec in (summary_resp.get("result") or {}).items():
        if variation_id == "uids" or not isinstance(rec, dict):
            continue
        germline = rec.get("germline_classification") or {}
        review_status = (germline.get("review_status") or "").strip().lower()
        stars = REVIEW_STATUS_TO_STARS.get(review_status, 0)

        # Location and SPDI are reported per allele on variation_set; the first allele is used
        variation = (rec.get("variation_set") or [{}])[0] or {}
        positions = set()
        for loc in (variation.get("variation_loc") or []):
            if str(loc.get("assembly_name") or "").startswith("GRCh38"):
                try:
                    positions.add(int(loc.get("start")))
                except (TypeError, ValueError):
                    pass
                break
        ref_allele, alt_allele = None, None
        parts = str(variation.get("canonical_spdi") or "").split(":")
        if len(parts) == 4:
            try:
                positions.add(int(parts[1]) + 1)
            except ValueError:
                pass
            ref_allele, alt_allele = parts[2], parts[3]

        entry = {"stars": stars, "ref": ref_allele, "alt": alt_allele, "variation_id": variation_id}
        for pos in positions:
            region.setdefault((chrom, pos), []).append(entry)

def get_ensembl_variant_stars(variant, clinvar_region):
    """
    Returns the ClinVar star rating for an Ensembl variant dict.
    Tries allele-exact match (chr+pos+ref+alt) first, then falls back to the
    max stars at that position for multi-allelic sites.
    Returns -1 when the variant is not present in clinvar_region.
    """
    if not clinvar_region:
        return -1
    chrom = str(variant.get("seq_region_name") or "").strip()
    pos = variant.get("start")
    if not chrom or pos is None:
        return -1
    key = (chrom, int(pos))
    entries = clinvar_region.get(key)
    if not entries:
        return -1

    ref = str(variant.get("ref_allele") or variant.get("ref") or "").upper()
    raw_alt = variant.get("alt_alleles") or variant.get("alt") or ""
    alt = (raw_alt[0] if isinstance(raw_alt, list) else str(raw_alt)).upper()

    if ref and alt:
        exact = [e["stars"] for e in entries
                 if (e.get("ref") or "").upper() == ref and (e.get("alt") or "").upper() == alt]
        if exact:
            return max(exact)

    return max(e["stars"] for e in entries)

def count_variants_in_region(variants_in_region, clinvar_region):
    """
    Counts pathogenic/benign variants in variants_in_region, applying MIN_REVIEW_STARS
    when clinvar_region is not None. When clinvar_region is None (NCBI API failure),
    the filter is skipped so the caller's result degrades gracefully rather than
    silently clearing all counts.

    Returns (counts, pathogenic_variant_links, splice_variant_links, star_filter_active).
    star_filter_active is False only when the NCBI lookup failed.
    """
    star_filter_active = clinvar_region is not None

    counts = {'missense': 0, 'inframe_del': 0, 'splice': 0, 'nonsense': 0, 'frameshift': 0, 'benign_splice': 0}
    pathogenic_variant_links = {'missense': [], 'inframe_del': [], 'splice': [], 'nonsense': [], 'frameshift': []}
    seen_pathogenic = {'missense': set(), 'inframe_del': set(), 'splice': set(), 'nonsense': set(), 'frameshift': set()}
    splice_variant_links = {"pathogenic": [], "benign": []}
    seen_splice = {"pathogenic": set(), "benign": set()}

    for v in variants_in_region:
        clclass = classify_variant_clinsig(v.get('clinical_significance'))
        if clclass not in ("pathogenic", "benign"):
            continue
        if star_filter_active:
            stars = get_ensembl_variant_stars(v, clinvar_region)
            if stars != -1 and stars < MIN_REVIEW_STARS:
                continue

        conseq = (v.get("consequence_type") or "").lower()
        vid = str(v.get('id') or v.get('variation_name') or '').strip()
        vsrc = (v.get('source') or '').strip() or "Ensembl"
        vlink = _build_variant_link(v) if vid else None

        if clclass == "pathogenic":
            vtype = None
            if "missense" in conseq: counts['missense'] += 1; vtype = 'missense'
            elif "inframe_deletion" in conseq: counts['inframe_del'] += 1; vtype = 'inframe_del'
            elif "splice_donor" in conseq or "splice_acceptor" in conseq: counts['splice'] += 1; vtype = 'splice'
            elif "stop_gained" in conseq: counts['nonsense'] += 1; vtype = 'nonsense'
            elif "frameshift" in conseq: counts['frameshift'] += 1; vtype = 'frameshift'
            if vtype and vid and vid not in seen_pathogenic[vtype]:
                pathogenic_variant_links[vtype].append({"id": vid, "source": vsrc, "url": vlink})
                seen_pathogenic[vtype].add(vid)
            if vtype == 'splice' and vid and vid not in seen_splice["pathogenic"]:
                splice_variant_links["pathogenic"].append({"id": vid, "source": vsrc, "url": vlink})
                seen_splice["pathogenic"].add(vid)
        elif clclass == 'benign' and ("splice_donor" in conseq or "splice_acceptor" in conseq):
            counts['benign_splice'] += 1
            if vid and vid not in seen_splice["benign"]:
                splice_variant_links["benign"].append({"id": vid, "source": vsrc, "url": vlink})
                seen_splice["benign"].add(vid)

    return counts, pathogenic_variant_links, splice_variant_links, star_filter_active

# --- TReXome DB (tandem repeat domain) access ---

def _trexome_request(endpoint: str, params: Dict[str, Any], max_retries: int = 3) -> Optional[Dict[str, Any]]:
    """GETs a TReXome DB API endpoint, retrying on transient errors. Returns the parsed JSON body, or None on failure."""
    url = f"{TREXOME_API_BASE}{endpoint}"
    backoff = 1.0
    _t0 = time.perf_counter()
    try:
        for attempt in range(max_retries):
            try:
                resp = requests.get(url, params=params, timeout=120)
                if resp.status_code == 200:
                    return resp.json()
                if resp.status_code not in (429, 500, 502, 503, 504) or attempt + 1 == max_retries:
                    return None
            except requests.RequestException:
                if attempt + 1 == max_retries:
                    return None
            time.sleep(backoff)
            backoff *= 2
        return None
    finally:
        _log_timing(f"TReXome: {endpoint}", _t0)

def _trexome_detail_url(uniprot_id: Optional[str]) -> Optional[str]:
    """Builds a link to a repeat's detail page on the TReXome DB frontend, keyed by UniProt accession."""
    if not uniprot_id:
        return None
    return f"{TREXOME_FRONTEND_BASE}/detail.html?id={uniprot_id}"

def _trexome_repeats_for_gene(gene_symbol: str) -> Tuple[Dict[str, Any], ...]:
    """Fetches and caches every repeat-exon match row for a gene from TReXome DB's /api/repeats (paginated).
    A failed page yields no rows rather than a partial list, cached for _TREXOME_CACHE's short TTL."""
    if not gene_symbol:
        return tuple()
    def fetch():
        all_rows: List[Dict[str, Any]] = []
        page = 1
        while True:
            data = _trexome_request("/api/repeats", {"gene": gene_symbol, "page": page, "page_size": 200})
            if not data:
                return tuple()
            all_rows.extend(data.get("results") or [])
            if page >= data.get("total_pages", 1):
                return tuple(all_rows)
            page += 1
    return _TREXOME_CACHE.get_or_set(('trexome_repeats', gene_symbol), fetch)

def _trexome_overlap_percentage(row: Dict[str, Any]) -> Optional[float]:
    """Extracts the exon/repeat overlap percentage for a repeat-exon match row. TReXome DB only
    exposes this figure nested inside the row's original source object (source.ensembl_exon_info.
    transcripts[].containing_exons[]), keyed by exon_id within the matching transcript."""
    try:
        transcripts = (((row.get('source') or {}).get('ensembl_exon_info') or {}).get('transcripts')) or []
        row_transcript_id = row.get('transcript_id')
        row_exon_id = row.get('exon_id')
        for t in transcripts:
            if row_transcript_id and t.get('transcript_id') != row_transcript_id:
                continue
            for ce in (t.get('containing_exons') or []):
                if ce.get('exon_id') == row_exon_id:
                    return ce.get('overlap_percentage')
    except Exception:
        pass
    return None

def get_trexome_repeats_for_exon(gene_symbol: Optional[str], chrom: Optional[str], exon_start: Optional[int], exon_end: Optional[int]) -> List[Dict[str, Any]]:
    """Returns TReXome DB tandem repeat domain records overlapping a specific exon, identified by
    genomic coordinates rather than Ensembl exon ID. Matching on the exon's own (chrom, start, end)
    against TReXome's flattened exon_start/exon_end is more robust than comparing exon_id strings,
    since TReXome may index a different (but coordinate-identical) transcript than the one chosen
    here, or a different Ensembl release where the same exon was assigned a different stable ID."""
    if not gene_symbol or not chrom or exon_start is None or exon_end is None:
        return []
    norm_chrom = str(chrom).strip().lower()
    if not norm_chrom.startswith('chr'):
        norm_chrom = f"chr{norm_chrom}"
    rows = _trexome_repeats_for_gene(gene_symbol)
    seen = set()
    repeats = []
    for r in rows:
        row_chrom = str(r.get('chrom') or '').strip().lower()
        row_start, row_end = r.get('exon_start'), r.get('exon_end')
        if row_chrom != norm_chrom or row_start is None or row_end is None:
            continue
        if max(exon_start, row_start) > min(exon_end, row_end):
            continue
        key = (r.get('uniprot_id'), r.get('protein_start'), r.get('protein_end'))
        if key in seen:
            continue
        seen.add(key)
        repeats.append({
            "repeat_type": r.get('repeat_type'),
            "uniprot_id": r.get('uniprot_id'),
            "status": r.get('status'),
            "protein_start": r.get('protein_start'),
            "protein_end": r.get('protein_end'),
            "overlap_percentage": _trexome_overlap_percentage(r),
            "url": _trexome_detail_url(r.get('uniprot_id')),
        })
    return repeats

def get_trexome_gene_repeat_summary(gene_symbol: str) -> Dict[str, Any]:
    """Fetches a gene-level tandem repeat domain summary (total repeat count across the gene's
    annotated proteins, repeat types present, and a representative TReXome DB link) via /api/proteins."""
    empty = {"repeat_count": 0, "repeat_types": [], "url": None}
    if not gene_symbol:
        return empty
    def fetch():
        data = _trexome_request("/api/proteins", {"gene": gene_symbol, "page_size": 200})
        results = (data or {}).get("results") or []
        if not results:
            return empty
        return {
            "repeat_count": sum(r.get('repeat_count') or 0 for r in results),
            "repeat_types": sorted({r.get('repeat_type') for r in results if r.get('repeat_type')}),
            "url": _trexome_detail_url(results[0].get('uniprot_id')),
        }
    return dict(_TREXOME_CACHE.get_or_set(('trexome_summary', gene_symbol), fetch))

def _base_transcript_id(transcript_id: Optional[str]) -> str:
    """Strips a trailing '.<version>' suffix (e.g. 'NM_000546.5' -> 'NM_000546') for comparison."""
    if not transcript_id:
        return ""
    return re.sub(r'\.\d+$', '', transcript_id.strip()).upper()

def _find_consequence_by_transcript(consequences: List[Dict[str, Any]], transcript_id: Optional[str]) -> Optional[Dict[str, Any]]:
    """
    Finds the consequence matching a user-supplied transcript ID, ignoring version suffix.
    VEP's `transcript_id` is always an Ensembl ID (ENST...); a RefSeq accession (NM_...) typed
    by the user only ever shows up cross-referenced via `mane_select` (when it's the MANE Select
    transcript) or `xref_refseq` (a best-effort RefSeq cross-reference Ensembl sometimes provides),
    so both are checked as well.
    """
    if not transcript_id:
        return None
    target_base = _base_transcript_id(transcript_id)
    for c in consequences:
        if _base_transcript_id(c.get('transcript_id')) == target_base:
            return c
    for c in consequences:
        if _base_transcript_id(c.get('mane_select')) == target_base:
            return c
    for c in consequences:
        xref = c.get('xref_refseq')
        if not xref:
            continue
        xref_ids = xref if isinstance(xref, list) else str(xref).split(',')
        if any(_base_transcript_id(x) == target_base for x in xref_ids):
            return c
    return None

def choose_best_consequence(consequences: List[Dict[str, Any]], gene_symbol_from_query: Optional[str] = None) -> Optional[Dict[str, Any]]:
    """
    Selects the most relevant transcript consequence from a VEP response (MANE Select ->
    Ensembl canonical -> first protein-coding transcript whose CDS the variant falls in),
    restricted to the queried gene if given.
    """
    if not consequences: return None

    if gene_symbol_from_query:
        targeted_consequences = [
            c for c in consequences if (c.get('gene_symbol') or '').upper() == gene_symbol_from_query.upper()
        ]
        if targeted_consequences:
            consequences = targeted_consequences

    mane_select = [c for c in consequences if c.get('mane_select')]
    if mane_select: return mane_select[0]

    canonical_hits = [c for c in consequences if c.get('canonical') == 1]
    if canonical_hits: return canonical_hits[0]

    # VEP's cds_start/cds_end are the variant's positions within each CDS, not the transcript's
    # length, so they can't rank transcripts; they only show that the variant hits the CDS
    coding_consequences = [c for c in consequences if c.get('biotype') == 'protein_coding' and c.get('cds_end')]
    return coding_consequences[0] if coding_consequences else consequences[0]

def extract_exons_from_transcript(transcript: Dict[str, Any]):
    exons_raw = sorted(transcript.get('Exon', []), key=lambda e: e['start'])
    if transcript.get('strand') == -1: exons_raw.reverse()
    
    translation_data = transcript.get('Translation', {})
    cds_start, cds_end = translation_data.get('start'), translation_data.get('end')
    seq_region = transcript.get('seq_region_name')
    
    normalized, coding_exon_count = [], 0
    for i, e in enumerate(exons_raw, 1):
        start, end = e['start'], e['end']
        cds_len_of_exon, is_coding = 0, False
        if cds_start and cds_end:
            overlap_start, overlap_end = max(start, cds_start), min(end, cds_end)
            if overlap_end >= overlap_start:
                cds_len_of_exon = overlap_end - overlap_start + 1
                is_coding = True
                coding_exon_count += 1
        
        normalized.append({
            'total_exon_number': i, 'coding_exon_number': coding_exon_count if is_coding else None,
            'exon_id': e.get('id'), 'start': start, 'end': end,
            'seq_region_name': seq_region, 'cds_length': cds_len_of_exon
        })
    return normalized

# --- ASO Strategy Assessment Logic ---

def _n1c_rows_for_gene(gene_symbol: Optional[str]) -> pd.DataFrame:
    """N1C registry rows for a gene (case-insensitive); empty if the registry or gene is unavailable."""
    df = n1c_variants_df  # one reference: a background refresh may swap the global mid-call
    if df is None or df.empty or not gene_symbol or 'Gene' not in df.columns:
        return pd.DataFrame()
    return df[df['Gene'].str.upper() == gene_symbol.upper()]

def _n1c_entry_link(row: pd.Series) -> Optional[str]:
    """Registry entry page for an N1C registry row, or None if the row has no ID."""
    nid = row.get('ID')
    if nid is None or pd.isna(nid) or not str(nid).strip() or str(nid).strip().lower() == 'nan':
        return None
    try:
        nid = int(float(nid))
    except (ValueError, TypeError):
        nid = str(nid).strip()
    return f"https://generegistry.n1collaborative.org/entry.html?id={quote(str(nid))}"

# --- N1C field interpretation ---
# Both N1C tables are read by their actual column names. Free-text columns (comments,
# indications, publications) are never searched for approach keywords or exon numbers.

def _row_text(row: pd.Series, column: str) -> str:
    """A table cell as stripped text; '' when the column is missing or the value is empty."""
    value = row.get(column)
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ''
    return str(value).strip()

def _row_is_aso_project(row: pd.Series) -> bool:
    """True if the registry row's Therapeutic Modality is an antisense oligonucleotide."""
    return bool(re.search(r'antisense|\baso\b', _row_text(row, 'Therapeutic Modality'), re.IGNORECASE))

def _row_is_exon_skipping(row: pd.Series) -> bool:
    """True if the registry row's Approach is exon skipping."""
    return bool(re.search(r'exon\s*skip', _row_text(row, 'Approach'), re.IGNORECASE))

def _row_is_knockdown(row: pd.Series) -> bool:
    """True if the registry row's Approach is knockdown (including allele-selective knockdown)."""
    return bool(re.search(r'knock\s*down|silenc', _row_text(row, 'Approach'), re.IGNORECASE))

# Registry Status values meaning the ASO exists, i.e. is in clinical use or has been developed
N1C_ESTABLISHED_STATUSES = ('clinical application', 'developed')

def _row_project_state(row: pd.Series) -> str:
    """How far a registry project is: 'established' (an ASO with Status 'Clinical Application' or
    'Developed'), 'under_development' (an ASO with Status 'Under Development'), or 'other' (not an
    ASO, e.g. gene editing or gene replacement, or any other status such as 'Currently Not
    Pursued for Development')."""
    if not _row_is_aso_project(row):
        return 'other'
    status = _row_text(row, 'Status').lower()
    if status in N1C_ESTABLISHED_STATUSES:
        return 'established'
    if status == 'under development':
        return 'under_development'
    return 'other'

_PROJECT_STATE_PRIORITY = {'established': 0, 'under_development': 1, 'other': 2}

def _cdna_changes(text: Optional[str]) -> List[str]:
    """The c. changes named in a text, normalised for exact comparison: lower-cased, cut before
    any parenthetical note, and without the optional deleted/duplicated bases that VEP omits
    ('c.696dupG' -> 'c.696dup'). A registry cell can name more than one change."""
    changes = []
    for token in re.findall(r'c\.[^\s;,()]+', text or '', re.IGNORECASE):
        token = re.sub(r'(del|dup)[acgt]+$', r'\1', token.lower())
        if len(token) > 2:
            changes.append(token)
    return changes

def _same_cdna_change(query: str, registry_cell: str) -> bool:
    """True if the registry cell names exactly the c. change of the query (after the
    normalisation in _cdna_changes). Substring matching is not enough: 'c.1087del' must not
    match 'c.1087delinsTT'."""
    query_changes = _cdna_changes(query)
    return bool(query_changes) and query_changes[0] in _cdna_changes(registry_cell)

def check_n1c_registry(gene_symbol: str, formatted_hgvs: str) -> Optional[Dict[str, Any]]:
    """
    Searches the pre-loaded N1C registry for a project on exactly this variant. An established
    ASO project makes the variant 'Eligible'; an ASO project under development is shown as
    'Under Development'. Any other entry (a non-ASO modality such as gene editing, or a project no
    longer pursued) is only shown for information: it is flagged 'informational' and the ASO
    assessment runs as usual. If several rows match, the most advanced ASO project is reported.
    """
    matches = [
        row for _, row in _n1c_rows_for_gene(gene_symbol).iterrows()
        if _same_cdna_change(formatted_hgvs, _row_text(row, 'Coding DNA change (c.)'))
    ]
    if not matches:
        return None
    row = min(matches, key=lambda r: _PROJECT_STATE_PRIORITY[_row_project_state(r)])
    state = _row_project_state(row)

    # Registry fields are third-party content rendered as HTML by the frontend - escape them
    registry_c_dot = _row_text(row, 'Coding DNA change (c.)')
    status = _row_text(row, 'Status') or 'N/A'
    modality = _row_text(row, 'Therapeutic Modality') or 'N/A'
    approach = _row_text(row, 'Approach')
    link = _n1c_entry_link(row)
    link_html = f"<a href='{escape(link)}' target='_blank' rel='noopener noreferrer'>Click here to view the N=1 Collaborative registry page.</a>" if link else ""
    reason = (
        f"<p>A direct match for a variant in <strong>{escape(gene_symbol)}</strong> was found in the N=1 Collaborative Projects Registry.</p>"
        f"<ul>"
        f"<li><strong>Matched Variant:</strong> {escape(registry_c_dot)}</li>"
        f"<li><strong>Status:</strong> {escape(status)}</li>"
        f"<li><strong>Therapeutic Modality:</strong> {escape(modality)}</li>"
        + (f"<li><strong>Approach:</strong> {escape(approach)}</li>" if approach else "")
        + f"</ul>"
    )
    if state == 'other':
        reason += "<p>This entry is not an existing or developing antisense oligonucleotide, so it does not establish ASO eligibility; the ASO assessment below applies.</p>"
    reason += link_html

    classification = {'established': "Eligible", 'under_development': "Under Development", 'other': "Informational"}[state]
    result = {"classification": classification, "reason": reason, "link": link}
    if state == 'under_development':
        result["under_development"] = True
    elif state == 'other':
        result["informational"] = True
    return result

# Curated eligibility values in the N1C assessed-variants table, compared after lower-casing and
# removing everything but letters. Spelling variants found in the table are listed explicitly
# instead of matching loosely, so a value such as 'Noteligible' can never be read as 'eligible';
# a value not listed here is shown verbatim and classified 'Unable to Assess'.
N1C_CURATED_ELIGIBILITY = {
    'eligible': 'Eligible',
    'likelyeligible': 'Likely Eligible',
    'unlikelyeligible': 'Unlikely Eligible',
    'noteligible': 'Not Eligible',
    'noteligiblle': 'Not Eligible',
    'noteligibe': 'Not Eligible',
    'unabletoassess': 'Unable to Assess',
    'underdevelopment': 'Under Development',
}

def _curated_eligibility(value: str) -> Optional[str]:
    """Maps a curated Eligibility value to an AVEC classification, or None if it is not recognised."""
    return N1C_CURATED_ELIGIBILITY.get(re.sub(r'[^a-z]', '', value.lower()))

def check_n1c_assessed_variants(gene_symbol: str, formatted_hgvs: str) -> Optional[Dict[str, Any]]:
    """Returns the N1C curated assessment of exactly this variant, if the assessed-variants table
    has one with an eligibility call."""
    assessed_df = n1c_assessed_df  # one reference: a background refresh may swap the global mid-call
    if assessed_df is None or assessed_df.empty or not gene_symbol or 'Gene' not in assessed_df.columns:
        return None
    if not _cdna_changes(formatted_hgvs):
        return None

    for _, row in assessed_df[assessed_df['Gene'].str.upper() == gene_symbol.upper()].iterrows():
        c_dot = _row_text(row, 'Coding DNA change (c.)')
        eligibility = _row_text(row, 'Eligibility')
        # Rows without an eligibility call carry no curated assessment
        if not eligibility or not _same_cdna_change(formatted_hgvs, c_dot):
            continue

        classification = _curated_eligibility(eligibility)
        link = None
        try:
            link = f"https://generegistry.n1collaborative.org/variant_entry.html?id={int(float(row.get('ID')))}"
        except (ValueError, TypeError):
            pass

        # Registry fields are third-party content rendered as HTML by the frontend - escape them
        reason = (
            f"<p>A curated assessment exists in the N=1 Collaborative Assessed Variants database for <strong>{escape(gene_symbol)}</strong>.</p>"
            f"<p><strong>Matched Variant:</strong> {escape(c_dot)}</p>"
        )
        if classification is None:
            reason += f"<p><strong>Curated eligibility:</strong> {escape(eligibility)} (not a value AVEC recognises; see the curated entry)</p>"
        approach = _row_text(row, 'Approach')
        if approach:
            reason += f"<p><strong>Approach:</strong> {escape(approach)}</p>"
        assessment_text = _row_text(row, 'Assessment text')
        if assessment_text:
            reason += f"<p><strong>Assessment:</strong> {escape(assessment_text)}</p>"
        if link:
            reason += f"<p><a href='{escape(link)}' target='_blank' rel='noopener noreferrer'>View N1C curated entry</a></p>"
        return {
            "classification": classification or "Unable to Assess",
            "reason": reason,
            "link": link
        }
    return None

def n1c_exon_skipping_variant_exon_map(gene_symbol: str, all_exons: List[Dict[str, Any]], client: EnsemblClient) -> Tuple[Dict[int, List[str]], Dict[int, bool]]:
    """
    Maps the gene's N1C exon-skipping ASO projects (RefSeq Transcript + c. change) to exon numbers
    of this transcript, using the same HGVS parsing as user input. Returns
    ({exon_number: [links]}, {exon_number: is_under_development}); an exon is a key of the first
    map even when its projects have no registry link. is_under_development[exon] is True only if
    every matching project for that exon is still under development. Projects that are neither
    established nor under development (see _row_project_state) are ignored.
    """
    exon_map: Dict[int, List[str]] = {}
    exon_status_map: Dict[int, set] = {}
    if not all_exons:
        return exon_map, {}

    # --- Pass 1: resolve each row's HGVS query (no network yet) ---
    row_queries: List[Tuple[pd.Series, bool, str]] = []
    for _, row in _n1c_rows_for_gene(gene_symbol).iterrows():
        state = _row_project_state(row)
        if not _row_is_exon_skipping(row) or state == 'other':
            continue
        refseq_val = _row_text(row, 'RefSeq Transcript')
        cdot_val = _row_text(row, 'Coding DNA change (c.)')
        if not refseq_val or not cdot_val:
            continue
        try:
            hgvs_query, _, _ = parse_hgvs_query(f"{refseq_val}:{cdot_val}")
        except Exception:
            continue
        if not hgvs_query:
            continue
        row_queries.append((row, state == 'under_development', hgvs_query))

    # --- Pass 2: fetch VEP for every distinct query concurrently. Rows are per-gene
    # registry entries so duplicate HGVS strings across rows are common; dedupe so each
    # one is only fetched once, regardless of how many rows reference it. ---
    unique_queries = {hgvs_query for _, _, hgvs_query in row_queries}
    query_futures = {q: _IO_EXECUTOR.submit(client.vep_hgvs, q) for q in unique_queries}

    # --- Pass 3: map each row's variant to the exon it falls in, in original row order ---
    for row, is_under_development, hgvs_query in row_queries:
        try:
            vep_data = query_futures[hgvs_query].result()
            if not vep_data or not isinstance(vep_data, list):
                continue
            vep_entry = vep_data[0]
            v_chrom, v_start, v_end = vep_entry.get('seq_region_name'), vep_entry.get('start'), vep_entry.get('end')
            if not all([v_chrom, v_start, v_end]):
                continue
            target_exon = next((ex for ex in all_exons if ex['seq_region_name'] == v_chrom and max(v_start, ex['start']) <= min(v_end, ex['end'])), None)
            if not target_exon:
                continue
            exon_num = target_exon.get('total_exon_number')
            if exon_num is None:
                continue
            link = _n1c_entry_link(row)
            exon_links = exon_map.setdefault(exon_num, [])  # key exists even without a link
            if link and link not in exon_links:
                exon_links.append(link)
            exon_status_map.setdefault(exon_num, set()).add(is_under_development)
        except Exception:
            continue

    under_development_map: Dict[int, bool] = {
        n: all(statuses) for n, statuses in exon_status_map.items()
    }
    return exon_map, under_development_map

_MOA_LABELS = {"GoF": "gain-of-function", "DN": "dominant-negative", "LoF": "loss-of-function"}

def n1c_gene_knockdown_entry(gene_symbol: str, moa: str) -> Optional[Dict[str, Any]]:
    """Returns the gene's N1C knockdown ASO project (Approach: Knockdown), preferring an
    established project over one under development; other entries are ignored. Only called when
    the selected mechanism (moa) and the gene's inheritance make knockdown applicable."""
    candidates = [
        row for _, row in _n1c_rows_for_gene(gene_symbol).iterrows()
        if _row_is_knockdown(row) and _row_project_state(row) != 'other'
    ]
    if not candidates:
        return None
    row = min(candidates, key=lambda r: _PROJECT_STATE_PRIORITY[_row_project_state(r)])
    link = _n1c_entry_link(row)
    modality = _row_text(row, 'Therapeutic Modality') or _row_text(row, 'Approach')
    status = _row_text(row, 'Status')
    moa_label = _MOA_LABELS.get(moa, moa)

    # A registry match that is still "Under Development" should not be treated as a final
    # "Eligible" result; surface it with the "Under Development" styling and let the regular
    # knockdown assessment run in addition.
    is_under_development = _row_project_state(row) == 'under_development'
    if is_under_development:
        classification = "Under Development"
        reason = (
            f"N1C registry lists a knockdown project for {gene_symbol} that is still under development. "
            f"With the selected {moa_label} mechanism and the gene's mode of inheritance, knockdown is an applicable strategy; "
            f"the regular knockdown assessment is shown below in addition to this registry match."
        )
    else:
        classification = "Eligible"
        reason = (
            f"N1C registry lists knockdown project(s) for {gene_symbol}. "
            f"With the selected {moa_label} mechanism and the gene's mode of inheritance, knockdown is an applicable strategy, so the variant is considered eligible for knockdown."
        )
    details = {}
    if modality:
        details["Registry Modality"] = modality
    if status:
        details["Registry Status"] = status
    if link:
        details["N1C Registry Link"] = link
    result = {
        "classification": classification,
        "reason": reason,
        "link": link,
        "details": details
    }
    if is_under_development:
        result["under_development"] = True
    return result

# ClinGen mode-of-inheritance terms that don't say whether the gene acts dominantly or recessively
# ('Mode of inheritance' is ClinGen's placeholder for an undetermined one)
UNSPECIFIC_MOI_TERMS = {'x-linked inheritance', 'semidominant inheritance', 'mode of inheritance'}

def _add_summary_note(result: Dict[str, Any], note: str) -> None:
    """Appends a note to a result's summary note; several notes are separated by <br>."""
    summary = result.get("summary")
    if isinstance(summary, dict):
        summary["note"] = f"{summary['note']}<br>{note}" if summary.get("note") else note

def get_gene_characteristics(gene_symbol: str) -> Dict[str, Any]:
    """
    Retrieves MOI, Haploinsufficiency and Triplosensitivity (ClinGen) and rCNV scores
    (pHaplo/pTriplo, Collins et al. 2022). The TReXome gene-level repeat summary starts out
    empty here; _assess_for_transcript fetches it concurrently and fills it in.
    The mechanism of action is always selected by the user, never looked up.
    """
    characteristics = {
        "moi": [],
        "haploinsufficiency": {"text": "Unknown", "url": None},
        "triplosensitivity": {"text": "Unknown", "url": None},
        "gene_url": None,
        "rcnv": {"pHaplo": "N/A", "pTriplo": "N/A", "url": "https://doi.org/10.1016/j.cell.2022.06.036"}, # Init rCNV
        "trexome_repeats": {"repeat_count": 0, "repeat_types": [], "url": None}
    }

    # Use the global rcnv_df loaded in load_databases()
    if rcnv_df is not None and gene_symbol:
        # Exact match (case-insensitive handled by load step)
        match = rcnv_df[rcnv_df['gene'] == gene_symbol.upper()]
        if not match.empty:
            ph = match.iloc[0].get('pHaplo')
            pt = match.iloc[0].get('pTriplo')
            # Format to 3 decimal places if numeric
            try:
                if pd.notna(ph): characteristics['rcnv']['pHaplo'] = f"{float(ph):.3f}"
                if pd.notna(pt): characteristics['rcnv']['pTriplo'] = f"{float(pt):.3f}"
            except (ValueError, TypeError):
                pass

    # --- ClinGen Lookup ---
    if clingen_df is not None and gene_symbol in clingen_df.index:
        # The report has one row per gene-disease curation; gene-level fields (URL, dosage) are
        # identical across a gene's rows, but mode of inheritance differs per disease
        gene_rows = clingen_df.loc[[gene_symbol]]
        gene_data_row = gene_rows.iloc[0]

        gene_url_val = gene_data_row.get('gene_url')
        if pd.notna(gene_url_val):
            characteristics['gene_url'] = gene_url_val

        # 1. Mode of Inheritance (MOI), collected across all of the gene's curations; a single
        # curation can list several, separated by ' | '
        all_mois = {moi.strip() for moi_val in gene_rows['mode_of_inheritance'].dropna() for moi in re.split(r'[,|]', str(moi_val)) if moi.strip()}
        characteristics["moi"] = sorted(all_mois)

        # 2. Haploinsufficiency with Link
        hap_assertion = gene_data_row.get('dosage_haploinsufficiency_assertion')
        hap_url = gene_data_row.get('dosage_report')
        if pd.notna(hap_assertion):
            hap_score_str = str(hap_assertion).strip()
            haplo_text = "Unknown"
            if hap_score_str.startswith('3 -'): haplo_text = "Sufficient evidence"
            elif hap_score_str.startswith('2 -'): haplo_text = "Emerging evidence"
            elif hap_score_str.startswith('1 -'): haplo_text = "Little evidence"
            elif hap_score_str.startswith('0 -'): haplo_text = "No evidence"
            elif hap_score_str.startswith('30 -'): haplo_text = "Gene associated with autosomal recessive phenotype"
            elif hap_score_str.startswith('40 -'): haplo_text = "Dosage sensitivity unlikely"

            characteristics['haploinsufficiency'] = {
                "text": haplo_text,
                "url": hap_url if pd.notna(hap_url) else None
            }

        # 3. Triplosensitivity with Link
        triplo_assertion = gene_data_row.get('dosage_triplosensitivity_assertion')
        triplo_url = gene_data_row.get('dosage_report')
        if pd.notna(triplo_assertion):
            triplo_score_str = str(triplo_assertion).strip()
            triplo_text = "Unknown"
            if triplo_score_str.startswith('3 -'): triplo_text = "Sufficient evidence"
            elif triplo_score_str.startswith('2 -'): triplo_text = "Emerging evidence"
            elif triplo_score_str.startswith('1 -'): triplo_text = "Little evidence"
            elif triplo_score_str.startswith('30 -'): triplo_text = "Gene associated with autosomal recessive phenotype"
            elif triplo_score_str.startswith('0 -'): triplo_text = "No evidence"
            elif triplo_score_str.startswith('40 -'): triplo_text = "Sensitivity unlikely"
            
            characteristics['triplosensitivity'] = {
                "text": triplo_text,
                "url": triplo_url if pd.notna(triplo_url) else None
            }

    return characteristics

def assess_knockdown(gene_characteristics: Dict[str, Any], overrides: Optional[Dict[str, bool]] = None) -> Dict[str, Any]:
    """Assesses eligibility for a knockdown strategy."""
    overrides = overrides or {}
    haplo_obj = gene_characteristics.get("haploinsufficiency", {"text": "Unknown"})
    haplo_status_text = haplo_obj.get("text", "Unknown")

    # --- Logic to track original vs overridden state ---
    check_name = "Gene is not haploinsufficient"
    original_passed = haplo_status_text in ["No evidence", "Little evidence", "Dosage sensitivity unlikely"]

    has_override = check_name in overrides
    final_passed = overrides.get(check_name, original_passed)

    aso_exists_final = overrides.get('ASO already exists', False)

    checks = {
        check_name: {
            "passed": final_passed,
            "original_passed": original_passed,
            "overridden": has_override
        },
        "ASO already exists": {
            "passed": aso_exists_final,
            "original_passed": False,
            "overridden": 'ASO already exists' in overrides
        }
    }
    # Determine if haploinsufficiency status is unknown/missing for warning purposes
    haplo_unknown = (haplo_status_text is None) or (str(haplo_status_text).strip().lower() in ["", "unknown", "n/a"])

    if aso_exists_final:
        return {
            "classification": "Eligible",
            "reason": "An ASO for this gene already exists, confirming eligibility for knockdown.",
            "checks": checks
        }

    # Unknown status is checked first: it also fails the check, so after the 'not final_passed'
    # branch it was unreachable and unknown genes were reported as haploinsufficient
    if haplo_unknown and not has_override: # Only show 'Unable to Assess' if status is unknown AND not overridden
        reason = "Gene haploinsufficiency status is not available in datasources. Therefore knockdown could lead to unintended consequences and AVEC cannot determine amenability to knockdown."
        return {
            "classification": "Unable to Assess",
            "reason": reason,
            "checks": checks
        }

    # Classification logic now depends on the final 'is_not_haploinsufficient' value
    elif not final_passed:
        return {
            "classification": "Unlikely Eligible",
            "reason": "Gene is associated with haploinsufficiency or Autosomal Dominant LoF, therefore knockdown could lead to unintended consequences. If the haploinsufficiency phenotype is mild compared to the current phenotype, knockdown may still be considered.",
            "checks": checks
        }
    else:
        # Base reason
        reason = "Gene is assessed as insensitive to haploinsufficiency by ClinGen or by manual assessment."

        return {
            "classification": "Likely Eligible",
            "reason": reason,
            "checks": checks
        }

def assess_wt_upregulation(client: EnsemblClient, gene_id: str, gene_symbol: str, overrides: Optional[Dict[str, bool]] = None) -> Dict[str, Any]:
    """
    Assesses for WT upregulation by checking for overlapping NATs and by
    searching for a conventional antisense gene name ([GENE]-AS1).
    Trusts the '-AS1' naming convention without a strict biotype check.
    """
    overrides = overrides or {}
    if not gene_id or not gene_symbol:
        return {"classification": "Unable to Assess", "reason": "Missing Gene ID or Symbol."}

    aso_exists_final = overrides.get('ASO already exists', False)
    if aso_exists_final:
        return {
            "classification": "Eligible",
            "reason": "An ASO for this gene already exists, confirming eligibility for WT upregulation.",
            "checks": {
                "ASO already exists": {"passed": True, "original_passed": False, "overridden": True}
            }
        }

    found_antisense_genes = {}

    # Collect curated gene-level features from the N1C supplementary table.
    # Only columns holding an actual (non-N/A) value for this gene are surfaced;
    # missing/N/A entries are omitted entirely rather than shown as "Unknown".
    supp_details: Dict[str, str] = {}
    try:
        if n1c_supp_df is not None and 'Gene' in n1c_supp_df.columns:
            matches = n1c_supp_df[n1c_supp_df['Gene'] == str(gene_symbol).strip()]
            if not matches.empty:
                row = matches.iloc[0]
                for col in N1C_SUPP_DISPLAY_COLUMNS:
                    if col not in n1c_supp_df.columns:
                        continue
                    val = row.get(col)
                    if pd.isna(val):
                        continue
                    val_str = str(val).strip()
                    val_upper = val_str.upper()
                    if val_upper in ("Y", "YES"):
                        val_str = "Yes"
                    elif val_upper in ("N", "NO"):
                        val_str = "No"
                    supp_details[col] = val_str
    except Exception:
        pass  # Fail-quiet: do not block assessment if the supplemental table can't be parsed

    try:
        # --- Method 1: Search by genomic coordinate overlap ---
        # A NAT is a non-coding gene overlapping this one on the opposite strand. Ensembl folded
        # the old 'antisense' biotype into 'lncRNA' (release 100), so match either, by strand.
        overlapping_genes = client.get_overlapping_genes(gene_id)
        own_strand = next((g.get('strand') for g in overlapping_genes if g.get('id') == gene_id), None)
        for gene in overlapping_genes:
            if (gene.get('id') and gene.get('id') != gene_id
                    and gene.get('biotype') in ('antisense', 'lncRNA')
                    and own_strand is not None and gene.get('strand') == -own_strand):
                found_antisense_genes[gene['id']] = gene

        # --- Method 2: Search by conventional name ([GENE_SYMBOL]-AS1) ---
        as_gene = client.lookup_symbol(f"{gene_symbol}-AS1")
        if as_gene and as_gene.get('id'):
            found_antisense_genes[as_gene['id']] = as_gene

        # --- Evaluate curated evidence (any non-N/A supplementary column) and Ensembl NAT search ---
        curated_features = list(supp_details.keys())

        # Overlap results carry 'external_name'; lookup_symbol results carry 'display_name'
        nat_list = list(found_antisense_genes.values())
        nat_label = lambda nat: nat.get('external_name') or nat.get('display_name') or nat['id']
        nat_names = ", ".join(nat_label(nat) for nat in nat_list)
        nat_ids: List[str] = [nat['id'] for nat in nat_list]

        # Start details with curated statuses and cite the original supplementary table
        details = dict(supp_details)
        details["Datasource"] = "N1C VARIANT 2.0 Supplementary Table 3"
        # Add Ensembl NAT links if any were found
        for nat in nat_list:
            details[nat_label(nat)] = f"https://www.ensembl.org/Homo_sapiens/Gene/Summary?g={nat['id']}"

        # Shaped like every other check ({passed, original_passed, overridden}) - the frontend
        # reads status.passed, so a bare bool always rendered as failed
        checks = {
            "ASO already exists": {"passed": False, "original_passed": False, "overridden": False},
            "Overlapping antisense transcript found": {"passed": bool(nat_list), "original_passed": bool(nat_list), "overridden": False},
        }

        if curated_features:
            reason = f"Curated evidence indicates available: {', '.join(curated_features)} for {gene_symbol}."
            if nat_names:
                reason += f" Additionally, Ensembl found antisense transcript(s): {nat_names}."
            return {
                "classification": "Potential possibilities identified",
                "reason": reason,
                "details": details,
                "antisense_gene_ids": nat_ids,
                "checks": checks
            }
        else:
            reason = "No curated evidence found in the N1C supplementary table for this gene."
            if nat_names:
                reason += f" However, Ensembl reports antisense transcript(s): {nat_names}."
            return {
                "classification": "No potential possibilities identified",
                "reason": reason,
                "details": details,
                "antisense_gene_ids": nat_ids,
                "checks": checks
            }

    except Exception:
        traceback.print_exc()
        return {
            "classification": "Unable to Assess",
            "reason": "An error occurred while searching for antisense transcripts.",
            "checks": {}
        }

def _with_aso_exists_check(result: Optional[Dict[str, Any]], overrides: Dict[str, bool]) -> Optional[Dict[str, Any]]:
    """Adds the user-overridable 'ASO already exists' check to a splice-switching result that lacks one."""
    if result is not None:
        result.setdefault('checks', {}).setdefault('ASO already exists', {
            "passed": overrides.get('ASO already exists', False),
            "original_passed": False,
            "overridden": 'ASO already exists' in overrides
        })
    return result

def assess_splice_switching(variant_hgvs: str, vep_data: Dict[str, Any], gene_symbol: str, exonic_pathogenic_user_input: Optional[str] = None, overrides: Optional[Dict[str, bool]] = None, target_consequence: Optional[Dict[str, Any]] = None) -> Optional[Dict[str, Any]]:
    """
    Assesses a variant for splice-switching potential, adding method and DOI link.
    If not found in the DB, it returns a prompt for user validation.
    """
    overrides = overrides or {}
    if splicevar_df is None or not variant_hgvs or not gene_symbol:
        return None

    aso_exists_final = overrides.get('ASO already exists', False)
    if aso_exists_final:
        return {
            "classification": "Eligible",
            "reason": "User confirms an ASO for this variant already exists, therefore it is eligible for splice switching.",
            "checks": {
                "ASO already exists": {"passed": True, "original_passed": False, "overridden": True}
            }
        }

    # Build genomic key from VEP data (matches SpliceVarDB 'hg38' column: chrom-pos-ref-alt;
    # both sides are upper-case, the column already at load time)
    genomic_key = _build_hg38_key_from_vep(vep_data) if vep_data else None
    # --- Primary lookup: match by genomic coordinates (hg38) ---
    matched_rows = pd.DataFrame()
    if genomic_key and 'hg38' in splicevar_df.columns:
        matched_rows = splicevar_df[splicevar_df['hg38'] == genomic_key]

    if not matched_rows.empty:
        for _, row in matched_rows.iterrows():
            method = str(row.get('method', 'N/A')).strip()
            classification = str(row.get('classification', 'N/A')).strip().lower()
            # 'method' can list several assays ('Minigene Assay,RT-PCR'); any RNA-level one counts
            has_valid_method = any(m.strip() in ('rna-seq', 'rt-pcr') for m in method.lower().split(','))
            is_splice_altering = classification == 'splice-altering'
            is_normal = classification == 'normal'
            is_conflicting_or_low_freq = classification in ('conflicting', 'low-frequency')

            details = {"Confirmation Method": method}
            genomic_hgvs = str(row.get('hgvs_genomic', '')).strip()
            if genomic_hgvs and genomic_hgvs.lower() not in ['na', 'n/a', 'nan', '']:
                details["Genomic HGVS"] = genomic_hgvs
            doi = str(row.get('doi', '')).strip()
            if doi and doi.lower() not in ['na', 'n/a', 'nan', '']:
                details["Publication"] = f"https://doi.org/{doi}"

            if has_valid_method and is_conflicting_or_low_freq:
                return {
                    "classification": "Unable to Assess",
                    "reason": "Experimental evidence has conflicting evidence or low frequency reports of splicing changes for this variant.",
                    "details": details,
                    "checks": {
                        "ASO already exists": {"passed": False, "original_passed": False, "overridden": False}
                    }
                }

            if has_valid_method and is_normal:
                return {
                    "classification": "Not Eligible",
                    "reason": "Experimental evidence confirms no splice change for this variant.",
                    "details": details,
                    "checks": {
                        "ASO already exists": {"passed": False, "original_passed": False, "overridden": False}
                    }
                }

            if not (has_valid_method and is_splice_altering):
                continue

            return _with_aso_exists_check(_evaluate_splice_variant_position(variant_hgvs, vep_data, details, exonic_pathogenic_user_input=exonic_pathogenic_user_input, overrides=overrides, target_consequence=target_consequence), overrides)

    # --- Fallback: check SSCVDB (Variant_ID upper-cased at load time, like the key) ---
    if sscvdb_df is not None and not sscvdb_df.empty and vep_data:
        variant_key = _format_sscvdb_variant_id_from_vep(vep_data)
        if variant_key and 'Variant_ID' in sscvdb_df.columns:
            if (sscvdb_df['Variant_ID'] == variant_key).any():
                details = {
                    "Source Database": "SSCVDB",
                    "Evidence": "Splice-altering reported in SSCVDB",
                    "SSCVDB Gene Page": f"https://sscvdb.io/gene/{gene_symbol}"
                }
                return _with_aso_exists_check(_evaluate_splice_variant_position(variant_hgvs, vep_data, details, exonic_pathogenic_user_input=exonic_pathogenic_user_input, overrides=overrides, target_consequence=target_consequence), overrides)

    # --- Not found in either database ---
    return {
        "classification": "Not in Database",
        "reason": "This variant was not found in SpliceVarDB or SSCVDB. If there is experimental validation (by qPCR or RNA-seq) that this variant is splice-altering, please confirm below.",
        "user_validation_prompt": True,
        "checks": {
            "ASO already exists": {"passed": False, "original_passed": False, "overridden": False}
        }
    }

PAIR_SPLICE_CHECK = 'Exons do not contain pathogenic splice variants leading to exon skipping (user decision required)'

def _cds_problem(cds_seq: Optional[str], total_cds_len: int) -> Optional[str]:
    """Why the exon-skipping sequence checks can't run on this CDS, or None if they can. The
    exon structure comes from Ensembl REST and the CDS from refget; the per-exon slices are only
    right if the CDS is exactly as long as the exons' coding parts together (it would not be if
    the two services ever served different Ensembl releases)."""
    if not cds_seq:
        return "The transcript's coding sequence could not be retrieved from Ensembl, so the premature stop codon check cannot be performed. Please try again later."
    if len(cds_seq) != total_cds_len:
        return (f"The transcript's coding sequence from Ensembl ({len(cds_seq)} nt) does not match its exon structure "
                f"({total_cds_len} nt), so the premature stop codon check cannot be performed reliably. Please try again later.")
    return None

def _exon_variant_window(exon: Dict[str, Any]) -> Tuple[str, int, int]:
    """(chrom, start, end) searched for an exon's variants: the exon widened by SPLICE_FLANK_BP
    intronic bases on each side, so its canonical splice-site variants are counted too."""
    return exon['seq_region_name'], max(1, exon['start'] - SPLICE_FLANK_BP), exon['end'] + SPLICE_FLANK_BP

def assess_consecutive_exons(client, transcript, all_exons, target_exon, overrides=None, gene_symbol=None):
    overrides = overrides or {}
    coding_exons = sorted([e for e in all_exons if e['cds_length'] > 0], key=lambda x: x['coding_exon_number'])
    
    target_idx = next((i for i, ex in enumerate(coding_exons) if ex['exon_id'] == target_exon['exon_id']), -1)
    if target_idx == -1: return None

    # Helper to assess a specific pair of exons
    def evaluate_pair(exon_a, exon_b):
        gene_id = transcript.get('Parent')
        transcript_id = transcript.get('id')
        protein_id = transcript.get("Translation", {}).get("id")
        protein_meta = client.get_protein_meta(protein_id) if protein_id else None
        
        total_cds_len = sum(e['cds_length'] for e in coding_exons)
        total_protein_len = protein_meta['length'] if (protein_meta and 'length' in protein_meta) else (total_cds_len / 3 if total_cds_len > 0 else 0)
        
        combined_cds_len = exon_a['cds_length'] + exon_b['cds_length']
        combined_aa_len = combined_cds_len / 3
        combined_len_nt = (exon_a['end'] - exon_a['start'] + 1) + (exon_b['end'] - exon_b['start'] + 1)
        
        # Stop codon check
        cds_seq = client.get_cds_sequence(transcript_id)
        cds_problem = _cds_problem(cds_seq, total_cds_len)
        if cds_problem:
            return {
                "classification": "Unable to Assess",
                "reason": cds_problem,
                "note": f"Assessment based on skipping two consecutive exons: Exon {exon_a['total_exon_number']} and Exon {exon_b['total_exon_number']}.",
                "coding_exon_number": f"{exon_a['coding_exon_number']}+{exon_b['coding_exon_number']}",
                "total_exon_number": f"{exon_a['total_exon_number']}+{exon_b['total_exon_number']} / {len(all_exons)}",
                "gene_id": gene_id, "transcript_id": transcript_id,
            }
        orig_cond2_no_stop = False
        orig_cond10_no_new_codon = True
        try:
            current_pos = 0
            cds_map = {}
            for ex in coding_exons:
                cds_map[ex['coding_exon_number']] = cds_seq[current_pos : current_pos + ex['cds_length']]
                current_pos += ex['cds_length']
            skipped_cds = "".join(cds_map[i] for i in sorted(cds_map.keys()) if i != exon_a['coding_exon_number'] and i != exon_b['coding_exon_number'])
            if skipped_cds:
                prot = str(Seq(skipped_cds).translate(to_stop=False))
                orig_cond2_no_stop = "*" not in prot[:-1]

            # Check for novel codon creation at junction
            idx_a = next((i for i, ex in enumerate(coding_exons) if ex['coding_exon_number'] == exon_a['coding_exon_number']), -1)

            if idx_a > 0 and idx_a + 1 < len(coding_exons) - 1:
                len_prev = sum(e['cds_length'] for i, e in enumerate(coding_exons) if i < idx_a)
                phase = len_prev % 3

                if phase != 0:
                    prev_seq = cds_map[coding_exons[idx_a-1]['coding_exon_number']]
                    next_seq = cds_map[coding_exons[idx_a+2]['coding_exon_number']]
                    target_seq = cds_map[exon_a['coding_exon_number']]
                    needed = 3 - phase
                    if len(prev_seq) >= phase and len(next_seq) >= needed and len(target_seq) >= needed:
                        new_codon = prev_seq[-phase:] + next_seq[:needed]
                        orig_codon = prev_seq[-phase:] + target_seq[:needed]
                        if str(Seq(new_codon).translate()) != str(Seq(orig_codon).translate()):
                            orig_cond10_no_new_codon = False
        except Exception: orig_cond2_no_stop = False

        orig_cond3_not_terminal = (exon_a['coding_exon_number'] != 1 and exon_b['coding_exon_number'] != len(coding_exons))
        orig_cond4_small = (Decimal(str(combined_aa_len)) / Decimal(str(total_protein_len))) < Decimal('0.1') if total_protein_len > 0 else False
        
        domains = client.get_domains(protein_id) if protein_id else []
        overlapping_domain_names, overlapping_domain_details = [], []
        if domains:
            cds_pos_start = sum(e['cds_length'] for e in coding_exons if e['coding_exon_number'] < exon_a['coding_exon_number'])
            exon_aa_start = (cds_pos_start // 3) + 1
            exon_aa_end = ((cds_pos_start + combined_cds_len - 1) // 3) + 1
            for d in domains:
                if d.get('start', 0) <= exon_aa_end and d.get('end', 0) >= exon_aa_start:
                    label = d.get('display_name') or d.get('description', d.get('id', 'Unknown Domain'))
                    overlapping_domain_names.append(label)
                    overlapping_domain_details.append({
                        "label": label, "name": d.get('description') or d.get('id') or "Domain",
                        "source": d.get('source') or d.get('type'), "id": d.get('id'), "url": d.get('url'),
                        "start": d.get('start'), "end": d.get('end')
                    })
        domain_count = len(overlapping_domain_names)
        orig_cond5_no_domain = domain_count == 0

        trexome_repeats_a = get_trexome_repeats_for_exon(gene_symbol, exon_a.get('seq_region_name'), exon_a.get('start'), exon_a.get('end'))
        trexome_repeats_b = get_trexome_repeats_for_exon(gene_symbol, exon_b.get('seq_region_name'), exon_b.get('start'), exon_b.get('end'))
        seen_trexome_keys = set()
        trexome_repeats = []
        for r in trexome_repeats_a + trexome_repeats_b:
            key = (r.get('uniprot_id'), r.get('protein_start'), r.get('protein_end'))
            if key in seen_trexome_keys:
                continue
            seen_trexome_keys.add(key)
            trexome_repeats.append(r)

        # Variants
        def _get_counts_for_exon(exon):
            # Variant window: the exon plus its canonical splice sites (see SPLICE_FLANK_BP)
            chrom, start, end = _exon_variant_window(exon)
            vars_in_reg = client.overlap_region_variation(chrom, start, end)
            clinvar_reg = fetch_clinvar_region_cached(chrom, start, end)
            return count_variants_in_region(vars_in_reg, clinvar_reg)

        _fut_a = _IO_EXECUTOR.submit(_get_counts_for_exon, exon_a)
        _fut_b = _IO_EXECUTOR.submit(_get_counts_for_exon, exon_b)
        c1, pv1, l1, sfa1 = _fut_a.result()
        c2, pv2, l2, sfa2 = _fut_b.result()
        counts = {k: c1.get(k, 0) + c2.get(k, 0) for k in c1}
        pathogenic_variant_links = {k: pv1[k] + pv2[k] for k in pv1}
        splice_variant_links = {"pathogenic": l1["pathogenic"] + l2["pathogenic"], "benign": l1["benign"] + l2["benign"]}
        star_filter_active = sfa1 and sfa2

        orig_cond6_missense = True  # always True by default; only False when user overrides
        orig_cond7_splice = True  # as for single exons: passes by default, only fails when the user decides so
        orig_cond8_no_inframe_del = counts['inframe_del'] == 0
        orig_cond_important_domain = True  # always True by default; only False when user overrides

        # Apply Overrides
        cond2_no_stop = overrides.get('No New Stop Codon', orig_cond2_no_stop)
        cond3_not_terminal = overrides.get('Not First/Last Exon', orig_cond3_not_terminal)
        cond4_small = overrides.get('Is <10% of Protein', orig_cond4_small)
        cond5_no_domain = overrides.get('No Domain Overlap', orig_cond5_no_domain)
        cond6_missense = overrides.get('exon does not contain a mutational hotspot (user input required)', orig_cond6_missense)
        cond7_splice = overrides.get(PAIR_SPLICE_CHECK, orig_cond7_splice)
        cond8_no_inframe_del = overrides.get('No Pathogenic In-Frame Deletions', orig_cond8_no_inframe_del)
        cond10_no_new_codon = overrides.get('No New Codon Formation', orig_cond10_no_new_codon)
        cond_important_domain = overrides.get('Important Domain', orig_cond_important_domain)

        classification, reason = "Undetermined", ""
        if not cond3_not_terminal: classification, reason = "Not Eligible", "One of the exons is the first or last coding exon."
        elif not cond2_no_stop: classification, reason = "Not Eligible", "Skipping these exons is predicted to create a premature stop codon."
        elif not cond7_splice: classification, reason = "Not Eligible", f"Exons contain {counts['splice']} pathogenic splice variant(s), User confirmed the splice site is functionally important and leads to pathogenic exon skipping."
        elif not cond4_small and domain_count > 1: classification, reason = "Not Eligible", f"Exons overlap with {domain_count} protein domains and constitute >10% of the coding region."
        elif not cond8_no_inframe_del: classification, reason = "Not Eligible", f"Exons contain {counts['inframe_del']} pathogenic in-frame deletion(s)."
        elif not cond_important_domain: classification, reason = "Not Eligible", "Exons harbor an important protein domain as indicated by user."
        elif not cond10_no_new_codon: classification, reason = "Unlikely Eligible", "Skipping these exons creates a novel codon at the junction that codes for a different amino acid."
        elif not cond6_missense: classification, reason = "Unlikely Eligible", f"Exons are a mutational hotspot with {counts['missense']} pathogenic missense variants."
        elif not cond4_small: classification, reason = "Unlikely Eligible", "Exons constitute >=10% of the protein."
        elif not cond5_no_domain: classification, reason = "Unlikely Eligible", f"Exons overlap with {domain_count} protein domain(s)."
        else: classification, reason = "Likely Eligible", "Exons meet the primary criteria for multi-exon skipping."
        
        note = f"Assessment based on skipping two consecutive exons: Exon {exon_a['total_exon_number']} and Exon {exon_b['total_exon_number']}."
        if combined_len_nt > 300:
            note += " Combined exon length exceeds 300 nucleotides. Proceed with caution."
        if trexome_repeats:
            repeat_types = ', '.join(sorted({r['repeat_type'] for r in trexome_repeats if r.get('repeat_type')}))
            note += f" These exons overlap {len(trexome_repeats)} tandem repeat domain(s) ({repeat_types}) annotated in TReXome DB."
        note += " Please note that some exceptions for exon skipping assessment apply. Please refer to the methods section for more details."

        return {
            "classification": classification, "reason": reason, "note": note, "exon_length": combined_len_nt,
            "frac_cds": f"{((Decimal(str(combined_aa_len)) / Decimal(str(total_protein_len))) * 100):.2f}%" if total_protein_len > 0 else "N/A",
            "pathogenic_variant_counts": counts, "domain_count": domain_count, "domain_names": overlapping_domain_names,
            "domain_details": overlapping_domain_details, "splice_variant_links": splice_variant_links, "pathogenic_variant_links": pathogenic_variant_links,
            "trexome_repeats": trexome_repeats,
            "coding_exon_number": f"{exon_a['coding_exon_number']}+{exon_b['coding_exon_number']}",
            "total_exon_number": f"{exon_a['total_exon_number']}+{exon_b['total_exon_number']} / {len(all_exons)}",
            "gene_id": gene_id, "transcript_id": transcript_id,
            "clinvar_url": f"https://www.ncbi.nlm.nih.gov/clinvar/?term=GRCh38%3A{exon_a['seq_region_name']}%3A{min(exon_a['start'], exon_b['start']) - SPLICE_FLANK_BP}-{max(exon_a['end'], exon_b['end']) + SPLICE_FLANK_BP}",
            "checks": {
                "Is In-Frame": {"passed": True, "original_passed": True, "overridden": 'Is In-Frame' in overrides},
                "No New Stop Codon": {"passed": cond2_no_stop, "original_passed": orig_cond2_no_stop, "overridden": 'No New Stop Codon' in overrides},
                "Not First/Last Exon": {"passed": cond3_not_terminal, "original_passed": orig_cond3_not_terminal, "overridden": 'Not First/Last Exon' in overrides},
                PAIR_SPLICE_CHECK: {"passed": cond7_splice, "original_passed": orig_cond7_splice, "overridden": PAIR_SPLICE_CHECK in overrides},
                "No Pathogenic In-Frame Deletions": {"passed": cond8_no_inframe_del, "original_passed": orig_cond8_no_inframe_del, "overridden": 'No Pathogenic In-Frame Deletions' in overrides},
                "Important Domain": {"passed": cond_important_domain, "original_passed": orig_cond_important_domain, "overridden": 'Important Domain' in overrides},
                "No Domain Overlap": {"passed": cond5_no_domain, "original_passed": orig_cond5_no_domain, "overridden": 'No Domain Overlap' in overrides},
                "No New Codon Formation": {"passed": cond10_no_new_codon, "original_passed": orig_cond10_no_new_codon, "overridden": 'No New Codon Formation' in overrides},
                "exon does not contain a mutational hotspot (user input required)": {"passed": cond6_missense, "original_passed": orig_cond6_missense, "overridden": 'exon does not contain a mutational hotspot (user input required)' in overrides},
                "Is <10% of Protein": {"passed": cond4_small, "original_passed": orig_cond4_small, "overridden": 'Is <10% of Protein' in overrides}
            },
            "min_review_stars_applied": MIN_REVIEW_STARS if star_filter_active else None,
        }

    candidates = []
    if target_idx > 0:
        prev = coding_exons[target_idx - 1]
        if (prev['cds_length'] + target_exon['cds_length']) % 3 == 0:
            candidates.append((prev, target_exon))
    if target_idx < len(coding_exons) - 1:
        nxt = coding_exons[target_idx + 1]
        if (target_exon['cds_length'] + nxt['cds_length']) % 3 == 0:
            candidates.append((target_exon, nxt))
            
    if not candidates: return None

    # Each candidate pair is an independent Ensembl/ClinVar lookup - evaluate them
    # concurrently instead of one after another.
    pair_futures = [_PAIR_EXECUTOR.submit(evaluate_pair, ex_a, ex_b) for ex_a, ex_b in candidates]
    return [f.result() for f in pair_futures]

def assess_single_exon(client, transcript, all_exons, target_exon, vep_entry: Dict[str, Any], overrides: Optional[Dict[str, bool]] = None, aso_exists: bool = False, moa_user_input: Optional[str] = None, gene_symbol: Optional[str] = None):
    overrides = overrides or {}
    # --- Step 1: Data Gathering and Calculations ---
    # Set precision for Decimal calculations
    getcontext().prec = 10

    gene_id = transcript.get('Parent') 
    transcript_id = transcript.get('id')
    protein_id = transcript.get("Translation", {}).get("id")
    protein_meta = client.get_protein_meta(protein_id) if protein_id else None
    
    cds_seq = client.get_cds_sequence(transcript_id)
    
    coding_exons = [e for e in all_exons if e['cds_length'] > 0]
    total_coding_exons = len(coding_exons)
    total_cds_len = sum(e['cds_length'] for e in coding_exons)
    if protein_meta and 'length' in protein_meta:
        total_protein_len = protein_meta['length']
    else:
        total_protein_len = total_cds_len / 3 if total_cds_len > 0 else 0

    if not target_exon.get('coding_exon_number'):
        return [{"classification": "Unable to Assess", "reason": f"The variant maps to exon {target_exon['total_exon_number']}, which is non-coding."}]
    # Without a CDS that fits the exon structure the stop-codon / junction-codon checks can't run.
    # Say so rather than letting cond2 default to False, which would read as "premature stop
    # codon created".
    cds_problem = _cds_problem(cds_seq, total_cds_len)
    if cds_problem:
        return [{"classification": "Unable to Assess", "reason": cds_problem}]

    # Variant window: the exon plus its canonical splice sites (see SPLICE_FLANK_BP)
    chrom, start, end = _exon_variant_window(target_exon)
    variants_in_region = client.overlap_region_variation(chrom, start, end)
    clinvar_region = fetch_clinvar_region_cached(chrom, start, end)

    clinvar_url = f"https://www.ncbi.nlm.nih.gov/clinvar/?term=GRCh38%3A{chrom}%3A{start}-{end}"

    counts, pathogenic_variant_links, splice_variant_links, star_filter_active = count_variants_in_region(variants_in_region, clinvar_region)
    
    exon_cds_len = target_exon['cds_length']
    exon_aa_len = exon_cds_len / 3
    coding_exon_number = target_exon['coding_exon_number']

    # --- Step 2: Condition Checks ---
    # These are the automatically calculated states
    original_conds = {}
    original_conds['cond1_inframe'] = (exon_cds_len % 3 == 0)
    original_conds['cond2_no_stop'] = False
    original_conds['cond10_no_new_codon'] = True
    try:
        cds_map, current_pos = {}, 0
        sorted_coding_exons = sorted(coding_exons, key=lambda x: x['coding_exon_number'])
        for exon in sorted_coding_exons:
            cds_map[exon['coding_exon_number']] = cds_seq[current_pos : current_pos + exon['cds_length']]
            current_pos += exon['cds_length']
        skipped_cds = "".join(cds_map[i] for i in sorted(cds_map.keys()) if i != coding_exon_number)
        if skipped_cds:
            prot = str(Seq(skipped_cds).translate(to_stop=False))
            original_conds['cond2_no_stop'] = "*" not in prot[:-1]

        # Check for novel codon creation at junction
        target_idx = next((i for i, ex in enumerate(sorted_coding_exons) if ex['coding_exon_number'] == coding_exon_number), -1)
        if original_conds['cond1_inframe'] and target_idx > 0 and target_idx < len(sorted_coding_exons) - 1:
            len_prev = sum(e['cds_length'] for i, e in enumerate(sorted_coding_exons) if i < target_idx)
            phase = len_prev % 3
            if phase != 0:
                prev_seq = cds_map[sorted_coding_exons[target_idx-1]['coding_exon_number']]
                next_seq = cds_map[sorted_coding_exons[target_idx+1]['coding_exon_number']]
                target_seq = cds_map[coding_exon_number]
                needed = 3 - phase
                if len(prev_seq) >= phase and len(next_seq) >= needed and len(target_seq) >= needed:
                    new_codon = prev_seq[-phase:] + next_seq[:needed]
                    orig_codon = prev_seq[-phase:] + target_seq[:needed]
                    if str(Seq(new_codon).translate()) != str(Seq(orig_codon).translate()):
                        original_conds['cond10_no_new_codon'] = False
    except Exception:
        original_conds['cond2_no_stop'] = False

    original_conds['cond3_not_terminal'] = (coding_exon_number is not None and coding_exon_number not in (1, total_coding_exons))
    original_conds['cond4_small'] = (Decimal(str(exon_aa_len)) / Decimal(str(total_protein_len))) < Decimal('0.1') if total_protein_len > 0 else False
    domains = client.get_domains(protein_id) if protein_id else []
    overlapping_domain_names, overlapping_domain_details = [], []
    if domains:
        cds_pos_start = sum(e['cds_length'] for e in sorted(coding_exons, key=lambda x: x['coding_exon_number']) if e['coding_exon_number'] < coding_exon_number)
        exon_aa_start, exon_aa_end = (cds_pos_start // 3) + 1, ((cds_pos_start + exon_cds_len -1) // 3) + 1
        for d in domains:
            if d.get('start', 0) <= exon_aa_end and d.get('end', 0) >= exon_aa_start:
                label = d.get('display_name') or d.get('description', d.get('id', 'Unknown Domain'))
                overlapping_domain_names.append(label)
                overlapping_domain_details.append({
                    "label": label,
                    "name": d.get('description') or d.get('id') or d.get('hit_id') or "Domain",
                    "source": d.get('source') or d.get('type'),
                    "id": d.get('id') or d.get('hit_id'),
                    "url": d.get('url'),
                    "start": d.get('start'),
                    "end": d.get('end')
                })

    domain_count = len(overlapping_domain_names)
    original_conds['cond5_no_domain'] = domain_count == 0
    original_conds['cond6_missense'] = True  # always True by default; only False when user overrides
    original_conds['cond7_splice'] = True
    original_conds['cond8_no_inframe_del'] = counts['inframe_del'] == 0
    original_conds['cond9_benign_splice'] = counts['benign_splice'] > 0
    original_conds['aso_exists'] = aso_exists
    original_conds['cond_important_domain'] = True  # always True by default; only False when user overrides

    trexome_repeats = get_trexome_repeats_for_exon(gene_symbol, target_exon.get('seq_region_name'), target_exon.get('start'), target_exon.get('end'))

    # --- Step 2.5: Apply Overrides to get FINAL states ---
    # Start with original values, then overwrite if an override exists
    cond1_inframe = overrides.get('Exon is in-frame', original_conds['cond1_inframe'])
    cond2_no_stop = overrides.get('Skipping exon does not create new stop codon', original_conds['cond2_no_stop'])
    cond3_not_terminal = overrides.get('Exon is not the first or last coding exon', original_conds['cond3_not_terminal'])
    cond4_small = overrides.get('Exon encodes for <10% of Protein', original_conds['cond4_small'])
    cond5_no_domain = overrides.get('Exon does not overlap with a domain', original_conds['cond5_no_domain'])
    cond6_missense = overrides.get('exon does not contain a mutational hotspot (user input required)', original_conds['cond6_missense'])
    cond7_splice = overrides.get('Exon does not contain pathogenic splice variants leading to exon skipping(user decision required)', original_conds['cond7_splice'])
    cond8_no_inframe_del = overrides.get('Exon does not contain pathogenic in-frame deletions', original_conds['cond8_no_inframe_del'])
    cond9_benign_splice = overrides.get('Benign splice variant leading to exon skipping found within exon', original_conds['cond9_benign_splice'])
    cond10_no_new_codon = overrides.get('No New Codon Formation', original_conds['cond10_no_new_codon'])
    aso_exists_final = overrides.get('ASO targeting the exon already exists', original_conds['aso_exists'])
    cond_important_domain = overrides.get('Exon does not contain an important domain (user decision required)', original_conds['cond_important_domain'])
    has_override = bool(overrides)

    # Check if any override actually changed the outcome
    is_changed = False
    if has_override:
        check_map = {
            'Exon is in-frame': cond1_inframe != original_conds['cond1_inframe'], 'Skipping exon does not create new stop codon': cond2_no_stop != original_conds['cond2_no_stop'],
            'Exon is not the first or last coding exon': cond3_not_terminal != original_conds['cond3_not_terminal'], 'Exon encodes for <10% of Protein': cond4_small != original_conds['cond4_small'],
            'Exon does not overlap with a domain': cond5_no_domain != original_conds['cond5_no_domain'], 'exon does not contain a mutational hotspot (user input required)': cond6_missense != original_conds['cond6_missense'],
            'Exon does not contain pathogenic splice variants leading to exon skipping(user decision required)': cond7_splice != original_conds['cond7_splice'], 'Exon does not contain pathogenic in-frame deletions': cond8_no_inframe_del != original_conds['cond8_no_inframe_del'],
            'Benign splice variant leading to exon skipping found within exon': cond9_benign_splice != original_conds['cond9_benign_splice'],
            'No New Codon Formation': cond10_no_new_codon != original_conds['cond10_no_new_codon'],
            'ASO targeting the exon already exists': aso_exists_final != original_conds['aso_exists'],
            'Exon does not contain an important domain (user decision required)': cond_important_domain != original_conds['cond_important_domain']
        }
        is_changed = any(check_map.get(check_name, False) for check_name in overrides)

    # --- Step 3: Classification Logic Chain ---
    if aso_exists_final:
        classification, reason = "Eligible", "An ASO targeting this exon already exists."
    elif not cond3_not_terminal:
        classification, reason = "Not Eligible", "Exon is the first or last coding exon."
    else:
        if not cond1_inframe:
            classification, reason = "Not Eligible", "Exon is out-of-frame, which would disrupt the reading frame. Multi-exon skipping was not found to be a viable in-frame alternative."
        elif not cond2_no_stop:
            classification, reason = "Not Eligible", "Skipping this exon is predicted to create a premature stop codon."
        elif not cond7_splice:
            classification, reason = "Not Eligible", f"Exon contains {counts['splice']} pathogenic splice variant(s), User confirmed the splice site is functionally important and leads to pathogenic exon skipping."
        elif not cond4_small and domain_count > 1:
            classification, reason = "Not Eligible", f"Exon overlaps with {domain_count} protein domains and constitutes >10% of the coding region."
        elif not cond8_no_inframe_del:
            classification, reason = "Not Eligible", f"Exon contains {counts['inframe_del']} pathogenic in-frame deletion(s)."
        elif not cond_important_domain:
            classification, reason = "Not Eligible", "Exon harbors an important protein domain as indicated by user."
        elif cond9_benign_splice:
            if counts['splice'] > 0 or counts['inframe_del'] > 0:
                classification, reason = "Unlikely Eligible", "Exon has benign splice variants but also pathogenic splice or in-frame deletion variants, suggesting caution."
            else:
                classification, reason = "Eligible", "Exon contains benign splice variants, suggesting it may be safely skipped."
        elif not cond10_no_new_codon:
            classification, reason = "Unlikely Eligible", "Skipping this exon creates a novel codon at the junction that codes for a different amino acid."
        elif not cond6_missense:
            classification, reason = "Unlikely Eligible", f"Exon is a mutational hotspot with {counts['missense']} pathogenic missense variants."
        elif not cond4_small:
            classification, reason = "Unlikely Eligible", "Exon constitutes >=10% of the protein, risking major functional loss."
        elif not cond5_no_domain:
            classification, reason = "Unlikely Eligible", f"Exon overlaps with {domain_count} protein domain(s)."
        else:
            classification, reason = "Likely Eligible", "Exon meets the primary criteria for a skippable exon."
    
    note_parts = []
    if is_changed:
        note_parts.append("Classification based on manual override of guideline checks.")
    exon_len_nt = target_exon['end'] - target_exon['start'] + 1
    if exon_len_nt > 300:
        note_parts.append("Exon length exceeds 300 nucleotides. Proceed with caution.")
    if trexome_repeats:
        repeat_types = ', '.join(sorted({r['repeat_type'] for r in trexome_repeats if r.get('repeat_type')}))
        note_parts.append(f"Exon overlaps {len(trexome_repeats)} tandem repeat domain(s) ({repeat_types}) annotated in TReXome DB.")
    note_parts.append("Please note that some exceptions for exon skipping apply. Please refer to the methods section for more details.")
    if moa_user_input not in ("GoF", "LoF", "DN"):
        note_parts.append("\n Exon skipping assessment without knowledge of pathomechanism would technically be classified as 'unable to assess'. We strongly recommend to identify the pathomechanism of the variant before continuing")
    
    note = " ".join(note_parts)

    # --- Step 4: Visualization Data Generation ---
    visualization_data = None
    try:
        v_chrom, v_start, v_end = vep_entry.get('seq_region_name'), vep_entry.get('start'), vep_entry.get('end')
        if not all([v_chrom, v_start, v_end]): 
            raise ValueError("Missing variant coordinates for visualization.")
        domain_features = []
        if protein_id and domains:
            cds_map = []
            cumulative_cds_len = 0
            is_reverse_strand = transcript.get('strand') == -1
            for exon in sorted(coding_exons, key=lambda x: x['coding_exon_number']):
                cds_len_of_exon = exon['cds_length']
                cds_map.append({
                    'chr': exon['seq_region_name'], 
                    'genomic_start': exon['start'], 
                    'genomic_end': exon['end'], 
                    'transcript_cds_start': cumulative_cds_len + 1, 
                    'transcript_cds_end': cumulative_cds_len + cds_len_of_exon
                })
                cumulative_cds_len += cds_len_of_exon
            
            for domain in domains:
                domain_cds_start, domain_cds_end = (domain['start'] - 1) * 3 + 1, domain['end'] * 3
                for exon_map_entry in cds_map:
                    overlap_start = max(domain_cds_start, exon_map_entry['transcript_cds_start'])
                    overlap_end = min(domain_cds_end, exon_map_entry['transcript_cds_end'])
                    
                    if overlap_start <= overlap_end:
                        offset_start = overlap_start - exon_map_entry['transcript_cds_start']
                        offset_end = overlap_end - exon_map_entry['transcript_cds_start']
                        
                        if not is_reverse_strand:
                            feat_start = exon_map_entry['genomic_start'] + offset_start
                            feat_end = exon_map_entry['genomic_start'] + offset_end
                        else: # On reverse strand, offsets are from the end
                            feat_start = exon_map_entry['genomic_end'] - offset_end
                            feat_end = exon_map_entry['genomic_end'] - offset_start
                            
                        domain_features.append({
                            "chr": exon_map_entry['chr'], 
                            "start": feat_start - 1, 
                            "end": feat_end, 
                            "name": domain.get('display_name', domain.get('description', domain.get('id', 'Domain')))
                        })

        padding = 1000
        visualization_data = {
            "locus": f"{v_chrom}:{max(1, v_start - padding)}-{v_end + padding}",
            "variantTrack": {"name": "Variant", "features": [{"chr": v_chrom, "start": v_start - 1, "end": v_end, "name": vep_entry.get('id', 'Variant')}]},
            "domainTrack": {"name": "Protein Domains", "features": domain_features} if domain_features else None
        }
    except Exception:
        traceback.print_exc()
        visualization_data = None

    # --- FINAL RETURN STATEMENT ---
    single_res = {
        "classification": classification,
        "reason": reason,
        "exon_length": exon_len_nt,
        "note": note,
        "frac_cds": (
            f"{((Decimal(str(exon_aa_len)) / Decimal(str(total_protein_len))) * 100):.2f}%"
            if total_protein_len > 0
            else "N/A"
        ),
        "pathogenic_variant_counts": counts, 
        "domain_count": domain_count,
        "domain_names": overlapping_domain_names,
        "domain_details": overlapping_domain_details,
        "trexome_repeats": trexome_repeats,
        "splice_variant_links": splice_variant_links,
        "pathogenic_variant_links": pathogenic_variant_links,
        "coding_exon_number": coding_exon_number,
        "total_exon_number": f"{target_exon['total_exon_number']} / {len(all_exons)}",
        "gene_id": gene_id,
        "transcript_id": transcript_id,
        "clinvar_url": clinvar_url,
        "checks": {
            "ASO targeting the exon already exists": {"passed": aso_exists_final, "original_passed": original_conds['aso_exists'], "overridden": 'ASO targeting the exon already exists' in overrides},
            "Benign splice variant leading to exon skipping found within exon": {"passed": cond9_benign_splice, "original_passed": original_conds['cond9_benign_splice'], "overridden": 'Benign splice variant leading to exon skipping found within exon' in overrides},
            "Exon is in-frame": {"passed": cond1_inframe, "original_passed": original_conds['cond1_inframe'], "overridden": 'Exon is in-frame' in overrides},
            "Skipping exon does not create new stop codon": {"passed": cond2_no_stop, "original_passed": original_conds['cond2_no_stop'], "overridden": 'Skipping exon does not create new stop codon' in overrides},
            "Exon is not the first or last coding exon": {"passed": cond3_not_terminal, "original_passed": original_conds['cond3_not_terminal'], "overridden": 'Exon is not the first or last coding exon' in overrides},
            "Exon does not contain pathogenic splice variants leading to exon skipping(user decision required)": {"passed": cond7_splice, "original_passed": original_conds['cond7_splice'], "overridden": 'Exon does not contain pathogenic splice variants leading to exon skipping(user decision required)' in overrides},
            "Exon does not contain pathogenic in-frame deletions": {"passed": cond8_no_inframe_del, "original_passed": original_conds['cond8_no_inframe_del'], "overridden": 'Exon does not contain pathogenic in-frame deletions' in overrides},
            "Exon does not contain an important domain (user decision required)": {"passed": cond_important_domain, "original_passed": original_conds['cond_important_domain'], "overridden": 'Exon does not contain an important domain (user decision required)' in overrides},
            "Exon does not overlap with a domain": {"passed": cond5_no_domain, "original_passed": original_conds['cond5_no_domain'], "overridden": 'Exon does not overlap with a domain' in overrides},
            "No New Codon Formation": {"passed": cond10_no_new_codon, "original_passed": original_conds['cond10_no_new_codon'], "overridden": 'No New Codon Formation' in overrides},
            "exon does not contain a mutational hotspot (user input required)": {"passed": cond6_missense, "original_passed": original_conds['cond6_missense'], "overridden": 'exon does not contain a mutational hotspot (user input required)' in overrides},
            "Exon encodes for <10% of Protein": {"passed": cond4_small, "original_passed": original_conds['cond4_small'], "overridden": 'Exon encodes for <10% of Protein' in overrides}
        },
        "min_review_stars_applied": MIN_REVIEW_STARS if star_filter_active else None,
        "visualization": visualization_data
    }

    results_list = [single_res]

    # Check consecutive exons if originally out of frame, regardless of override status of 'Is In-Frame'
    if not original_conds['cond1_inframe']:
        consecutive_results = assess_consecutive_exons(client, transcript, all_exons, target_exon, overrides, gene_symbol=gene_symbol)
        if consecutive_results:
            results_list.extend(consecutive_results)

    return results_list

def _assess_for_transcript(
    target_consequence: Dict[str, Any],
    vep_entry: Dict[str, Any],
    hgvs_query: str,
    splice_user_input: Optional[str],
    moa_user_input: Optional[str],
    exonic_pathogenic_user_input: Optional[str],
    exon_skipping_overrides: Dict,
    knockdown_overrides: Dict,
    wt_upregulation_overrides: Dict,
    gene_characteristics_overrides: Optional[Dict],
    splice_switching_overrides: Dict,
    client: EnsemblClient
) -> Dict[str, Any]:
    """
    Runs the full ASO-strategy assessment for a single, already-resolved transcript consequence.
    """
    gene_symbol = target_consequence['gene_symbol']
    definitive_transcript_id = target_consequence['transcript_id']
    gene_id = target_consequence.get('gene_id')

    # Only depends on target_consequence, already known - decided early so the exon-
    # skipping transcript lookup below can be kicked off now and overlap with gene-
    # characteristics/knockdown/WT-upregulation/splice-switching work instead of only
    # starting once execution reaches the exon-skipping block.
    consequence_terms = set(target_consequence.get('consequence_terms', []))
    exonic_terms = {'missense_variant', 'stop_gained', 'frameshift_variant', 'synonymous_variant', 'inframe_deletion', 'inframe_insertion','splice_donor_variant', 'splice_acceptor_variant'}
    is_exonic = any(term in consequence_terms for term in exonic_terms)
    transcript_future = _IO_EXECUTOR.submit(client.lookup_id_expand, definitive_transcript_id) if is_exonic else None

    # TReXome only annotates the summary and exon notes, but its host can take tens of seconds
    # to wake up - start its lookups now so they overlap with everything below. The per-exon
    # prefetch warms the cache that assess_single_exon/assess_consecutive_exons read.
    trexome_summary_future = _TREXOME_EXECUTOR.submit(get_trexome_gene_repeat_summary, gene_symbol)
    if is_exonic:
        _TREXOME_EXECUTOR.submit(_trexome_repeats_for_gene, gene_symbol)

    # --- Initialize Result & Gene Characteristics ---
    gene_characteristics = get_gene_characteristics(gene_symbol)
    if gene_characteristics_overrides:
        gene_characteristics.update(gene_characteristics_overrides)

    protein_effect = None

    # Attempt 1: Check the chosen target_consequence directly
    if 'hgvsp' in target_consequence:
        protein_effect = target_consequence['hgvsp']

    # Format HGVSp to include parentheses for predicted effect: ID:p.(Effect)
    if protein_effect and ':p.' in protein_effect and ':p.(' not in protein_effect:
        parts = protein_effect.split(':p.')
        if len(parts) == 2:
            protein_effect = f"{parts[0]}:p.({parts[1]})"

    # Attempt 2: Fallback to Consequence Type
    if not protein_effect:
        terms = target_consequence.get('consequence_terms', [])
        if terms:
            protein_effect = ", ".join(terms).replace('_', ' ').title()
        else:
            protein_effect = "Unable to calculate"

    final_result = {
        "summary": {
            "gene": gene_symbol,
            "transcript_id": definitive_transcript_id,
            "protein_effect": protein_effect,
            "varaico_url": _build_varaico_url(gene_symbol, vep_entry, client),
            "litvar": _build_litvar_info(vep_entry),
            "variant_notations": _build_variant_notations(target_consequence, protein_effect, vep_entry),
            **gene_characteristics
        },
        "assessments": {}
    }

    def finish() -> Dict[str, Any]:
        """Collects the TReXome gene summary started above; used on every return path."""
        try:
            final_result["summary"]["trexome_repeats"] = trexome_summary_future.result()
        except Exception:
            pass  # keep the empty default from get_gene_characteristics
        return final_result

    # Both N1C datasets are keyed by c. notation; for genomic (g.) or other non-c. input, fall
    # back to VEP's cDNA HGVS for this transcript instead of silently never matching
    registry_hgvs = hgvs_query if re.search(r'c\.', hgvs_query, re.IGNORECASE) else (target_consequence.get('hgvsc') or '')

    # --- N1C Assessed Variants (curated) Check (shown first, standard assessment still runs) ---
    assessed_match = check_n1c_assessed_variants(gene_symbol, registry_hgvs)
    if assessed_match:
        final_result["assessments"]["N1C_Assessed_Variants"] = assessed_match

    # --- N1C Registry Check (exit early only for an established ASO project on this variant;
    # under-development and informational entries are shown alongside the assessment) ---
    n1c_result = check_n1c_registry(gene_symbol, registry_hgvs)
    if n1c_result:
        final_result["assessments"]["N1C_Registry_Check"] = n1c_result
        if n1c_result.get("classification") == "Eligible":
            return finish()

    # --- Gene-Level Strategies (Knockdown, WT Upregulation) ---
    moi = set(gene_characteristics.get("moi", []))
    # ClinGen also uses terms that don't say whether the gene acts dominantly or recessively,
    # which the strategy selection below depends on (e.g. plain 'X-linked inheritance')
    unspecific_moi = sorted(m for m in moi if str(m).strip().lower() in UNSPECIFIC_MOI_TERMS)
    if unspecific_moi:
        _add_summary_note(final_result, (
            f"The mode of inheritance is given as {', '.join(unspecific_moi)}, which does not say whether the gene acts "
            f"dominantly or recessively. Gene-level strategies (knockdown, WT upregulation) that depend on it may be missing; "
            f"set the mode of inheritance (AD/AR/XD/XR) to assess them."
        ))

    # Determine resolved MoA: user selection takes precedence; otherwise a single known MoA
    # Force user selection for MoA; do not auto-resolve from database
    resolved_moa = moa_user_input if moa_user_input in ("GoF", "LoF", "DN") else None
    # Expose resolved_moa for frontend rendering
    final_result["summary"]["resolved_moa"] = resolved_moa
    if resolved_moa:
        if resolved_moa in("GoF", "DN"):
            is_ad = any("autosomal dominant" in str(m).lower() or "x-linked dominant" in str(m).lower() for m in moi)
            if is_ad:
                n1c_gene_kd = n1c_gene_knockdown_entry(gene_symbol, resolved_moa)
                if n1c_gene_kd:
                    final_result["assessments"]["N1C_Gene_Knockdown"] = n1c_gene_kd
                    if not n1c_gene_kd.get("under_development"):
                        return finish()

                # Only run regular knockdown if there is no N1C gene-level knockdown entry,
                # or the matching entry is still "Under Development"
                if not n1c_gene_kd or n1c_gene_kd.get("under_development"):
                    knockdown_assessment = assess_knockdown(gene_characteristics, overrides=knockdown_overrides)
                    if resolved_moa == "DN":
                        knockdown_assessment["note"] = "Dominant Negative mechanism; allele-specific knockdown is necessary."

                    if resolved_moa == "GoF":
                        knockdown_assessment["note"] = "Dominant mode of inheritance; allele-specific knockdown is recommended."
                    final_result["assessments"]["Knockdown"] = knockdown_assessment

        if resolved_moa == "GoF":
            # Skipped when the dominant branch above already produced the knockdown assessment
            # (e.g. X-linked dominant matches both), so its note isn't overwritten
            is_xlinked_or_ar = any("x-linked" in str(m).lower() or "autosomal recessive" in str(m).lower() for m in moi)
            if is_xlinked_or_ar and "Knockdown" not in final_result["assessments"]:
                n1c_gene_kd = n1c_gene_knockdown_entry(gene_symbol, resolved_moa)
                if n1c_gene_kd:
                    final_result["assessments"]["N1C_Gene_Knockdown"] = n1c_gene_kd
                    if not n1c_gene_kd.get("under_development"):
                        return finish()

                if not n1c_gene_kd or n1c_gene_kd.get("under_development"):
                    knockdown_assessment = assess_knockdown(gene_characteristics, overrides=knockdown_overrides)
                    final_result["assessments"]["Knockdown"] = knockdown_assessment

        elif resolved_moa == "LoF":
            # Show WT Upregulation only when LoF and Autosomal Dominant MOI (case-insensitive match)
            is_lof_ad = any("autosomal dominant inheritance" in str(m).lower() for m in moi)
            if is_lof_ad:
                final_result["assessments"]["WT_Upregulation"] = assess_wt_upregulation(client, gene_id, gene_symbol, overrides=wt_upregulation_overrides)
            is_x_linked_dominant = any("x-linked dominant inheritance" in str(m).lower() for m in moi)
            if is_x_linked_dominant:
                wt_upregulation_assessment = assess_wt_upregulation(client, gene_id, gene_symbol, overrides=wt_upregulation_overrides)
                wt_upregulation_assessment["note"] = "WT-upregulation only possible if more than one X-Chromosome is available"
                final_result["assessments"]["WT_Upregulation"] = wt_upregulation_assessment
                _add_summary_note(final_result, "WT-upregulation only possible if more than one X-Chromosome is available")

    # --- Variant-Specific Strategies (Splice & Exon Skipping) ---
    # (variant type / is_exonic already resolved above, so transcript_future could be
    # started as early as possible)

    # Run Splice Switching Assessment
    variant_identifier_from_vep = vep_entry.get('input')
    if variant_identifier_from_vep:
        splice_assessment = None
        if splice_user_input == 'yes':
            details = {"Confirmation Method": "User-provided validation (qPCR/RNA-seq/cDNA)"}
            splice_assessment = _with_aso_exists_check(_evaluate_splice_variant_position(variant_identifier_from_vep, vep_entry, details, exonic_pathogenic_user_input=exonic_pathogenic_user_input, overrides=splice_switching_overrides, target_consequence=target_consequence), splice_switching_overrides)
        elif splice_user_input == 'no':
            splice_assessment = {"classification": "Unable to Assess", "reason": "User confirmed no known splice-altering effect."}
        else:
            splice_assessment = assess_splice_switching(variant_identifier_from_vep, vep_entry, gene_symbol, exonic_pathogenic_user_input=exonic_pathogenic_user_input, overrides=splice_switching_overrides, target_consequence=target_consequence)

        if splice_assessment:
            final_result["assessments"]["Splice_Switching"] = splice_assessment

    # Run Exon Skipping Assessment *if* variant is exonic/splice
    if is_exonic and transcript_future is not None:
        exon_skip_assessment_added = False
        transcript_data = transcript_future.result()
        if transcript_data:
            all_exons = extract_exons_from_transcript(transcript_data)
            v_start, v_end = vep_entry['start'], vep_entry['end']
            target_exon = next((ex for ex in all_exons if ex['seq_region_name'] == vep_entry['seq_region_name'] and max(v_start, ex['start']) <= min(v_end, ex['end'])), None)

            if target_exon:
                # N1C registry exon-skipping support: if an N1C exon-skipping ASO project's variant
                # lies in this exon, the exon is marked skippable
                try:
                    n1c_exon_map, n1c_exon_ud_map = n1c_exon_skipping_variant_exon_map(gene_symbol, all_exons, client)
                except Exception:
                    n1c_exon_map, n1c_exon_ud_map = {}, {}

                exon_num = target_exon.get('total_exon_number')
                n1c_exon_match = exon_num is not None and exon_num in n1c_exon_map
                matched_links: List[str] = n1c_exon_map.get(exon_num, []) if n1c_exon_match else []
                # If the only matching N1C entries for this exon are still "Under Development",
                # don't let them short-circuit the regular exon-skipping assessment.
                n1c_exon_under_development = n1c_exon_match and n1c_exon_ud_map.get(exon_num, False)

                exon_skip_results = assess_single_exon(client, transcript_data, all_exons, target_exon, vep_entry, overrides=exon_skipping_overrides, aso_exists=(n1c_exon_match and not n1c_exon_under_development), moa_user_input=moa_user_input, gene_symbol=gene_symbol)

                # Only the single-exon result carries visualization data; lift it to the top level
                visualization = exon_skip_results[0].pop("visualization", None)
                if visualization:
                    final_result["visualization"] = visualization

                if n1c_exon_match:
                    # Add note to the single exon result (first element)
                    target_res = exon_skip_results[0]
                    prev_reason = target_res.get('reason', '')
                    det = target_res.get('details', {}) if isinstance(target_res.get('details'), dict) else {}
                    if n1c_exon_under_development:
                        note = f"N1C registry contains a project that maps to exon {exon_num} in {gene_symbol}, but it is currently listed as 'Under Development'. The regular exon-skipping eligibility assessment for this exon is shown above."
                        det['N1C Registry Status'] = 'Under Development'
                    else:
                        note = f"N1C registry contains a project that maps to exon {exon_num} in {gene_symbol}. The input variant lies in the same exon and is therefore considered skippable."
                        target_res['classification'] = 'Eligible'
                    target_res['reason'] = (prev_reason + ' ' + note).strip()
                    if matched_links: det['N1C Exon Skipping (same exon)'] = ', '.join(sorted(set(matched_links)))
                    det['Supported Exon'] = str(exon_num)
                    target_res['details'] = det

                final_result["assessments"]["Exon_Skipping"] = exon_skip_results
                exon_skip_assessment_added = True

        if not exon_skip_assessment_added:
            # Add this block if the logic above fails to add an assessment
            final_result["assessments"]["Exon_Skipping"] = [{
                "classification": "Unable to Assess",
                "reason": "Variant is exonic or in a splice region, but the target exon could not be determined (e.g., VEP/Ensembl data issue)."
            }]

    # Final Fallback
    # This will now only trigger for non-exonic, non-splice, non-gene-strategy variants.
    if not final_result["assessments"]:
        final_result["assessments"]["General_Assessment"] = {
            "classification": "Unable to Assess",
            "reason": "Could not determine a primary ASO strategy. The variant is not exonic and no gene-level strategies (Knockdown, WT Upregulation) are applicable."
        }

    return finish()


def _build_transcript_info(
    requested_transcript_id: Optional[str],
    input_consequence: Optional[Dict[str, Any]],
    default_consequence: Optional[Dict[str, Any]],
    available_alternate: bool
) -> Optional[Dict[str, Any]]:
    """
    Builds metadata describing the user's requested transcript vs. the MANE Select/canonical
    alternative, for rendering a warning + switch control in the UI. Returns None when the
    user didn't request a specific transcript (plain gene query - default behavior unchanged).
    """
    if not requested_transcript_id:
        return None

    if default_consequence and default_consequence.get('mane_select'):
        fallback_label = "MANE Select"
    elif default_consequence and default_consequence.get('canonical') == 1:
        fallback_label = "Ensembl Canonical"
    else:
        fallback_label = "Most Relevant Transcript"

    if not input_consequence:
        # Ensembl's VEP response didn't cross-reference the requested transcript to any of
        # the transcripts overlapping this variant (this can happen for non-MANE RefSeq
        # transcripts, which VEP doesn't always expose a RefSeq<->Ensembl mapping for).
        default_id = default_consequence.get('transcript_id') if default_consequence else None
        return {
            "selected": "default",
            "requested_transcript_id": requested_transcript_id,
            "input_transcript_id": None,
            "is_input_mane": None,
            "default_transcript_id": default_id,
            "default_label": fallback_label,
            "available_alternate": False,
            "warning": f"Could not confirm transcript '{escape(requested_transcript_id)}' among the transcripts overlapping this variant; showing the {fallback_label} transcript ({default_id}) instead.",
            "version_note": None,
        }

    is_input_mane = bool(input_consequence.get('mane_select'))
    default_id = default_consequence.get('transcript_id') if default_consequence else None
    default_label = fallback_label

    warning = None
    if not is_input_mane:
        warning = "This is not the MANE Select transcript."
        if default_id:
            warning += f" {default_label}: {default_id}."

    version_note = None
    input_transcript_id = input_consequence.get('transcript_id')
    if input_transcript_id and requested_transcript_id.strip().upper() != input_transcript_id.strip().upper():
        version_note = f"You entered {escape(requested_transcript_id)}, but Ensembl's current annotation uses {input_transcript_id}."

    return {
        "selected": "input",
        "requested_transcript_id": requested_transcript_id,
        "input_transcript_id": input_transcript_id,
        "is_input_mane": is_input_mane,
        "default_transcript_id": default_id,
        "default_label": default_label,
        "available_alternate": available_alternate,
        "warning": warning,
        "version_note": version_note,
    }


def process_single_variant(query: str, client: EnsemblClient, splice_user_input: Optional[str] = None, moa_user_input: Optional[str] = None, exonic_pathogenic_user_input: Optional[str] = None, exon_skipping_overrides: Optional[Dict] = None, knockdown_overrides: Optional[Dict] = None, wt_upregulation_overrides: Optional[Dict] = None, gene_characteristics_overrides: Optional[Dict] = None, splice_switching_overrides: Optional[Dict] = None) -> Dict[str, Any]:
    """
    Contains the complete assessment logic for a single variant query.
    """
    exon_skipping_overrides = exon_skipping_overrides or {}
    knockdown_overrides = knockdown_overrides or {}
    wt_upregulation_overrides = wt_upregulation_overrides or {}
    splice_switching_overrides = splice_switching_overrides or {}
    _refresh_n1c_if_stale()
    try:
        # --- 1. VEP and Consequence Selection ---
        hgvs_query, gene_symbol_from_query, requested_transcript_id = parse_hgvs_query(query)
        if not hgvs_query:
            return {"classification": "Error", "reason": "Invalid input format. Please use a recognized HGVS format (e.g., 'GENE c.123A>G')."}

        vep_data = client.vep_hgvs(hgvs_query)
        if not vep_data or not isinstance(vep_data, list):
            return {"classification": "Unable to Assess", "reason": f"VEP analysis failed for '{escape(hgvs_query)}'. The variant may be invalid or not found."}

        vep_entry = vep_data[0]
        all_consequences = vep_entry.get('transcript_consequences', [])

        # Default pick (MANE Select -> canonical -> longest coding), and the consequence
        # matching the transcript the user explicitly typed, if any.
        default_consequence = choose_best_consequence(all_consequences, gene_symbol_from_query=gene_symbol_from_query)
        input_consequence = _find_consequence_by_transcript(all_consequences, requested_transcript_id) if requested_transcript_id else None

        # The user's input transcript drives the assessment by default; fall back to the
        # MANE Select/canonical pick when no transcript was explicitly requested, or when the
        # requested one couldn't be cross-referenced among VEP's transcript consequences.
        target_consequence = input_consequence or default_consequence

        if not target_consequence:
            reason = "VEP did not return a consequence."
            if gene_symbol_from_query:
                reason += f" No valid consequence was found for the specified gene '{escape(gene_symbol_from_query)}'."
            return {"classification": "Unable to Assess", "reason": reason}

        common_args = dict(
            vep_entry=vep_entry, hgvs_query=hgvs_query, splice_user_input=splice_user_input,
            moa_user_input=moa_user_input, exonic_pathogenic_user_input=exonic_pathogenic_user_input,
            exon_skipping_overrides=exon_skipping_overrides, knockdown_overrides=knockdown_overrides,
            wt_upregulation_overrides=wt_upregulation_overrides, gene_characteristics_overrides=gene_characteristics_overrides,
            splice_switching_overrides=splice_switching_overrides, client=client
        )

        _t0 = time.perf_counter()
        primary_result = _assess_for_transcript(target_consequence, **common_args)
        _log_timing("Primary transcript assessment (total)", _t0)

        # Offer a switch to the MANE Select/canonical transcript only when it's actually
        # a different transcript of the same gene (avoids toggling into an unrelated
        # overlapping gene's transcript in rare ambiguous-locus cases).
        alternate_consequence = None
        if (
            input_consequence and default_consequence
            and _base_transcript_id(default_consequence.get('transcript_id')) != _base_transcript_id(input_consequence.get('transcript_id'))
            and default_consequence.get('gene_symbol') == target_consequence.get('gene_symbol')
        ):
            alternate_consequence = default_consequence

        if alternate_consequence:
            _t0 = time.perf_counter()
            primary_result["alternative_result"] = _assess_for_transcript(alternate_consequence, **common_args)
            _log_timing("Alternate transcript assessment (total)", _t0)

        transcript_info = _build_transcript_info(requested_transcript_id, input_consequence, default_consequence, bool(alternate_consequence))
        if transcript_info:
            primary_result["transcript_info"] = transcript_info

        # Without the N1C tables, existing-ASO matches silently can't be found - say so
        if _n1c_unavailable():
            for result in (primary_result, primary_result.get("alternative_result")):
                if result:
                    _add_summary_note(result, N1C_UNAVAILABLE_NOTE)

        return primary_result

    except Exception:
        traceback.print_exc()
        return {"classification": "Error", "reason": "An unexpected server error occurred. Please try again later."}

# --- Main Flask Routes ---
# No automatic static folder: DOCS_DIR also holds docs/data, which must not be downloadable
app = Flask(__name__, static_folder=None)
# Requests carry a short query plus a few override flags; anything larger is rejected with 413
app.config['MAX_CONTENT_LENGTH'] = 64 * 1024

load_databases()

def _send_asset(directory, filename):
    """Serves an image file directly inside directory; nothing in subdirectories or of other types."""
    if '/' in filename or '\\' in filename or os.path.splitext(filename)[1].lower() not in STATIC_ASSET_EXTENSIONS:
        abort(404)
    return send_from_directory(directory, filename)

@app.route('/static/<path:filename>')
def static_assets(filename):
    """Serves the logo and figures in the DOCS_DIR root."""
    return _send_asset(DOCS_DIR, filename)

@app.route('/guidelines/images/<path:filename>')
def guideline_images(filename):
    """Serves the figures of the guidelines page."""
    return _send_asset(GUIDELINES_IMAGES_DIR, filename)

@app.route('/')
def index(): return send_from_directory(DOCS_DIR,'index.html')
@app.route('/usage')
def usage(): return send_from_directory(DOCS_DIR,'usage.html')
@app.route('/guidelines')
def guidelines(): return send_from_directory(DOCS_DIR,'guidelines.html')
@app.route('/about')
def about(): return send_from_directory(DOCS_DIR,'about.html')
@app.route('/cite')
def cite(): return send_from_directory(DOCS_DIR,'cite.html')

@app.route('/api_docs')
def api_docs():
    """Serves the API documentation page."""
    return send_from_directory(DOCS_DIR,'api_docs.html')

@app.route('/api/v1/assess', methods=['POST'])
def api_assess():
    """
    Handles a single variant assessment via a POST request for programmatic access.
    Returns the full assessment data as JSON.
    """
    data = request.get_json(silent=True)
    if not isinstance(data, dict) or not isinstance(data.get('query'), str):
        return jsonify({"error": "Request body must be JSON with a string 'query' key."}), 400

    query = data['query']
    if len(query) > MAX_QUERY_LENGTH:
        return jsonify({"error": f"Query exceeds the maximum length of {MAX_QUERY_LENGTH} characters."}), 400
    if not parse_hgvs_query(query)[0]:
        return jsonify({"error": "Invalid input format. Please use a recognized HGVS format (e.g., 'GENE c.123A>G')."}), 400

    client = EnsemblClient()
    _t0 = time.perf_counter()
    result = process_single_variant(query, client)
    _log_timing("Total assessment", _t0)

    # Provide specific HTTP status codes based on the outcome
    classification = result.get("classification")
    if classification == "Error":
        return jsonify({"error": result.get("reason", "An internal server error occurred.")}), 500
    if classification == "Unable to Assess":
        return jsonify({"error": result.get("reason", "Could not assess the provided variant.")}), 404
    
    return jsonify(result)

@app.route('/assess', methods=['POST'])
def assess():
    """Handles a single variant assessment request from the frontend."""
    data = request.get_json(silent=True)
    if not isinstance(data, dict) or not isinstance(data.get('query'), str):
        return jsonify({"classification": "Error", "reason": "No query provided."}), 400
    if len(data['query']) > MAX_QUERY_LENGTH:
        return jsonify({"classification": "Error", "reason": f"Query exceeds the maximum length of {MAX_QUERY_LENGTH} characters."}), 400
    # Overrides are looked up with .get(); anything that isn't a JSON object is ignored
    as_dict = lambda value: value if isinstance(value, dict) else {}
    query = data['query']
    splice_input = data.get('splice_user_input')
    exonic_pathogenic_user_input = data.get('exonic_pathogenic_user_input')
    moa_input = data.get('moa_user_input')
    es_overrides = as_dict(data.get('exon_skipping_overrides'))
    kd_overrides = as_dict(data.get('knockdown_overrides'))
    wt_overrides = as_dict(data.get('wt_upregulation_overrides'))
    gene_char_overrides = as_dict(data.get('gene_characteristics_overrides')) or None
    ss_overrides = as_dict(data.get('splice_switching_overrides'))
    client = EnsemblClient()
    _t0 = time.perf_counter()
    result = process_single_variant(query, client, splice_user_input=splice_input, moa_user_input=moa_input, exonic_pathogenic_user_input=exonic_pathogenic_user_input, exon_skipping_overrides=es_overrides, knockdown_overrides=kd_overrides, wt_upregulation_overrides=wt_overrides, gene_characteristics_overrides=gene_char_overrides, splice_switching_overrides=ss_overrides)
    _log_timing("Total assessment", _t0)

    return jsonify(result)

if __name__ == '__main__':
    app.run()
