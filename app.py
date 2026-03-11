import os
from flask import Flask, request, jsonify, render_template, send_file
import time
import re
import requests
from functools import lru_cache
from decimal import Decimal, getcontext
from typing import Dict, Any, Optional, Tuple, List
from Bio.Seq import Seq 
import pandas as pd
import io

# ===== CONFIGURATION =====
ENSEMBL_REST = "https://rest.ensembl.org"
HEADERS = {"Content-Type": "application/json", "Accept": "application/json"}
N1C_API_URL = "https://gene-registry.onrender.com/api/data?table=N1C_projects" 
N1C_API_ASSESSED_URL = "https://gene-registry.onrender.com/api/data?table=assessed_variants"
# Global DataFrames to be loaded at startup
clingen_df: Optional[pd.DataFrame] = None
goflof_df: Optional[pd.DataFrame] = None
splicevar_df: Optional[pd.DataFrame] = None
n1c_variants_df: Optional[pd.DataFrame] = None 
n1c_assessed_df: Optional[pd.DataFrame] = None 
sscvdb_df: Optional[pd.DataFrame] = None 
rcnv_df: Optional[pd.DataFrame] = None

# Supplementary N1C table with gene-level features (uORF, NAT, PE)
n1c_supp_df: Optional[pd.DataFrame] = None 

# --- Data Loading ---
def load_databases():
    """
    Loads all necessary data files and fetches N1C registry data
    into global pandas DataFrames.
    """
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    DATA_DIR = os.path.join(BASE_DIR, 'data')

    clingen_path = os.path.join(DATA_DIR, 'Clingen-Curation-Activity-Summary-Report-2025-10-15.csv')
    goflof_path = os.path.join(DATA_DIR, 'goflof_HGMD2019_v032021_allfeat.csv')
    splicevar_path = os.path.join(DATA_DIR, 'splicevardb.xlsx')
    sscvdb_path = os.path.join(DATA_DIR, 'SSCVDB.xlsx')
    n1c_supp_path = os.path.join(DATA_DIR, 'N1C_Variant_Supp_Table.xlsx')

    global clingen_df, goflof_df, splicevar_df, n1c_variants_df, n1c_assessed_df, sscvdb_df, n1c_supp_df
    try:
        clingen_df = pd.read_csv(clingen_path).set_index('gene_symbol')
        goflof_df = pd.read_csv(goflof_path).set_index('GENE')
        
        # Load SpliceVarDB from Excel
        splicevar_df = pd.read_excel(splicevar_path)
        # Sanitize SpliceVarDB data (critical for lookups)
        splicevar_df.columns = splicevar_df.columns.str.strip()
        if 'hgvs' in splicevar_df.columns and 'gene' in splicevar_df.columns:
            splicevar_df['hgvs'] = splicevar_df['hgvs'].astype(str).str.strip()
            splicevar_df['gene'] = splicevar_df['gene'].astype(str).str.strip()

        # Load SSCVDB from Excel
        try:
            sscvdb_df = pd.read_excel(sscvdb_path)
            sscvdb_df.columns = sscvdb_df.columns.str.strip()
            if 'Variant ID' in sscvdb_df.columns:
                sscvdb_df['Variant ID'] = sscvdb_df['Variant ID'].astype(str).str.strip()
        except Exception as e:
            print(f"Warning: Could not load SSCVDB.xlsx: {e}")

        # Load N1C supplementary table (uORF / NAT / PE per gene)
        try:
            n1c_supp_df = pd.read_excel(n1c_supp_path)
            n1c_supp_df.columns = n1c_supp_df.columns.str.strip()
            for col in ['Gene', 'uORF', 'NAT', 'PE']:
                if col in n1c_supp_df.columns:
                    n1c_supp_df[col] = n1c_supp_df[col].astype(str).str.strip()

        except Exception as e:
            print(f"Warning: Could not load N1C_Variant_Supp_Table.xlsx: {e}")

        # Fetch and load N1C variants data
        response = requests.get(N1C_API_URL, timeout=30)
        response.raise_for_status() # Will raise an error if the request fails
        n1c_data = response.json()
        n1c_variants_df = pd.DataFrame(n1c_data)
        
        # Sanitize the N1C columns we will search on
        if 'Gene' in n1c_variants_df.columns:
            n1c_variants_df['Gene'] = n1c_variants_df['Gene'].astype(str).str.strip()
        if 'Coding DNA change (c.)' in n1c_variants_df.columns:
            n1c_variants_df['Coding DNA change (c.)'] = n1c_variants_df['Coding DNA change (c.)'].astype(str).str.strip()

        # Fetch and load N1C assessed variants (curated) data
        response2 = requests.get(N1C_API_ASSESSED_URL, timeout=30)
        response2.raise_for_status()
        n1c_assessed_data = response2.json()
        n1c_assessed_df = pd.DataFrame(n1c_assessed_data)
        # Sanitize columns commonly used for matching
        if 'Gene' in n1c_assessed_df.columns:
            n1c_assessed_df['Gene'] = n1c_assessed_df['Gene'].astype(str).str.strip()
        # Normalize any potential c. notation columns to a unified helper accessor later
        for col in list(n1c_assessed_df.columns):
            if isinstance(col, str):
                n1c_assessed_df[col] = n1c_assessed_df[col].astype(str).str.strip()

        rcnv_path = os.path.join(DATA_DIR, 'rCNV.gene_scores.tsv')
        global rcnv_df
        if os.path.exists(rcnv_path):
            rcnv_df = pd.read_csv(rcnv_path, sep='\t')
            # Standardize gene column for matching
            if 'gene' in rcnv_df.columns:
                rcnv_df['gene'] = rcnv_df['gene'].astype(str).str.strip().str.upper()
        else:
            print(f"Warning: rCNV file not found at {rcnv_path}")

    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"An error occurred during database loading: {e}")
        exit(1)
# --- Add this method to EnsemblClient ---
def get_protein_metadata(self, protein_id):
    """Fetches metadata for a protein, including length."""
    if not protein_id: return None
    return self._get(f"/lookup/id/{protein_id}")

def _format_sscvdb_variant_id_from_vep(vep_entry: Dict[str, Any]) -> Optional[str]:
    """Formats a VEP entry to SSCVDB Variant ID style: chr<chrom>-<pos>-<ref>-<alt>.
    Returns None if required fields are missing or allele string is ambiguous."""
    try:
        chrom = str(vep_entry.get('seq_region_name', '')).strip()
        pos = vep_entry.get('start')
        allele = vep_entry.get('allele_string')
        if not chrom or not pos or not allele:
            return None
        if '/' not in allele:
            return None
        ref, alt = allele.split('/', 1)
        chrom_prefixed = f"chr{chrom}" if not chrom.lower().startswith('chr') else chrom
        return f"{chrom_prefixed}-{pos}-{ref}-{alt}"
    except Exception:
        return None

class EnsemblClient:
    def __init__(self, base_url=ENSEMBL_REST, headers=HEADERS, delay=0.1):
        self.base_url = base_url.rstrip('/')
        self.session = requests.Session()
        self.session.headers.update(headers)
        self.delay = delay

    def _get(self, path, params=None, max_retries=5):
        url = f"{self.base_url}{path}"
        backoff = 1.0
        for attempt in range(max_retries):
            time.sleep(self.delay)
            try:
                resp = self.session.get(url, params=params, timeout=30)
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

    def lookup_id_expand(self, identifier): return self._get(f"/lookup/id/{identifier}", params={'expand': '1'})
    def vep_hgvs(self, hgvs_string): return self._get(f"/vep/human/hgvs/{hgvs_string.strip()}/?hgvs=1", params={'variant_class': 1})
    def get_cds_sequence(self, transcript_id):
        data = self._get(f"/sequence/id/{transcript_id}", params={"type": "cds"})
        return data.get("seq") if isinstance(data, dict) else None
    def get_domains(self, protein_id):
        all_features = self._get(f"/overlap/translation/{protein_id}", params={"feature": "protein_feature"})
        if not all_features or not isinstance(all_features, list): return []
        domain_sources = {src.lower() for src in (
            'CDD','Pfam','SMART','PROSITE profiles','DisProt','PROSITE patterns','PRINTS','TIGRFAM','ProDom', 'mobidb-lite', 'RepeatsDB'
        )}
        repeat_disorder_sources = {
            'mobidb-lite'
        }
        allowed_sources = domain_sources | repeat_disorder_sources
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
            if src == 'tigrfam' and identifier:
                return f"https://www.ebi.ac.uk/interpro/entry/tigrfams/{ident_upper}/"
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
    def overlap_region_variation(self, chrom, start, end):
        data = self._get(f"/overlap/region/human/{chrom}:{start}-{end}", params={'feature': 'variation'})
        return data if isinstance(data, list) else []
    def get_overlapping_genes(self, gene_id):
        """Fetches all genes that overlap with a given Ensembl Gene ID."""
        data = self._get(f"/overlap/id/{gene_id}", params={"feature": "gene"})
        return data if isinstance(data, list) else []
    def lookup_symbol(self, symbol):
        """Fetches gene data for a given symbol."""
        data = self._get(f"/lookup/symbol/human/{symbol}", params={'expand': '0'})
        return data if isinstance(data, dict) else None

    def lookup_symbol_expand(self, symbol):
        """Fetches expanded gene data for a given symbol."""
        data = self._get(f"/lookup/symbol/human/{symbol}", params={'expand': '1'})
        return data if isinstance(data, dict) else None
    
# --- Helper & Parsing Functions ---

def _evaluate_splice_variant_position(
    variant_hgvs: str, 
    vep_data: Dict[str, Any], 
    details: Dict[str, Any],
    exonic_pathogenic_user_input: Optional[str] = None,
    moa_user_input: Optional[str] = None
) -> Optional[Dict[str, Any]]:
    """
    Assesses a validated splice-altering variant based on its genomic position.
    This logic is shared by database-found variants and user-validated variants.
    """
    client = EnsemblClient()
    # Find the most relevant consequence to get gene and transcript
    all_consequences = vep_data.get('transcript_consequences', [])
    target_consequence = choose_best_consequence(all_consequences)
    if not target_consequence:
        return {"classification": "Unable to Assess", "reason": "Could not determine a target transcript from VEP response."}

    gene_symbol = target_consequence.get('gene_symbol')
    transcript_id = target_consequence.get('transcript_id')
    if not gene_symbol or not transcript_id:
        return {"classification": "Unable to Assess", "reason": "Could not determine gene symbol or transcript ID."}

    core_hgvs_match = re.search(r'(c\..*)', variant_hgvs, re.IGNORECASE)
    if not core_hgvs_match: 
        return None
    core_canonical_hgvs = core_hgvs_match.group(1).lower()

    result = {"details": details}
    consequence_terms = set(target_consequence.get('consequence_terms', []))
    
    is_intronic_by_consequence = 'intron_variant' in consequence_terms or 'splice_acceptor_variant' in consequence_terms or 'splice_donor_variant' in consequence_terms
    is_intronic_by_notation = '+' in core_canonical_hgvs or '-' in core_canonical_hgvs
    is_intronic = is_intronic_by_consequence or is_intronic_by_notation

    if is_intronic:
        dist_match = re.search(r'[+-](\d+)', core_canonical_hgvs)
        if not dist_match: 
            return None 
        dist = int(dist_match.group(1))
        result['details']['Distance to Canonical Splice Site'] = f"{dist} bp"
        
        if '+' in core_canonical_hgvs: # Downstream (e.g., c.123+1G>A)
            if dist <= 5:
                result.update({"classification": "Not Eligible", "reason": "Variant is a splice-altering variant located too close to the canonical splice site (<=+5bp). If there is functional evidence variant results in leaky exon skipping, consider exon inclusion."})
            elif 6 <= dist <= 50:
                result.update({"classification": "Unlikely Eligible", "reason": "Variant is a splice-altering variant located near the canonical splice site (+6-+50bp). AVEC only assesses dependent on location. If you have experimental validation that the branch point is not weakened, this should be considered 'likely eligible"})
            else:
                result.update({"classification": "Likely Eligible", "reason": "Variant is a splice-altering variant in a favorable deep-intronic position (>+50bp)."})
        
        elif '-' in core_canonical_hgvs: # Upstream (e.g., c.124-2A>G)
            if dist <= 5:
                result.update({"classification": "Not Eligible", "reason": "Variant is a splice-altering variant located too close to the canonical splice site (>=-5bp). If there is functional evidence variant results in leaky exon skipping, consider exon inclusion."})
            elif 6 <= dist <= 100:
                result.update({"classification": "Unlikely Eligible", "reason": "Variant is a splice-altering variant located near the canonical splice site (-6-(-100b)p). AVEC only assesses dependent on location. If you have experimental validation that the branch point is not weakened, this should be considered 'likely eligible'."})
            else:
                result.update({"classification": "Likely Eligible", "reason": "Variant is a splice-altering variant in a favorable deep-intronic position (<-100bp)."})
        
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
        result.update({"classification": "Unable to assess", "reason": "Something went wrong."})

    return result

def parse_hgvs_query(query: str) -> Tuple[Optional[str], Optional[str]]:
    """
    Parses a query string into a VEP-compatible HGVS string and an optional gene symbol.
    """
    query = query.strip()
    
    # Pattern 1: Handles formats like "NM_12345.6:c.123A>G" or "GENE:c.123A>G"
    match_colon = re.search(r'([^:]+):([cgnmp]\..*)', query, re.IGNORECASE)
    if match_colon:
        identifier = match_colon.group(1).strip()
        variant = match_colon.group(2).strip()
        hgvs_string = f"{identifier}:{variant}"
        # If the identifier is a transcript, we don't have a gene symbol from the query
        if identifier.startswith("NM_") or identifier.startswith("ENST"):
            return hgvs_string, None
        else:
            return hgvs_string, identifier

    # Pattern 2: Handles formats like "GENE c.123A>G"
    match_space = re.search(r'([A-Z0-9\-_]+)\s+([cgnmp]\..*)', query, re.IGNORECASE)
    if match_space:
        gene = match_space.group(1).strip()
        variant = match_space.group(2).strip()
        return f"{gene}:{variant}", gene
        
    return None, None

def classify_variant_clinsig(clinsig_field):
    if clinsig_field is None: return 'other'
    vals = [v.lower() for v in (clinsig_field if isinstance(clinsig_field, list) else [clinsig_field]) if isinstance(v, str)]
    if any('pathogenic' in v for v in vals) and not any('likely' in v for v in vals): return 'pathogenic'
    if any('likely pathogenic' in v for v in vals): return 'pathogenic'
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

def choose_best_consequence(consequences: List[Dict[str, Any]], canonical_id: Optional[str] = None, gene_symbol_from_query: Optional[str] = None) -> Optional[Dict[str, Any]]:
    """
    Selects the most relevant transcript consequence from a VEP response.
    """
    if not consequences: return None

    if gene_symbol_from_query:
        targeted_consequences = [
            c for c in consequences if c.get('gene_symbol', '').upper() == gene_symbol_from_query.upper()
        ]
        if targeted_consequences:
            consequences = targeted_consequences

    mane_select = [c for c in consequences if c.get('mane_select')]
    if mane_select: return mane_select[0]
    
    if canonical_id:
        canonical_id_base = canonical_id.split('.')[0]
        for c in consequences:
            if c.get('transcript_id', '').startswith(canonical_id_base): return c
            
    coding_consequences = sorted(
        [c for c in consequences if c.get('biotype') == 'protein_coding' and c.get('cds_end')],
        key=lambda c: c['cds_end'] - c['cds_start'] if c.get('cds_start') else -1, reverse=True
    )
    
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

def check_n1c_registry(gene_symbol: str, original_query: str, formatted_hgvs: str) -> Optional[Dict[str, Any]]:
    """
    Searches the pre-loaded N1C registry DataFrame for a matching variant.
    If a fully developed therapy exists, it returns an 'Eligible' assessment.
    If it's 'under development', it returns a special status to allow further assessment.
    """
    if n1c_variants_df is None or n1c_variants_df.empty or not gene_symbol:
        return None

    core_hgvs_match = re.search(r'(c\..*)', formatted_hgvs, re.IGNORECASE)
    if not core_hgvs_match:
        return None
    core_hgvs = core_hgvs_match.group(1)

    gene_matches = n1c_variants_df[n1c_variants_df['Gene'].str.upper() == gene_symbol.upper()]
    if gene_matches.empty:
        return None
        
    for index, row in gene_matches.iterrows():
        registry_c_dot = row.get('Coding DNA change (c.)')
        if registry_c_dot and core_hgvs.lower() in registry_c_dot.lower():
            status = str(row.get('Status', 'N/A')).strip()
            modality = row.get('Therapeutic Modality', 'N/A')
            therapypublication = row.get('Therapy Publication', '')
            n1cid = row.get('ID')

            try:
                if pd.notna(n1cid):
                    n1cid = int(float(n1cid))
            except (ValueError, TypeError):
                pass
            
            link = f"https://generegistry.n1collaborative.org/entry.html?id={n1cid}" if n1cid else None

            is_under_development = 'under development' in status.lower()
            has_publication = therapypublication and pd.notna(therapypublication) and str(therapypublication).strip().lower() not in ['nan', 'na', 'n/a', '']

            if is_under_development and not has_publication:
                return {
                    "classification": "Under Development",
                    "reason": (
                        f"<p>An ASO for a variant in <strong>{gene_symbol}</strong> is listed as 'Under Development' in the N=1 Collaborative Projects Registry.</p>"
                        f"<p>Standard assessment will be performed below.</p>"
                        f"<ul>"
                        f"<li><strong>Matched Variant:</strong> {registry_c_dot}</li>"
                        f"<li><strong>Status:</strong> {status}</li>"
                        f"<li><strong>Therapeutic Modality:</strong> {modality}</li>"
                        f"</ul>"
                        f"<a href='{link}' target='_blank' rel='noopener noreferrer'>Click here to view the N=1 Collaborative registry page.</a>"
                    ),
                    "link": link,
                    "continue_assessment": True
                }
            
            return {
                "classification": "Eligible",
                "reason": (
                    f"<p>A direct match for a variant in <strong>{gene_symbol}</strong> was found in the N=1 Collaborative Projects Registry.</p>"
                    f"<ul>"
                    f"<li><strong>Matched Variant:</strong> {registry_c_dot}</li>"
                    f"<li><strong>Status:</strong> {status}</li>"
                    f"<li><strong>Therapeutic Modality:</strong> {modality}</li>"
                    f"</ul>"
                    f"<a href='{link}' target='_blank' rel='noopener noreferrer'>Click here to view the N=1 Collaborative registry page.</a>"
                ),
                "link": link,
                "continue_assessment": False
            }
    return None

def _get_c_notation_from_row(row: pd.Series) -> Optional[str]:
    """Best-effort extraction of c. notation from an assessed variants row."""
    candidate_cols = [
        'Variant (c.)', 'Coding DNA change (c.)', 'c_dot', 'HGVS', 'HGVS.c', 'Variant', 'Variant_c'
    ]
    for col in candidate_cols:
        if col in row:
            val = row.get(col)
            if isinstance(val, str) and 'c.' in val:
                return val
    # As a fallback, try to find any field containing 'c.' pattern
    for col in row.index:
        try:
            val = row[col]
            if isinstance(val, str) and 'c.' in val:
                return val
        except Exception:
            continue
    return None

def check_n1c_assessed_variants(gene_symbol: str, formatted_hgvs: str) -> Optional[Dict[str, Any]]:
    """Checks the N1C assessed variants dataset for a curated match and returns a curated assessment."""
    if n1c_assessed_df is None or n1c_assessed_df.empty or not gene_symbol or not formatted_hgvs:
        return None

    core_hgvs_match = re.search(r'(c\..*)', formatted_hgvs, re.IGNORECASE)
    if not core_hgvs_match:
        return None
    core_hgvs = core_hgvs_match.group(1)

    gene_matches = n1c_assessed_df[n1c_assessed_df['Gene'].str.upper() == gene_symbol.upper()] if 'Gene' in n1c_assessed_df.columns else n1c_assessed_df
    if gene_matches.empty:
        return None

    def _normalize_curated_classification(row: pd.Series) -> Tuple[str, Optional[str]]:
        """Attempts to return (classification, reason) based on row fields."""
        # Known target classes for coloring
        allowed = {"Eligible", "Likely Eligible", "Unlikely Eligible", "Not Eligible", "Unable to Assess"}
        # Candidate columns that may encode eligibility
        class_cols = [
            'Eligibility', 'ASO Eligibility', 'Classification', 'classification', 'eligibility', 'Status'
        ]
        val = None
        for col in class_cols:
            if col in row and isinstance(row[col], str) and row[col].strip():
                val = row[col].strip()
                break
        classification = None
        if val:
            v = val.lower()
            if 'not eligible' in v:
                classification = 'Not Eligible'
            elif 'unlikely' in v:
                classification = 'Unlikely Eligible'
            elif 'likely' in v:
                classification = 'Likely Eligible'
            elif 'eligible' in v:
                classification = 'Eligible'
            elif 'unable' in v:
                classification = 'Unable to Assess'
        if classification is None:
            classification = 'Eligible' if 'curated' in (row.get('Tags','') or '').lower() else 'Unable to Assess'
        # Reason/rationale
        reason_cols = ['Notes', 'Rationale', 'Reason', 'Comment']
        reason = None
        for col in reason_cols:
            if col in row and isinstance(row[col], str) and row[col].strip():
                reason = row[col].strip()
                break
        return classification, reason

    for _, row in gene_matches.iterrows():
        c_not = _get_c_notation_from_row(row)
        if c_not and core_hgvs.lower() in str(c_not).lower():
            # Build link using ID when present; prefer variant_entry.html
            link = None
            variant_id_val = row.get('ID')
            if variant_id_val not in (None, "", float('nan')):
                try:
                    if pd.notna(variant_id_val):
                        variant_id_val = int(float(variant_id_val))
                    link = f"https://generegistry.n1collaborative.org/variant_entry.html?id={variant_id_val}"
                except Exception:
                    link = None
            if not link and 'Link' in row:
                link = row['Link']

            classification, curated_reason = _normalize_curated_classification(row)
            reason = (
                f"<p>A curated assessment exists in the N=1 Collaborative Assessed Variants database for <strong>{gene_symbol}</strong>.</p>"
                f"<p><strong>Matched Variant:</strong> {c_not}</p>"
            )
            if curated_reason:
                reason += f"<p><strong>Notes:</strong> {curated_reason}</p>"
            if link:
                reason += f"<p><a href='{link}' target='_blank' rel='noopener noreferrer'>View N1C curated entry</a></p>"
            return {
                "classification": classification,
                "reason": reason,
                "link": link
            }
    return None

def _extract_exon_numbers_from_text(text: str) -> List[int]:
    """Extracts exon numbers (including ranges) from free text."""
    if not isinstance(text, str) or not text:
        return []
    exons: List[int] = []
    try:
        # Ranges like "exons 45-55" or "exon 2-3"
        for m in re.finditer(r"exons?\s*(\d{1,3})\s*[-–to]+\s*(\d{1,3})", text, flags=re.IGNORECASE):
            a, b = int(m.group(1)), int(m.group(2))
            if a <= b:
                exons.extend(list(range(a, b + 1)))
            else:
                exons.extend(list(range(b, a + 1)))
        # Singles like "exon 51" or "exon-51"
        for m in re.finditer(r"exons?[-\s]*(\d{1,3})", text, flags=re.IGNORECASE):
            try:
                val = int(m.group(1))
                exons.append(val)
            except Exception:
                pass
    except Exception:
        return []
    # Deduplicate while preserving order
    seen = set()
    out: List[int] = []
    for e in exons:
        if e not in seen:
            seen.add(e)
            out.append(e)
    return out

def _row_is_exon_skipping(row: pd.Series) -> bool:
    """Checks if a registry row describes an exon-skipping approach."""
    approach_cols = ['Approach', 'Therapeutic Modality', 'Therapeutic approach', 'Modality']
    for col in approach_cols:
        if col in row and isinstance(row.get(col), str):
            if re.search(r"exon\s*skip", row.get(col), re.IGNORECASE):
                return True
    # Fallback: search any string field for exon skip
    for col in row.index:
        val = row[col]
        if isinstance(val, str) and re.search(r"exon\s*skip", val, re.IGNORECASE):
            return True
    return False

def _row_is_knockdown(row: pd.Series) -> bool:
    """Checks if a registry row describes a knockdown approach."""
    approach_cols = ['Approach', 'Therapeutic Modality', 'Therapeutic approach', 'Modality']
    for col in approach_cols:
        if col in row and isinstance(row.get(col), str):
            if re.search(r"knock\s*down|silenc", row.get(col), re.IGNORECASE):
                return True
    for col in row.index:
        val = row[col]
        if isinstance(val, str) and re.search(r"knock\s*down|silenc", val, re.IGNORECASE):
            return True
    return False

def n1c_exon_skipping_exon_numbers_for_gene(gene_symbol: str) -> Tuple[set, List[str], Dict[int, List[str]]]:
    """
    Returns (set_of_exon_numbers, list_of_links, map_exon_to_links) for N1C registry rows that
    indicate exon skipping for the given gene. Filters to Approach: Exon Skipping when available.
    """
    exon_set: set = set()
    links: set = set()
    exon_link_map: Dict[int, set] = {}
    if n1c_variants_df is None or n1c_variants_df.empty or not gene_symbol:
        return exon_set, [], {}
    try:
        if 'Gene' in n1c_variants_df.columns:
            gene_matches = n1c_variants_df[n1c_variants_df['Gene'].astype(str).str.upper() == gene_symbol.upper()]
        else:
            gene_matches = n1c_variants_df
    except Exception:
        gene_matches = n1c_variants_df
    if gene_matches is None or len(gene_matches) == 0:
        return exon_set, [], {}

    candidate_cols = [
        'Target Exon', 'Targeted Exon', 'Exon', 'Exons', 'Exon(s)', 'Exon(s) Targeted',
        'Exon Target', 'Exon Targets', 'Intervention Description', 'Therapeutic Modality',
        'Approach', 'Description', 'Notes', 'Summary', 'Treatment Details'
    ]

    for _, row in gene_matches.iterrows():
        if not _row_is_exon_skipping(row):
            continue

        texts: List[str] = []
        for col in candidate_cols:
            if col in row and isinstance(row.get(col), str) and row.get(col).strip():
                texts.append(row.get(col))

        if not texts:
            try:
                texts = [val for val in row if isinstance(val, str)]
            except Exception:
                texts = []

        nums: List[int] = []
        for t in texts:
            nums.extend(_extract_exon_numbers_from_text(t))
        # Deduplicate in order
        seen = set()
        dedup_nums = []
        for n in nums:
            if n not in seen:
                dedup_nums.append(n)
                seen.add(n)

        if not dedup_nums:
            continue

        link = None
        try:
            nid = row.get('ID')
            if nid is not None and str(nid).strip() and str(nid).strip().lower() != 'nan':
                link = f"https://generegistry.n1collaborative.org/entry.html?id={str(nid).strip()}"
                try:
                    nid = int(float(nid))
                except (ValueError, TypeError):
                    pass
                link = f"https://generegistry.n1collaborative.org/entry.html?id={nid}"
        except Exception:
            link = None

        for n in dedup_nums:
            exon_set.add(n)
            if link:
                links.add(link)
                exon_link_map.setdefault(n, set()).add(link)

    # Normalize map values to lists
    normalized_map: Dict[int, List[str]] = {k: sorted(v) for k, v in exon_link_map.items()}
    return exon_set, sorted(links), normalized_map

def n1c_exon_skipping_variant_exon_map(gene_symbol: str, all_exons: List[Dict[str, Any]], client: EnsemblClient) -> Dict[int, List[str]]:
    """
    Maps N1C registry exon-skipping variants (using RefSeq Transcript + c. change) to exon numbers,
    using the same HGVS parsing as user input. Returns {exon_number: [links]}.
    """
    exon_map: Dict[int, List[str]] = {}
    if n1c_variants_df is None or n1c_variants_df.empty or not gene_symbol or not all_exons:
        return exon_map

    try:
        if 'Gene' in n1c_variants_df.columns:
            gene_matches = n1c_variants_df[n1c_variants_df['Gene'].astype(str).str.upper() == gene_symbol.upper()]
        else:
            gene_matches = n1c_variants_df
    except Exception:
        gene_matches = n1c_variants_df

    if gene_matches is None or len(gene_matches) == 0:
        return exon_map

    for _, row in gene_matches.iterrows():
        if not _row_is_exon_skipping(row):
            continue

        refseq_val = None
        cdot_val = None
        for cand in ['RefSeq Transcript', 'Transcript', 'RefSeq', 'Transcript ID']:
            if cand in row and isinstance(row.get(cand), str) and row.get(cand).strip():
                refseq_val = row.get(cand).strip()
                break
        for cand in ['Coding DNA change (c.)', 'Variant (c.)', 'c_dot', 'c.']:
            if cand in row and isinstance(row.get(cand), str) and row.get(cand).strip():
                cdot_val = row.get(cand).strip()
                break

        if not refseq_val or not cdot_val:
            continue

        hgvs_input = f"{refseq_val}:{cdot_val}"
        try:
            hgvs_query, _ = parse_hgvs_query(hgvs_input)
            if not hgvs_query:
                continue
            vep_data = client.vep_hgvs(hgvs_query)
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
            link = None
            try:
                nid = row.get('ID')
                if nid is not None and str(nid).strip() and str(nid).strip().lower() != 'nan':
                    link = f"https://generegistry.n1collaborative.org/entry.html?id={str(nid).strip()}"
                    try:
                        nid = int(float(nid))
                    except (ValueError, TypeError):
                        pass
                    link = f"https://generegistry.n1collaborative.org/entry.html?id={nid}"
            except Exception:
                link = None
            if link:
                exon_map.setdefault(exon_num, []).append(link)
            else:
                exon_map.setdefault(exon_num, [])  # ensure key exists even without link
        except Exception:
            continue

    # Deduplicate links for each exon
    for ex_num, links in list(exon_map.items()):
        dedup = []
        seen = set()
        for l in links:
            if l not in seen:
                seen.add(l)
                dedup.append(l)
        exon_map[ex_num] = dedup
    return exon_map

def n1c_gene_knockdown_entry(gene_symbol: str) -> Optional[Dict[str, Any]]:
    """Returns a knockdown registry match at the gene level (Approach: Knockdown)."""
    if n1c_variants_df is None or n1c_variants_df.empty or not gene_symbol:
        return None
    try:
        if 'Gene' in n1c_variants_df.columns:
            gene_matches = n1c_variants_df[n1c_variants_df['Gene'].astype(str).str.upper() == gene_symbol.upper()]
        else:
            gene_matches = n1c_variants_df
    except Exception:
        gene_matches = n1c_variants_df
    if gene_matches is None or len(gene_matches) == 0:
        return None
    for _, row in gene_matches.iterrows():
        if not _row_is_knockdown(row):
            continue
        link = None
        try:
            nid = row.get('ID')
            if nid is not None and str(nid).strip() and str(nid).strip().lower() != 'nan':
                link = f"https://generegistry.n1collaborative.org/entry.html?id={str(nid).strip()}"
                try:
                    nid = int(float(nid))
                except (ValueError, TypeError):
                    pass
                link = f"https://generegistry.n1collaborative.org/entry.html?id={nid}"
        except Exception:
            link = None
        modality = row.get('Therapeutic Modality') if isinstance(row.get('Therapeutic Modality'), str) else row.get('Approach')
        reason = (
            f"N1C registry lists knockdown project(s) for {gene_symbol}. "
            f"Because the variant is Autosomal Dominant with GoF mechanism, it is considered eligible for knockdown."
        )
        details = {}
        if modality:
            details["Registry Modality"] = str(modality)
        if link:
            details["N1C Registry Link"] = link
        return {
            "classification": "Eligible",
            "reason": reason,
            "link": link,
            "details": details
        }
    return None

def get_gene_characteristics(gene_symbol: str) -> Dict[str, Any]:
    """
    Retrieves MOI, Haploinsufficiency, Triplosensitivity, and MOA from loaded dataframes.
    Includes rCNV scores (pHaplo/pTriplo) from Collins et al. 2022.
    """
    characteristics = {
        "moi": [],
        "haploinsufficiency": {"text": "Unknown", "url": None},
        "triplosensitivity": {"text": "Unknown", "url": None}, # Added init
        "moa": [],
        "gene_url": None,
        "rcnv": {"pHaplo": "N/A", "pTriplo": "N/A", "url": "https://doi.org/10.1016/j.cell.2022.06.036"} # Init rCNV
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
        gene_data_row = clingen_df.loc[gene_symbol]
        if isinstance(gene_data_row, pd.DataFrame):
            gene_data_row = gene_data_row.iloc[0]

        gene_url_val = gene_data_row.get('gene_url')
        if pd.notna(gene_url_val):
            characteristics['gene_url'] = gene_url_val

        # 1. Mode of Inheritance (MOI)
        moi_val = gene_data_row.get('mode_of_inheritance')
        if pd.notna(moi_val):
            all_mois = [moi.strip() for moi in str(moi_val).split(',')]
            characteristics["moi"] = sorted(list(set(all_mois)))

        # 2. Haploinsufficiency with Link
        hap_assertion = gene_data_row.get('dosage_haploinsufficiency_assertion')
        hap_url = gene_data_row.get('dosage_report')
        if pd.notna(hap_assertion):
            hap_score_str = str(hap_assertion).strip()
            haplo_text = "Unknown"
            if hap_score_str.startswith('3 -'): haplo_text = "Sufficient evidence"
            elif hap_score_str.startswith('1 -'): haplo_text = "Little evidence"
            elif hap_score_str.startswith('30 -'): haplo_text = "Gene associated with autosomal recessive phenotype"
            elif hap_score_str.startswith('40 -'): haplo_text = "Dosage sensitivity unlikely"
            else: haplo_text = "No evidence"
            
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

    # --- GOF/LOF Lookup ---
    if goflof_df is not None and gene_symbol in goflof_df.index:
        gene_data = goflof_df.loc[[gene_symbol]]
        if not gene_data[gene_data['LABEL'].str.contains('GOF', na=False)].empty:
            characteristics["moa"].append("GoF")
        if not gene_data[gene_data['LABEL'].str.contains('LOF', na=False)].empty:
            characteristics["moa"].append("LoF")
        if characteristics["moa"]:
            characteristics["moa"] = sorted(list(set(characteristics["moa"])))

    return characteristics

def assess_knockdown(gene_characteristics: Dict[str, Any], overrides: Dict[str, bool] = {}) -> Dict[str, Any]:
    """Assesses eligibility for a knockdown strategy."""
    haplo_obj = gene_characteristics.get("haploinsufficiency", {"text": "Unknown"})
    haplo_status_text = haplo_obj.get("text", "Unknown")
    
    moi = gene_characteristics.get("moi", [])
    # --- Logic to track original vs overridden state ---
    check_name = "Gene is not haploinsufficient"
    original_passed = haplo_status_text in ["No evidence", "Little evidence", "Dosage sensitivity unlikely"]
    
    has_override = check_name in overrides
    final_passed = overrides.get(check_name, original_passed)

    checks = {
        check_name: {
            "passed": final_passed,
            "original_passed": original_passed,
            "overridden": has_override
        }
    }
    # Determine if haploinsufficiency status is unknown/missing for warning purposes
    haplo_unknown = (haplo_status_text is None) or (str(haplo_status_text).strip().lower() in ["", "unknown", "n/a"])
    
    # Classification logic now depends on the final 'is_not_haploinsufficient' value
    
    if haplo_unknown and not has_override:
        reason = "Gene haploinsufficiency status is not available in datasources. Therefore knockdown could lead to unintended consequences and AVEC cannot determine amenability to knockdown."
        reason = "ClinGen haploinsufficiency status is unknown, therefore AVEC is unable to assess the variant for knockdown. If there is evidence for haploinsufficiency status use the override to reassess"
        return {
            "classification": "Unable to Assess",
            "reason": reason,
            "checks": checks
        }
    elif not final_passed:
        reason = "Gene is associated with haploinsufficiency or pathogenic autosomal dominant loss-of-function variants."
        return {
            "classification": "Unlikely Eligible",
            "reason": reason,
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
    
def get_overlapping_genes(self, gene_id):
        """Fetches all genes that overlap with a given Ensembl Gene ID."""
        data = self._get(f"/overlap/id/{gene_id}", params={"feature": "gene"})
        return data if isinstance(data, list) else [] # This is a duplicate method

def assess_wt_upregulation(client: EnsemblClient, gene_id: str, gene_symbol: str) -> Dict[str, Any]:
    """
    Assesses for WT upregulation by checking for overlapping NATs and by
    searching for a conventional antisense gene name ([GENE]-AS1).
    Trusts the '-AS1' naming convention without a strict biotype check.
    """
    if not gene_id or not gene_symbol:
        return {"classification": "Unable to Assess", "reason": "Missing Gene ID or Symbol."}

    found_antisense_genes = {} 

    # Collect curated gene-level features from N1C supplementary table
    supp_details: Dict[str, str] = {}
    try:
        if 'n1c_supp_df' in globals() and n1c_supp_df is not None and 'Gene' in n1c_supp_df.columns:
            matches = n1c_supp_df[n1c_supp_df['Gene'].astype(str).str.strip() == str(gene_symbol).strip()]
            def _norm(val: Any) -> str:
                v = str(val).strip().upper()
                if v in ("Y", "YES"): return "Available"
                if v in ("N", "NO"): return "Not available"
                return "Unknown"
            if not matches.empty:
                row = matches.iloc[0]
                supp_details["uORF"] = _norm(row.get('uORF', 'N/A'))
                supp_details["NAT (curated)"] = _norm(row.get('NAT', 'N/A'))
                supp_details["Poison exon (PE)"] = _norm(row.get('PE', 'N/A'))
            else:
                supp_details["uORF"] = "Unknown"
                supp_details["NAT (curated)"] = "Unknown"
                supp_details["Poison exon (PE)"] = "Unknown"
        else:
            supp_details["uORF"] = "Unknown"
            supp_details["NAT (curated)"] = "Unknown"
            supp_details["Poison exon (PE)"] = "Unknown"
    except Exception:
        # Fail-quiet: do not block assessment if supplemental table can't be parsed
        supp_details["uORF"] = "Unknown"
        supp_details["NAT (curated)"] = "Unknown"
        supp_details["Poison exon (PE)"] = "Unknown"

    try:
        # --- Method 1: Search by genomic coordinate overlap ---
        overlapping_genes = client.get_overlapping_genes(gene_id)
        for gene in overlapping_genes:
            if gene.get('biotype') == 'antisense' and gene.get('id') != gene_id:
                found_antisense_genes[gene['id']] = gene

        # --- Method 2: Search by conventional name ([GENE_SYMBOL]-AS1) ---
        antisense_symbol = f"{gene_symbol}-AS1"
        as_gene = client.lookup_symbol(antisense_symbol)
        
        if as_gene:
            found_antisense_genes[as_gene['id']] = as_gene

        # --- Evaluate curated evidence (uORF / NAT / PE) and Ensembl NAT search ---
        has_uorf = supp_details.get("uORF") == "Available"
        has_nat_cur = supp_details.get("NAT (curated)") == "Available"
        has_pe = supp_details.get("Poison exon (PE)") == "Available"
        curated_features = [name for ok, name in [
            (has_uorf, "uORF"),
            (has_nat_cur, "NAT"),
            (has_pe, "Poison exon")
        ] if ok]

        nat_list = []
        nat_names = ""
        nat_ids: List[str] = []
        if found_antisense_genes:
            nat_list = list(found_antisense_genes.values())
            nat_names = ", ".join([nat.get('external_name', nat['id']) for nat in nat_list])
            nat_ids = [nat['id'] for nat in nat_list]

        # Start details with curated statuses and cite the original guideline table
        details = dict(supp_details)
        details["Datasource: Table S2 of N1C VARIANT Guidelines"] = "https://www.sciencedirect.com/science/article/pii/S0002929725000643#app2"
        # Add Ensembl NAT links if any were found
        if nat_list:
            for nat in nat_list:
                nat_name = nat.get('external_name', nat['id'])
                ensembl_link = f"https://www.ensembl.org/Homo_sapiens/Gene/Summary?g={nat['id']}"
                details[nat_name] = ensembl_link

        if curated_features:
            reason = f"Curated evidence indicates available: {', '.join(curated_features)} for {gene_symbol}."
            if nat_names:
                reason += f" Additionally, Ensembl found antisense transcript(s): {nat_names}."
            return {
                "classification": "Potential possibilities identified",
                "reason": reason,
                "details": details,
                "antisense_gene_ids": nat_ids,
                "checks": {"Overlapping antisense transcript found": bool(nat_list)}
            }
        else:
            reason = "No curated evidence for uORF, NAT or poison exon availability."
            if nat_names:
                reason += f" However, Ensembl reports antisense transcript(s): {nat_names}."
            return {
                "classification": "No potential possibilities identified",
                "reason": reason,
                "details": details,
                "antisense_gene_ids": nat_ids,
                "checks": {"Overlapping antisense transcript found": bool(nat_list)}
            }

    except Exception as e:
        return {
            "classification": "Unable to Assess",
            "reason": f"An error occurred while searching for antisense transcripts: {e}",
            "checks": {}
        }

def assess_splice_switching(variant_hgvs: str, vep_data: Dict[str, Any], gene_symbol: str, exonic_pathogenic_user_input: Optional[str] = None) -> Optional[Dict[str, Any]]:
    """
    Assesses a variant for splice-switching potential, adding method and DOI link.
    If not found in the DB, it returns a prompt for user validation.
    """
    if splicevar_df is None or not variant_hgvs or not gene_symbol:
        return None 

    clean_gene = gene_symbol.strip().upper()
    full_canonical_hgvs = variant_hgvs.strip()
    core_hgvs_match = re.search(r'(c\..*)', full_canonical_hgvs, re.IGNORECASE)
    if not core_hgvs_match: return None
    core_canonical_hgvs = core_hgvs_match.group(1).lower()

    gene_rows = splicevar_df[splicevar_df['gene'].str.strip().str.contains(clean_gene, case=False, na=False)]
    if gene_rows.empty:
        # Not found in SpliceVarDB for this gene -> check SSCVDB fallback
        if sscvdb_df is not None and not sscvdb_df.empty and vep_data:
            variant_key = _format_sscvdb_variant_id_from_vep(vep_data)
            if variant_key and 'Variant ID' in sscvdb_df.columns:
                if not sscvdb_df[sscvdb_df['Variant ID'].str.strip().str.lower() == variant_key.strip().lower()].empty:
                    details = {
                        "Source Database": "SSCVDB",
                        "Evidence": "Splice-altering reported in SSCVDB",
                        "SSCVDB Gene Page": f"https://sscvdb.io/gene/{gene_symbol}"
                    }
                    return _evaluate_splice_variant_position(variant_hgvs, vep_data, details, exonic_pathogenic_user_input=exonic_pathogenic_user_input)
        # --- Return prompt if not in either database ---
        return {
            "classification": "Not in Database",
            "reason": "This variant was not found in SpliceVarDB or SSCVDB. If there is experimental validation (by qPCR or RNA-seq) that this variant is splice-altering, please confirm below.",
            "user_validation_prompt": True
        }

    for index, row in gene_rows.iterrows():
        db_hgvs = str(row.get('hgvs', '')).strip()
        db_hgvs_lower = db_hgvs.lower()
        
        if db_hgvs_lower.endswith(core_canonical_hgvs):
            splice_info = row
            method = str(splice_info.get('method', 'N/A')).strip()
            classification = str(splice_info.get('classification', 'N/A')).strip().lower()
            has_valid_method = method.lower() in ['rna-seq', 'rt-pcr', 'minigene']
            is_splice_altering = classification == 'splice-altering'

            if not (has_valid_method and is_splice_altering):
                continue

            details = {"Confirmation Method": method}
            doi = str(splice_info.get('doi', '')).strip()
            if doi and doi.lower() not in ['na', 'n/a', '']:
                details["Publication"] = f"https://doi.org/{doi}"

            return _evaluate_splice_variant_position(variant_hgvs, vep_data, details, exonic_pathogenic_user_input=exonic_pathogenic_user_input)

    # --- If no exact HGVS match was found in SpliceVarDB, try SSCVDB before prompting ---
    if sscvdb_df is not None and not sscvdb_df.empty and vep_data:
        variant_key = _format_sscvdb_variant_id_from_vep(vep_data)
        if variant_key and 'Variant ID' in sscvdb_df.columns:
            if not sscvdb_df[sscvdb_df['Variant ID'].str.strip().str.lower() == variant_key.strip().lower()].empty:
                details = {
                    "Source Database": "SSCVDB",
                    "Evidence": "Splice-altering reported in SSCVDB",
                    "SSCVDB Gene Page": f"https://sscvdb.io/gene/{gene_symbol}"
                }
                return _evaluate_splice_variant_position(variant_hgvs, vep_data, details, exonic_pathogenic_user_input=exonic_pathogenic_user_input)

    # --- Return prompt object if still not found ---
    return {
        "classification": "Not in Database",
        "reason": "This variant was not found in SpliceVarDB. If there is experimental validation (by qPCR or RNA-seq) that this variant is splice-altering, please confirm below.",
        "user_validation_prompt": True
    }

def assess_consecutive_exons(client, transcript, all_exons, target_exon, vep_entry, refseq_id_for_viewer, overrides={}):
    coding_exons = sorted([e for e in all_exons if e['cds_length'] > 0], key=lambda x: x['coding_exon_number'])
    
    target_idx = next((i for i, ex in enumerate(coding_exons) if ex['exon_id'] == target_exon['exon_id']), -1)
    if target_idx == -1: return None

    # Helper to assess a specific pair of exons
    def evaluate_pair(exon_a, exon_b):
        gene_id = transcript.get('Parent')
        transcript_id = transcript.get('id')
        protein_id = transcript.get("Translation", {}).get("id")
        protein_meta = client._get(f"/lookup/id/{protein_id}") if protein_id else None
        
        total_cds_len = sum(e['cds_length'] for e in coding_exons)
        total_protein_len = protein_meta['length'] if (protein_meta and 'length' in protein_meta) else (total_cds_len / 3 if total_cds_len > 0 else 0)
        
        combined_cds_len = exon_a['cds_length'] + exon_b['cds_length']
        combined_aa_len = combined_cds_len / 3
        combined_len_nt = (exon_a['end'] - exon_a['start'] + 1) + (exon_b['end'] - exon_b['start'] + 1)
        
        # Stop codon check
        cds_seq = client.get_cds_sequence(transcript_id)
        orig_cond2_no_stop = False
        orig_cond10_no_new_codon = True
        if cds_seq:
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
        
        # Variants
        def get_counts(chrom, start, end):
            vars_in_reg = client.overlap_region_variation(chrom, start, end)
            c = {'missense': 0, 'inframe_del': 0, 'splice': 0, 'nonsense': 0, 'frameshift': 0, 'benign_splice':0}
            links = {"pathogenic": [], "benign": []}
            seen = {"pathogenic": set(), "benign": set()}
            for v in vars_in_reg:
                clclass = classify_variant_clinsig(v.get('clinical_significance'))
                conseq = (v.get("consequence_type") or "").lower()
                vid = str(v.get('id') or v.get('variation_name') or '').strip()
                vsrc = (v.get('source') or '').strip() or "Ensembl"
                vlink = _build_variant_link(v) if vid else None
                if clclass == "pathogenic":
                    if "missense" in conseq: c['missense'] += 1
                    elif "inframe_deletion" in conseq: c['inframe_del'] += 1
                    elif "splice_donor" in conseq or "splice_acceptor" in conseq: c['splice'] += 1
                    elif "stop_gained" in conseq: c['nonsense'] += 1
                    elif "frameshift" in conseq: c['frameshift'] += 1
                    if ("splice_donor" in conseq or "splice_acceptor" in conseq) and vid and vid not in seen["pathogenic"]:
                        links["pathogenic"].append({"id": vid, "source": vsrc, "url": vlink})
                        seen["pathogenic"].add(vid)
                elif clclass == 'benign' and ("splice_donor" in conseq or "splice_acceptor" in conseq):
                    c['benign_splice'] += 1
                    if vid and vid not in seen["benign"]:
                        links["benign"].append({"id": vid, "source": vsrc, "url": vlink})
                        seen["benign"].add(vid)
            return c, links

        c1, l1 = get_counts(exon_a['seq_region_name'], exon_a['start'], exon_a['end'])
        c2, l2 = get_counts(exon_b['seq_region_name'], exon_b['start'], exon_b['end'])
        counts = {k: c1.get(k, 0) + c2.get(k, 0) for k in c1}
        splice_variant_links = {"pathogenic": l1["pathogenic"] + l2["pathogenic"], "benign": l1["benign"] + l2["benign"]}
        
        orig_cond6_missense = counts['missense'] < 3 + counts['frameshift'] + counts['nonsense']
        orig_cond7_splice = counts['splice'] == 0
        orig_cond8_no_inframe_del = counts['inframe_del'] == 0
        
        # Apply Overrides
        cond2_no_stop = overrides.get('No New Stop Codon', orig_cond2_no_stop)
        cond3_not_terminal = overrides.get('Not First/Last Exon', orig_cond3_not_terminal)
        cond4_small = overrides.get('Is <10% of Protein', orig_cond4_small)
        cond5_no_domain = overrides.get('No Domain Overlap', orig_cond5_no_domain)
        cond6_missense = overrides.get('Low Missense Count', orig_cond6_missense)
        cond7_splice = overrides.get('No Pathogenic Splice Variants', orig_cond7_splice)
        cond8_no_inframe_del = overrides.get('No Pathogenic In-Frame Deletions', orig_cond8_no_inframe_del)
        cond10_no_new_codon = overrides.get('No New Codon Formation', orig_cond10_no_new_codon)

        classification, reason = "Undetermined", ""
        if not cond3_not_terminal: classification, reason = "Not Eligible", "One of the exons is the first or last coding exon."
        elif not cond2_no_stop: classification, reason = "Not Eligible", "Skipping these exons is predicted to create a premature stop codon."
        elif not cond7_splice: classification, reason = "Not Eligible", f"Exons contain {counts['splice']} pathogenic splice variant(s)."
        elif not cond4_small and domain_count > 1: classification, reason = "Not Eligible", f"Exons overlap with {domain_count} protein domains and constitute >10% of the coding region."
        elif not cond8_no_inframe_del: classification, reason = "Not Eligible", f"Exons contain {counts['inframe_del']} pathogenic in-frame deletion(s)."
        elif not cond10_no_new_codon: classification, reason = "Unlikely Eligible", "Skipping these exons creates a novel codon at the junction that codes for a different amino acid."
        elif not cond6_missense: classification, reason = "Unlikely Eligible", f"Exons are a mutational hotspot with {counts['missense']} pathogenic missense variants."
        elif not cond4_small: classification, reason = "Unlikely Eligible", "Exons constitute >=10% of the protein."
        elif not cond5_no_domain: classification, reason = "Unlikely Eligible", f"Exons overlap with {domain_count} protein domain(s)."
        else: classification, reason = "Likely Eligible", "Exons meet the primary criteria for multi-exon skipping."
        
        note = f"Assessment based on skipping two consecutive exons: Exon {exon_a['total_exon_number']} and Exon {exon_b['total_exon_number']}."
        if combined_len_nt > 300:
            note += " Combined exon length exceeds 300 nucleotides. Proceed with caution."
        note += " Please note that some exceptions for exon skipping assessment apply. Please refer to the methods section for more details."
            
        visualization_data = None
        try:
            v_chrom, v_start, v_end = vep_entry.get('seq_region_name'), vep_entry.get('start'), vep_entry.get('end')
            padding = 1000
            visualization_data = {
                "locus": f"{v_chrom}:{max(1, v_start - padding)}-{v_end + padding}",
                "variantTrack": {"name": "Variant", "features": [{"chr": v_chrom, "start": v_start - 1, "end": v_end, "name": vep_entry.get('id', 'Variant')}]},
                "domainTrack": None
            }
        except: pass

        return {
            "classification": classification, "reason": reason, "note": note, "exon_length": combined_len_nt,
            "frac_cds": f"{((Decimal(str(combined_aa_len)) / Decimal(str(total_protein_len))) * 100):.2f}%" if total_protein_len > 0 else "N/A",
            "pathogenic_variant_counts": counts, "domain_count": domain_count, "domain_names": overlapping_domain_names,
            "domain_details": overlapping_domain_details, "splice_variant_links": splice_variant_links,
            "coding_exon_number": f"{exon_a['coding_exon_number']}+{exon_b['coding_exon_number']}",
            "total_exon_number": f"{exon_a['total_exon_number']}+{exon_b['total_exon_number']}",
            "gene_id": gene_id, "transcript_id": transcript_id,
            "clinvar_url": f"https://www.ncbi.nlm.nih.gov/clinvar/?term=GRCh38%3A{exon_a['seq_region_name']}%3A{min(exon_a['start'], exon_b['start'])}-{max(exon_a['end'], exon_b['end'])}",
            "checks": {
                "Is In-Frame": {"passed": True, "original_passed": True, "overridden": 'Is In-Frame' in overrides},
                "No New Stop Codon": {"passed": cond2_no_stop, "original_passed": orig_cond2_no_stop, "overridden": 'No New Stop Codon' in overrides},
                "Not First/Last Exon": {"passed": cond3_not_terminal, "original_passed": orig_cond3_not_terminal, "overridden": 'Not First/Last Exon' in overrides},
                "No Pathogenic Splice Variants": {"passed": cond7_splice, "original_passed": orig_cond7_splice, "overridden": 'No Pathogenic Splice Variants' in overrides},
                "No Pathogenic In-Frame Deletions": {"passed": cond8_no_inframe_del, "original_passed": orig_cond8_no_inframe_del, "overridden": 'No Pathogenic In-Frame Deletions' in overrides},
                "No Domain Overlap": {"passed": cond5_no_domain, "original_passed": orig_cond5_no_domain, "overridden": 'No Domain Overlap' in overrides},
                "No New Codon Formation": {"passed": cond10_no_new_codon, "original_passed": orig_cond10_no_new_codon, "overridden": 'No New Codon Formation' in overrides},
                "Low Missense Count": {"passed": cond6_missense, "original_passed": orig_cond6_missense, "overridden": 'Low Missense Count' in overrides},
                "Is <10% of Protein": {"passed": cond4_small, "original_passed": orig_cond4_small, "overridden": 'Is <10% of Protein' in overrides}
            },
            "visualization": visualization_data
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
    
    results = []
    for ex_a, ex_b in candidates:
        results.append(evaluate_pair(ex_a, ex_b))
        
    def rank(res):
        c = res['classification'].lower()
        if 'not eligible' in c: return 0
        if 'unlikely' in c: return 1
        if 'likely' in c: return 2
        if 'eligible' in c: return 3
        return -1
        
    results.sort(key=rank, reverse=True)
    best_res = results[0]
    
    if len(results) > 1:
        best_rank = rank(best_res)
        other_res = results[1]
        other_rank = rank(other_res)
        if best_rank > 0 and best_rank == other_rank:
            note = best_res.get('note', '')
            note += f" Skipping Exon {other_res['total_exon_number']} is also {other_res['classification']}."
            best_res['note'] = note.strip()
            
    return best_res

def assess_single_exon(client, original_query, transcript, all_exons, target_exon, vep_entry: Dict[str, Any], overrides: Dict[str, bool] = {}, refseq_id_for_viewer: Optional[str] = None, aso_exists: bool = False):
    # --- Step 1: Data Gathering and Calculations ---
    # Set precision for Decimal calculations
    getcontext().prec = 10

    gene_id = transcript.get('Parent') 
    transcript_id = transcript.get('id')
    protein_id = transcript.get("Translation", {}).get("id")
    protein_meta = client._get(f"/lookup/id/{protein_id}") if protein_id else None
    
    cds_seq = client.get_cds_sequence(transcript_id)
    
    coding_exons = [e for e in all_exons if e['cds_length'] > 0]
    total_coding_exons = len(coding_exons)
    total_cds_len = sum(e['cds_length'] for e in coding_exons)
    if protein_meta and 'length' in protein_meta:
        total_protein_len = protein_meta['length']
    else:
        total_protein_len = total_cds_len / 3 if total_cds_len > 0 else 0

    if not target_exon.get('coding_exon_number'):
        return {"classification": "Unable to Assess", "reason": f"The variant maps to exon {target_exon['total_exon_number']}, which is non-coding."}
    
    chrom, start, end = target_exon['seq_region_name'], target_exon['start'], target_exon['end']
    variants_in_region = client.overlap_region_variation(chrom, start, end)
    
    clinvar_url = f"https://www.ncbi.nlm.nih.gov/clinvar/?term=GRCh38%3A{chrom}%3A{start}-{end}"

    counts = {'missense': 0, 'inframe_del': 0, 'splice': 0, 'nonsense': 0, 'frameshift': 0, 'benign_splice':0}
    splice_variant_links = {"pathogenic": [], "benign": []}
    seen_splice_variants = {"pathogenic": set(), "benign": set()}
    for v in variants_in_region:
        clclass = classify_variant_clinsig(v.get('clinical_significance'))
        conseq = (v.get("consequence_type") or "").lower()
        variant_id = str(v.get('id') or v.get('variation_name') or '').strip()
        variant_source = (v.get('source') or '').strip() or "Ensembl"
        variant_link = _build_variant_link(v) if variant_id else None
        if clclass == "pathogenic":
            if "missense" in conseq: counts['missense'] += 1
            elif "inframe_deletion" in conseq: counts['inframe_del'] += 1
            elif "splice_donor" in conseq or "splice_acceptor" in conseq: counts['splice'] += 1
            elif "stop_gained" in conseq: counts['nonsense'] += 1
            elif "frameshift" in conseq: counts['frameshift'] += 1
            if ("splice_donor" in conseq or "splice_acceptor" in conseq) and variant_id and variant_id not in seen_splice_variants["pathogenic"]:
                splice_variant_links["pathogenic"].append({
                    "id": variant_id,
                    "source": variant_source,
                    "url": variant_link
                })
                seen_splice_variants["pathogenic"].add(variant_id)
        elif clclass == 'benign' and ("splice_donor" in conseq or "splice_acceptor" in conseq):
            counts['benign_splice'] += 1
            if variant_id and variant_id not in seen_splice_variants["benign"]:
                splice_variant_links["benign"].append({
                    "id": variant_id,
                    "source": variant_source,
                    "url": variant_link
                })
                seen_splice_variants["benign"].add(variant_id)
    
    exon_cds_len = target_exon['cds_length']
    exon_aa_len = exon_cds_len / 3
    coding_exon_number = target_exon['coding_exon_number']
    
    # --- Step 1.5: Apply Overrides ---
    overridden_checks = {}
    has_override = False

    # --- Step 2: Condition Checks ---
    # These are the automatically calculated states
    original_conds = {}
    original_conds['cond1_inframe'] = (exon_cds_len % 3 == 0)
    original_conds['cond2_no_stop'] = False
    original_conds['cond10_no_new_codon'] = True
    if cds_seq:
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
    original_conds['cond6_missense'] = counts['missense'] < 3 + counts['frameshift'] + counts['nonsense']
    original_conds['cond7_splice'] = counts['splice'] == 0
    original_conds['cond8_no_inframe_del'] = counts['inframe_del'] == 0
    original_conds['cond9_benign_splice'] = counts['benign_splice'] > 0
    original_conds['aso_exists'] = aso_exists

    # --- Step 2.5: Apply Overrides to get FINAL states ---
    # Start with original values, then overwrite if an override exists
    cond1_inframe = (exon_cds_len % 3 == 0)
    cond2_no_stop = False
    if cds_seq:
        try:
            cds_map, current_pos = {}, 0
            sorted_coding_exons = sorted(coding_exons, key=lambda x: x['coding_exon_number'])
            for exon in sorted_coding_exons:
                cds_map[exon['coding_exon_number']] = cds_seq[current_pos : current_pos + exon['cds_length']]
                current_pos += exon['cds_length']
            skipped_cds = "".join(cds_map[i] for i in sorted(cds_map.keys()) if i != coding_exon_number)
            if skipped_cds:
                prot = str(Seq(skipped_cds).translate(to_stop=False))
                cond2_no_stop = "*" not in prot[:-1]
        except Exception: cond2_no_stop = False

    cond3_not_terminal = (coding_exon_number is not None and coding_exon_number not in (1, total_coding_exons))    
    cond4_small = (Decimal(str(exon_aa_len)) / Decimal(str(total_protein_len))) < Decimal('0.1') if total_protein_len > 0 else False    
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
    cond5_no_domain = domain_count == 0
    cond6_missense = counts['missense'] < 3 + counts['frameshift'] + counts['nonsense']
    cond7_splice = counts['splice'] == 0
    cond8_no_inframe_del = counts['inframe_del'] == 0
    cond9_benign_splice = counts['benign_splice'] > 0

    # --- Step 2.5: Apply Overrides to get FINAL states ---
    # Start with original values, then overwrite if an override exists
    cond1_inframe = overrides.get('Is In-Frame', original_conds['cond1_inframe'])
    cond2_no_stop = overrides.get('No New Stop Codon', original_conds['cond2_no_stop'])
    cond3_not_terminal = overrides.get('Not First/Last Exon', original_conds['cond3_not_terminal'])
    cond4_small = overrides.get('Is <10% of Protein', original_conds['cond4_small'])
    cond5_no_domain = overrides.get('No Domain Overlap', original_conds['cond5_no_domain'])
    cond6_missense = overrides.get('Low Missense Count', original_conds['cond6_missense'])
    cond7_splice = overrides.get('No Pathogenic Splice Variants', original_conds['cond7_splice'])
    cond8_no_inframe_del = overrides.get('No Pathogenic In-Frame Deletions', original_conds['cond8_no_inframe_del'])
    cond9_benign_splice = overrides.get('Benign splice variant found', original_conds['cond9_benign_splice'])
    cond10_no_new_codon = overrides.get('No New Codon Formation', original_conds['cond10_no_new_codon'])
    aso_exists_final = overrides.get('ASO already exists', original_conds['aso_exists'])
    has_override = bool(overrides)
    
    # Check if any override actually changed the outcome
    is_changed = False
    if has_override:
        check_map = {
            'Is In-Frame': cond1_inframe != original_conds['cond1_inframe'], 'No New Stop Codon': cond2_no_stop != original_conds['cond2_no_stop'],
            'Not First/Last Exon': cond3_not_terminal != original_conds['cond3_not_terminal'], 'Is <10% of Protein': cond4_small != original_conds['cond4_small'],
            'No Domain Overlap': cond5_no_domain != original_conds['cond5_no_domain'], 'Low Missense Count': cond6_missense != original_conds['cond6_missense'],
            'No Pathogenic Splice Variants': cond7_splice != original_conds['cond7_splice'], 'No Pathogenic In-Frame Deletions': cond8_no_inframe_del != original_conds['cond8_no_inframe_del'],
            'Benign splice variant found': cond9_benign_splice != original_conds['cond9_benign_splice'],
            'No New Codon Formation': cond10_no_new_codon != original_conds['cond10_no_new_codon'],
            'ASO already exists': aso_exists_final != original_conds['aso_exists']
        }
        is_changed = any(check_map.get(check_name, False) for check_name in overrides)

    # --- Step 3: Classification Logic Chain ---
    if aso_exists_final:
        classification, reason = "Eligible", "An ASO targeting this exon already exists."
    elif not cond3_not_terminal:
        classification, reason = "Not Eligible", "Exon is the first or last coding exon."
    else:
        # Check consecutive exons if originally out of frame, regardless of override status of 'Is In-Frame'
        if not original_conds['cond1_inframe']:
            consecutive_res = assess_consecutive_exons(client, transcript, all_exons, target_exon, vep_entry, refseq_id_for_viewer, overrides)
            if consecutive_res: return consecutive_res

        if not cond1_inframe:
            classification, reason = "Not Eligible", "Exon is out-of-frame, which would disrupt the reading frame. Multi-exon skipping was not found to be a viable in-frame alternative."
        elif not cond2_no_stop:
            classification, reason = "Not Eligible", "Skipping this exon is predicted to create a premature stop codon."
        elif not cond7_splice:
            classification, reason = "Not Eligible", f"Exon contains {counts['splice']} pathogenic splice variant(s), indicating exon loss is pathogenic."
        elif not cond4_small and domain_count > 1:
            classification, reason = "Not Eligible", f"Exon overlaps with {domain_count} protein domains and constitutes >10% of the coding region."
        elif not cond8_no_inframe_del:
            classification, reason = "Not Eligible", f"Exon contains {counts['inframe_del']} pathogenic in-frame deletion(s)."
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
    note_parts.append("Please note that some exceptions for exon skipping apply. Please refer to the methods section for more details.")
    
    note = " ".join(note_parts) if note_parts else None

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
    except Exception as e:
        import traceback; traceback.print_exc()
        visualization_data = None
        
    # --- FINAL RETURN STATEMENT ---
    return {
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
        "splice_variant_links": splice_variant_links,
        "coding_exon_number": coding_exon_number,
        "total_exon_number": target_exon['total_exon_number'], 
        "gene_id": gene_id,
        "transcript_id": transcript_id,
        "clinvar_url": clinvar_url,
        "checks": {
            "ASO already exists": {"passed": aso_exists_final, "original_passed": original_conds['aso_exists'], "overridden": 'ASO already exists' in overrides},
            "Benign splice variant found": {"passed": cond9_benign_splice, "original_passed": original_conds['cond9_benign_splice'], "overridden": 'Benign splice variant found' in overrides},
            "Is In-Frame": {"passed": cond1_inframe, "original_passed": original_conds['cond1_inframe'], "overridden": 'Is In-Frame' in overrides},
            "No New Stop Codon": {"passed": cond2_no_stop, "original_passed": original_conds['cond2_no_stop'], "overridden": 'No New Stop Codon' in overrides},
            "Not First/Last Exon": {"passed": cond3_not_terminal, "original_passed": original_conds['cond3_not_terminal'], "overridden": 'Not First/Last Exon' in overrides},
            "No Pathogenic Splice Variants": {"passed": cond7_splice, "original_passed": original_conds['cond7_splice'], "overridden": 'No Pathogenic Splice Variants' in overrides},
            "No Pathogenic In-Frame Deletions": {"passed": cond8_no_inframe_del, "original_passed": original_conds['cond8_no_inframe_del'], "overridden": 'No Pathogenic In-Frame Deletions' in overrides},
            "No Domain Overlap": {"passed": cond5_no_domain, "original_passed": original_conds['cond5_no_domain'], "overridden": 'No Domain Overlap' in overrides},
            "No New Codon Formation": {"passed": cond10_no_new_codon, "original_passed": original_conds['cond10_no_new_codon'], "overridden": 'No New Codon Formation' in overrides},
            "Low Missense Count": {"passed": cond6_missense, "original_passed": original_conds['cond6_missense'], "overridden": 'Low Missense Count' in overrides},
            "Is <10% of Protein": {"passed": cond4_small, "original_passed": original_conds['cond4_small'], "overridden": 'Is <10% of Protein' in overrides}
        },
        "visualization": visualization_data
    }

def process_single_variant(query: str, client: EnsemblClient, splice_user_input: Optional[str] = None, moa_user_input: Optional[str] = None, exonic_pathogenic_user_input: Optional[str] = None, exon_skipping_overrides: Dict = {}, knockdown_overrides: Dict = {}, wt_upregulation_overrides: Dict = {}, gene_characteristics_overrides: Optional[Dict] = None) -> Dict[str, Any]:
    """
    Contains the complete assessment logic for a single variant query.
    This version is more robust and handles potential unpacking errors.
    """
    try:
        # --- 1. VEP and Consequence Selection ---
        parsed_output = parse_hgvs_query(query)
        if not isinstance(parsed_output, tuple) or len(parsed_output) != 2:
            return {"classification": "Error", "reason": f"Could not parse the input query: '{query}'. Please check the format."}
            
        hgvs_query, gene_symbol_from_query = parsed_output
        if not hgvs_query:
            return {"classification": "Error", "reason": "Invalid input format. Please use a recognized HGVS format (e.g., 'GENE c.123A>G')."}

        vep_data = client.vep_hgvs(hgvs_query)
        if not vep_data or not isinstance(vep_data, list):
            return {"classification": "Unable to Assess", "reason": f"VEP analysis failed for '{hgvs_query}'. The variant may be invalid or not found."}
        
        vep_entry = vep_data[0]
        all_consequences = vep_entry.get('transcript_consequences', [])
        target_consequence = choose_best_consequence(all_consequences, gene_symbol_from_query=gene_symbol_from_query)
        
        if not target_consequence:
            reason = "VEP did not return a consequence."
            if gene_symbol_from_query:
                reason += f" No valid consequence was found for the specified gene '{gene_symbol_from_query}'."
            return {"classification": "Unable to Assess", "reason": reason}

        gene_symbol = target_consequence['gene_symbol']
        definitive_transcript_id = target_consequence['transcript_id']
        gene_id = target_consequence.get('gene_id')

        # --- 2. Get RefSeq ID (for viewer) ---
        refseq_id_for_viewer = None
        mane_consequence = next((c for c in all_consequences if c.get('mane_select')), None)
        if mane_consequence:
            refseq_match = re.search(r'(NM_[0-9]+\.[0-9]+)', mane_consequence['mane_select'])
            if refseq_match:
                refseq_id_for_viewer = refseq_match.group(1)
        if not refseq_id_for_viewer:
            for c in all_consequences:
                if c.get('transcript_id', '').startswith('NM_'):
                    refseq_id_for_viewer = c['transcript_id']
                    break

        # --- 3. Initialize Result & Gene Characteristics ---
        gene_characteristics = get_gene_characteristics(gene_symbol)
        if gene_characteristics_overrides:
            gene_characteristics.update(gene_characteristics_overrides)

        protein_effect = None

        # Attempt 1: Check the chosen target_consequence directly
        if 'hgvsp' in target_consequence:
            protein_effect = target_consequence['hgvsp'].split(':')[-1]
            protein_effect = target_consequence['hgvsp']
        print(protein_effect)
        
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
                **gene_characteristics
            },
            "assessments": {}
        }
        
        # --- 4. N1C Assessed Variants (curated) Check (Exit early if matched) ---
        assessed_match = check_n1c_assessed_variants(gene_symbol, hgvs_query)
        if assessed_match:
            final_result["assessments"]["N1C_Assessed_Variants"] = assessed_match
            return final_result

        # --- 4b. N1C Registry Check ---
        n1c_result = check_n1c_registry(gene_symbol, query, hgvs_query)
        if n1c_result:
            final_result["assessments"]["N1C_Registry_Check"] = n1c_result
            if not n1c_result.get("continue_assessment", False):
                return final_result

        # --- 5. Gene-Level Strategies (Knockdown, WT Upregulation) ---
        # Determine resolved MoA: user selection takes precedence; otherwise a single known MoA
        # Force user selection for MoA; do not auto-resolve from database
        resolved_moa = moa_user_input if moa_user_input in ("GoF", "LoF", "DN") else None
        # Expose resolved_moa for frontend rendering while preserving original Known Mechanism
        final_result["summary"]["resolved_moa"] = resolved_moa
        if resolved_moa:
            moi = set(gene_characteristics.get("moi", []))
            
            if resolved_moa in("GoF", "DN"):
                is_ad = any("autosomal dominant" in str(m).lower() or "x-linked dominant" in str(m).lower() for m in moi)
                if is_ad:
                    n1c_gene_kd = n1c_gene_knockdown_entry(gene_symbol)
                    if n1c_gene_kd:
                        final_result["assessments"]["N1C_Gene_Knockdown"] = n1c_gene_kd
                        return final_result  
                                        
                    # Only run regular knockdown if no N1C gene-level knockdown entry
                    knockdown_assessment = assess_knockdown(gene_characteristics, overrides=knockdown_overrides)
                    if resolved_moa == "DN":
                        knockdown_assessment["note"] = "Dominant Negative mechanism; allele-specific knockdown is necessary."

                    if resolved_moa == "GoF":
                        knockdown_assessment["note"] = "autosomal dominant mode of inheritance; allele-specific knockdown is recommended."
                    final_result["assessments"]["Knockdown"] = knockdown_assessment

            if resolved_moa in("GoF"):

                is_xlinked_or_ar = any("x-linked" in str(m).lower() or "autosomal recessive" in str(m).lower() for m in moi)
                if is_xlinked_or_ar:
                    n1c_gene_kd = n1c_gene_knockdown_entry(gene_symbol)
                    if n1c_gene_kd:
                        final_result["assessments"]["N1C_Gene_Knockdown"] = n1c_gene_kd
                        return final_result
                    
                    knockdown_assessment = assess_knockdown(gene_characteristics, overrides=knockdown_overrides)
                    final_result["assessments"]["Knockdown"] = knockdown_assessment
            

            
            elif resolved_moa == "LoF":
                # Show WT Upregulation only when LoF and Autosomal Dominant MOI (case-insensitive match)
                is_lof_ad = any("autosomal dominant inheritance" in str(m).lower() for m in moi)
                if is_lof_ad:
                    final_result["assessments"]["WT_Upregulation"] = assess_wt_upregulation(client, gene_id, gene_symbol)
                is_x_linked_dominant = any("x-linked dominant inheritance" in str(m).lower() for m in moi)
                if is_x_linked_dominant:
                    wt_upregulation_assessment = assess_wt_upregulation(client, gene_id, gene_symbol)
                    wt_upregulation_assessment["note"] = "WT-upregulation only possible if more than one X-Chromosome is available"
                    final_result["assessments"]["WT_Upregulation"] = wt_upregulation_assessment
                    final_result["summary"]["note"] = "WT-upregulation only possible if more than one X-Chromosome is available"

        # --- 6. Variant-Specific Strategies (Splice & Exon Skipping) ---
        
        # Define variant type
        consequence_terms = set(target_consequence.get('consequence_terms', []))
        exonic_terms = {'missense_variant', 'stop_gained', 'frameshift_variant', 'synonymous_variant', 'inframe_deletion', 'inframe_insertion','splice_donor_variant', 'splice_acceptor_variant'}
        splice_terms = {'splice_region_variant', }
        is_exonic = any(term in consequence_terms for term in exonic_terms)
        is_splice_region = any(term in consequence_terms for term in splice_terms)
        
        # Run Splice Switching Assessment
        variant_identifier_from_vep = vep_entry.get('input')
        if variant_identifier_from_vep:
            splice_assessment = None
            if splice_user_input == 'yes':
                details = {"Confirmation Method": "User-provided validation (qPCR/RNA-seq/cDNA)"}
                splice_assessment = _evaluate_splice_variant_position(variant_identifier_from_vep, vep_entry, details, exonic_pathogenic_user_input=exonic_pathogenic_user_input)
            elif splice_user_input == 'no':
                splice_assessment = {"classification": "Unable to Assess", "reason": "User confirmed no known splice-altering effect."}
            else:
                splice_assessment = assess_splice_switching(variant_identifier_from_vep, vep_entry, gene_symbol, exonic_pathogenic_user_input=exonic_pathogenic_user_input)
            
            if splice_assessment:
                final_result["assessments"]["Splice_Switching"] = splice_assessment

        # Run Exon Skipping Assessment *if* variant is exonic/splice
        if is_exonic:
            exon_skip_assessment_added = False
            transcript_data = client.lookup_id_expand(definitive_transcript_id)
            if transcript_data:
                all_exons = extract_exons_from_transcript(transcript_data)
                v_start, v_end = vep_entry['start'], vep_entry['end']
                target_exon = next((ex for ex in all_exons if ex['seq_region_name'] == vep_entry['seq_region_name'] and max(v_start, ex['start']) <= min(v_end, ex['end'])), None)
                
                if target_exon:                    
                    # N1C registry exon-skipping support: if N1C lists exon skipping for this exon, mark eligible
                    try:
                        n1c_exons, n1c_links, n1c_exon_link_map = n1c_exon_skipping_exon_numbers_for_gene(gene_symbol)
                    except Exception:
                        n1c_exons, n1c_links, n1c_exon_link_map = set(), [], {}
                    try:
                        n1c_variant_exon_map = n1c_exon_skipping_variant_exon_map(gene_symbol, all_exons, client)
                    except Exception:
                        n1c_variant_exon_map = {}

                    exon_num = target_exon.get('total_exon_number')
                    matched_links: List[str] = []
                    variant_links = n1c_variant_exon_map.get(exon_num, []) if exon_num is not None else []
                    has_text_match = exon_num in n1c_exons if exon_num is not None else False
                    if variant_links:
                        matched_links.extend(variant_links)
                    elif has_text_match:
                        matched_links.extend(n1c_exon_link_map.get(exon_num, []))
                        if not matched_links and n1c_links:
                            matched_links.extend(n1c_links)

                    aso_exists_for_exon_skip = bool(variant_links or has_text_match)

                    exon_skip_result = assess_single_exon(client, query, transcript_data, all_exons, target_exon, vep_entry, overrides=exon_skipping_overrides, refseq_id_for_viewer=refseq_id_for_viewer, aso_exists=aso_exists_for_exon_skip)
                    if "visualization" in exon_skip_result and exon_skip_result["visualization"]:
                        final_result["visualization"] = exon_skip_result.pop("visualization")

                    if variant_links or has_text_match:
                        exon_skip_result = dict(exon_skip_result)
                        note = f"N1C registry contains a project that maps to exon {exon_num} in {gene_symbol}. The input variant lies in the same exon and is therefore considered skippable."
                        prev_reason = exon_skip_result.get('reason', '')
                        exon_skip_result['reason'] = (prev_reason + ' ' + note).strip()
                        exon_skip_result['classification'] = 'Eligible'
                        det = exon_skip_result.get('details', {}) if isinstance(exon_skip_result.get('details'), dict) else {}
                        if matched_links:
                            det['N1C Exon Skipping (same exon)'] = ', '.join(sorted(set(matched_links)))
                        elif n1c_links:
                            det['N1C Exon Skipping Entries'] = ', '.join(n1c_links)
                        det['Supported Exon'] = str(exon_num)
                        exon_skip_result['details'] = det
                    final_result["assessments"]["Exon_Skipping"] = exon_skip_result
                    exon_skip_assessment_added = True
            
            if not exon_skip_assessment_added:
                # Add this block if the logic above fails to add an assessment
                final_result["assessments"]["Exon_Skipping"] = {
                    "classification": "Unable to Assess",
                    "reason": "Variant is exonic or in a splice region, but the target exon could not be determined (e.g., VEP/Ensembl data issue)."
                }

        # 7. Final Fallback
        # This will now only trigger for non-exonic, non-splice, non-gene-strategy variants.
        if not final_result["assessments"]:
            final_result["assessments"]["General_Assessment"] = { 
                "classification": "Unable to Assess", 
                "reason": "Could not determine a primary ASO strategy. The variant is not exonic and no gene-level strategies (Knockdown, WT Upregulation) are applicable." 
            }
        
        return final_result

    except Exception as e:
        import traceback; traceback.print_exc()
        return {"classification": "Error", "reason": f"An unexpected server error occurred: {str(e)}"}

# --- Main Flask Routes ---
app = Flask(__name__)

load_databases()

@app.route('/')
def index(): return render_template('index.html', title="Tool")
@app.route('/usage')
def usage(): return render_template('usage.html', title="Usage")
@app.route('/about')
def about(): return render_template('about.html', title="About/Methods")
@app.route('/cite')
def cite(): return render_template('cite.html', title="How to Cite")

@app.route('/api_docs')
def api_docs():
    """Serves the API documentation page."""
    return render_template('api_docs.html', title="API Documentation")

@app.route('/api/v1/assess', methods=['POST'])
def api_assess():
    """
    Handles a single variant assessment via a POST request for programmatic access.
    Returns the full assessment data as JSON.
    """
    data = request.get_json()
    if not data or 'query' not in data:
        return jsonify({"error": "Request body must be JSON with a 'query' key."}), 400

    query = data['query']

    client = EnsemblClient()
    result = process_single_variant(query, client)
    
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
    data = request.get_json()
    if not data or 'query' not in data:
        return jsonify({"classification": "Error", "reason": "No query provided."}), 400
    query = data['query']
    splice_input = data.get('splice_user_input', None)
    exonic_pathogenic_user_input = data.get('exonic_pathogenic_user_input', None)
    moa_input = data.get('moa_user_input', None)
    es_overrides = data.get('exon_skipping_overrides', {})
    kd_overrides = data.get('knockdown_overrides', {})
    wt_overrides = data.get('wt_upregulation_overrides', {})
    gene_char_overrides = data.get('gene_characteristics_overrides', None)
    client = EnsemblClient()
    result = process_single_variant(query, client, splice_user_input=splice_input, moa_user_input=moa_input, exonic_pathogenic_user_input=exonic_pathogenic_user_input, exon_skipping_overrides=es_overrides, knockdown_overrides=kd_overrides, wt_upregulation_overrides=wt_overrides, gene_characteristics_overrides=gene_char_overrides)
    
    return jsonify(result)
@app.route('/batch_assess', methods=['POST'])
def batch_assess():
    if 'file' not in request.files:
        return jsonify({"error": "No file part"}), 400
    file = request.files['file']
    if file.filename == '':
        return jsonify({"error": "No selected file"}), 400

    # Helper to normalize Y/N
    def norm_yn(val):
        if pd.isna(val): return None
        s = str(val).strip().upper()
        if s.startswith('Y'): return 'yes'
        if s.startswith('N'): return 'no'
        return None

    # Helper to normalize MoA
    def norm_moa(val):
        if pd.isna(val): return None
        s = str(val).strip().upper()
        if 'GOF' in s: return 'GoF'
        if 'LOF' in s: return 'LoF'
        if 'DN' in s: return 'DN'
        return None

    # Helper to normalize MOI
    def norm_moi(val):
        if pd.isna(val): return None
        s = str(val).strip().upper()
        res = []
        if 'AD' in s: res.append('Autosomal dominant inheritance')
        if 'AR' in s: res.append('Autosomal recessive inheritance')
        if 'XD' in s: res.append('X-linked dominant inheritance')
        if 'XR' in s: res.append('X-linked recessive inheritance')
        return res if res else None

    # Helper to normalize Haplo
    def norm_haplo(val):
        if pd.isna(val): return None
        s = str(val).strip().upper()
        if s.startswith('Y'): return {"text": "Sufficient evidence", "url": None}
        if s.startswith('N'): return {"text": "No evidence", "url": None}
        return None

    try:
        if file.filename.endswith('.xlsx'):
            df = pd.read_excel(file, header=None)
        else:
            df = pd.read_csv(file, header=None)
        
        # Check first row values to detect template
        first_row = df.iloc[0].astype(str).str.strip().tolist()
        is_template = 'HGVS_Variant' in first_row
        
        if is_template:
            df.columns = first_row
            df = df.iloc[1:].reset_index(drop=True)
            rows_to_process = df.to_dict('records')
        else:
            variants = df[0].dropna().astype(str).tolist()
            rows_to_process = [{'HGVS_Variant': v} for v in variants]

    except Exception as e:
        return jsonify({"error": f"Error reading file: {e}"}), 400

    client = EnsemblClient()
    output_rows = []

    for row_data in rows_to_process:
        variant = str(row_data.get('HGVS_Variant', '')).strip()
        if not variant or variant.lower() == 'nan': continue

        # Capture raw inputs for reporting
        raw_moa = row_data.get('Pathomechanism(GOF/LOF/DN)')
        raw_splice = row_data.get('Splicing effects (Y/N)')
        raw_exonic = row_data.get('Splicing-independent effects(Y/N)')
        raw_moi = row_data.get('Mode of Inheritance (AD/AR/XD/XR)')
        raw_haplo = row_data.get('Haploinsufficiency(Y/N)')
        
        user_inputs = []
        if not pd.isna(raw_moa): user_inputs.append(f"MoA: {raw_moa}")
        if not pd.isna(raw_splice): user_inputs.append(f"Splice: {raw_splice}")
        if not pd.isna(raw_exonic): user_inputs.append(f"Splice-Indep: {raw_exonic}")
        if not pd.isna(raw_moi): user_inputs.append(f"MOI: {raw_moi}")
        if not pd.isna(raw_haplo): user_inputs.append(f"Haplo: {raw_haplo}")
        user_input_str = "; ".join(user_inputs) if user_inputs else "None"

        # Extract overrides
        moa_input = norm_moa(raw_moa)
        splice_input = norm_yn(raw_splice)
        exonic_pathogenic_input = norm_yn(raw_exonic)
        
        moi_override = norm_moi(raw_moi)
        haplo_override = norm_haplo(raw_haplo)

        gene_char_overrides = {}
        if moi_override: gene_char_overrides['moi'] = moi_override
        if haplo_override: gene_char_overrides['haploinsufficiency'] = haplo_override

        result = process_single_variant(
            variant, client, 
            splice_user_input=splice_input, 
            moa_user_input=moa_input, 
            exonic_pathogenic_user_input=exonic_pathogenic_input,
            gene_characteristics_overrides=gene_char_overrides
        )
        
        row = {"Variant": variant}
        row["User Input"] = user_input_str
        summary = result.get("summary", {})
        assessments = result.get("assessments", {})

        row["Gene"] = summary.get("gene", "N/A")
        row["ClinGen MOI"] = ', '.join(summary.get("moi", []))
        row["ClinGen MOA"] = ', '.join(summary.get("moa", []))
        # Ensembl Transcript ID and link
        transcript_id = summary.get("transcript_id")
        row["Ensembl Transcript"] = transcript_id or "N/A"
        row["Ensembl Transcript Link"] = (
            f"https://www.ensembl.org/Homo_sapiens/Transcript/Summary?t={transcript_id}" if transcript_id else "N/A"
        )
        
        haplo_info = summary.get("haploinsufficiency", {})
        row["Haploinsufficiency"] = haplo_info.get("text", "N/A")
        row["ClinGen Link"] = haplo_info.get("url", "N/A")

        # 1. Assess if an ASO exists and add the N1C link(s)
        n1c_registry = assessments.get("N1C_Registry_Check", {}) or {}
        n1c_assessed = assessments.get("N1C_Assessed_Variants", {}) or {}
        if n1c_registry or n1c_assessed:
            row["Existing ASO (N1C)"] = "Yes"
            row["N1C Registry Link"] = n1c_registry.get("link", "N/A")
            row["N1C Assessed (Curated) Link"] = n1c_assessed.get("link", "N/A")
        else:
            row["Existing ASO (N1C)"] = "No"
            row["N1C Registry Link"] = "N/A"
            row["N1C Assessed (Curated) Link"] = "N/A"
            
        # 2. Get the Antisense Transcript ID for WT Upregulation
        wt_up = assessments.get("WT_Upregulation", {})
        antisense_ids = wt_up.get("antisense_gene_ids", [])
        row["Antisense Transcript ID"] = ", ".join(antisense_ids) if antisense_ids else "N/A"
        # Exon Skipping
        skip = assessments.get("Exon_Skipping", {})
        row["Exon Skipping Assessment"] = skip.get("classification", "NA")
        for check, status in skip.get("checks", {}).items():
            if isinstance(status, dict):
                row[f"ES Check: {check}"] = status.get("original_passed", "N/A")
            else:
                row[f"ES Check: {check}"] = status
        # Ensembl exon view and domains (if available)
        if skip:
            gid = skip.get("gene_id")
            tid = skip.get("transcript_id")
            if gid and tid:
                row["Ensembl Exon View Link"] = f"https://www.ensembl.org/Homo_sapiens/Transcript/Exons?db=core;g={gid};t={tid}"
            else:
                row["Ensembl Exon View Link"] = "N/A"
            domain_names = skip.get("domain_names") or []
            domain_details = skip.get("domain_details") or []
            if domain_details:
                formatted_domains = []
                for d in domain_details:
                    label = d.get("label") or d.get("name") or "Domain"
                    url = d.get("url")
                    formatted_domains.append(f"{label} [{url}]" if url else label)
                row["Domains"] = "; ".join(formatted_domains)
            else:
                row["Domains"] = ", ".join(domain_names) if domain_names else "N/A"

        # Splice Correction
        splice = assessments.get("Splice_Switching", {})
        splice_class = splice.get("classification", "Unable to Assess")
        if splice_class == "Awaiting User Input":
            splice_class = "Splicing-independent effects unknown. Manual/Single Variant processing needed"
        row["Splice Correction Assessment"] = splice_class
        row["Splicing Validation DOI"] = splice.get("details", {}).get("Publication DOI", "NA")
        # Splicing DB/Source links if available
        splice_details = splice.get("details", {}) if isinstance(splice.get("details", {}), dict) else {}
        splicing_db_link = splice_details.get("SSCVDB Gene Page") or splice_details.get("Publication") or "N/A"
        row["Splicing DB Link"] = splicing_db_link

        # WT Upregulation and Knockdown (assessments remain)
        row["WT-Upregulation"] = wt_up.get("classification", "NA")
        row["Knockdown"] = assessments.get("Knockdown", {}).get("classification", "NA")

        # Manual validation needs and dual MoA assessment when MoA unresolved
        manual_needs = []
        summary_moa_list = summary.get("moa", []) or []
        resolved_moa = summary.get("resolved_moa")
        # Splice manual validation prompt
        if splice.get("user_validation_prompt") or (splice.get("classification") in ("Not in Database", "Unable to Assess")):
            manual_needs.append("Splice validation (Variant was not found in SpliceVarDB/SSCVDB and therefore requires user confirmation)")
        # MoA unclear -> warn
        if not resolved_moa and (len(summary_moa_list) != 1):
            manual_needs.append("Mechanism selection (DN vs GoF vs LoF)")

        row["Manual Validations Needed"] = "; ".join(manual_needs) if manual_needs else "None"

        # Overall Eligibility: highest across all assessment classifications
        def _normalize_class(c: Optional[str]) -> str:
            if not c:
                return "Unable to Assess"
            s = str(c).strip().lower().replace('-', ' ')
            if 'not eligible' in s:
                return 'Not Eligible'
            if 'likely eligible' in s:
                return 'Likely Eligible'
            if 'unlikely eligible' in s:
                return 'Unlikely Eligible'
            if 'eligible' in s:
                return 'Eligible'
            if 'unable to assess' in s or 'not in database' in s:
                return 'Unable to Assess'
            return 'Unable to Assess'

        rank_order = {
            'Not Eligible': 1,
            'Unable to Assess': 2,
            'Unlikely Eligible': 3,
            'Likely Eligible': 4,
            'Eligible': 5,
        }
        best_label = 'Unable to Assess'
        best_score = 0
        for akey, aval in assessments.items():
            if not isinstance(aval, dict):
                continue
            label = _normalize_class(aval.get('classification'))
            score = rank_order.get(label, 2)
            if score > best_score:
                best_score = score
                best_label = label
        row["Overall Eligibility"] = best_label
        
        output_rows.append(row)

    # --- Create and send the Excel file (no changes needed below this line) ---
    if not output_rows:
        return jsonify({"error": "No variants found in file."}), 400
        
    output_df = pd.DataFrame(output_rows)
    output = io.BytesIO()
    with pd.ExcelWriter(output, engine='openpyxl') as writer:
        output_df.to_excel(writer, index=False, sheet_name='AVEC_Batch_Results')
    output.seek(0)
    
    return send_file(
        output,
        as_attachment=True,
        download_name='avec_batch_results.xlsx',
        mimetype='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'
    )

@app.route('/download_batch_template')
def download_batch_template():
    """Serves the batch processing template file."""
    try:
        BASE_DIR = os.path.dirname(os.path.abspath(__file__))
        DATA_DIR = os.path.join(BASE_DIR, 'data')
        template_path = os.path.join(DATA_DIR, 'batch_template.xlsx')
        return send_file(
            template_path,
            as_attachment=True,
            download_name='AVEC_batch_template.xlsx',
            mimetype='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'
        )
    except FileNotFoundError:
        return "Template file not found on server.", 404
