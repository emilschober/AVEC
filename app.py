import os
from flask import Flask, request, jsonify, render_template, send_file
import time
import re
import requests
from functools import lru_cache
from typing import Dict, Any, Optional, Tuple, List
from Bio.Seq import Seq
import pandas as pd
import io

# ===== CONFIGURATION =====
ENSEMBL_REST = "https://rest.ensembl.org"
HEADERS = {"Content-Type": "application/json", "Accept": "application/json"}
N1C_API_URL = "https://gene-registry.onrender.com/api/data?table=N1C_projects" 
# Global DataFrames to be loaded at startup
clingen_df: Optional[pd.DataFrame] = None
goflof_df: Optional[pd.DataFrame] = None
splicevar_df: Optional[pd.DataFrame] = None
n1c_variants_df: Optional[pd.DataFrame] = None 

def load_templates():
    with open('templates/base.html', 'w', encoding='utf-8') as f: f.write(base_html)
    with open('templates/index.html', 'w', encoding='utf-8') as f: f.write(index_html)
    with open('templates/about.html', 'w', encoding='utf-8') as f: f.write(about_html)
    with open('templates/cite.html', 'w', encoding='utf-8') as f: f.write(cite_html)
    with open('templates/api_docs.html', 'w', encoding='utf-8') as f: f.write(api_docs_html)

load_templates()

# --- Data Loading ---
def load_databases():
    """
    Loads all necessary data files and fetches N1C registry data
    into global pandas DataFrames.
    """
    global clingen_df, goflof_df, splicevar_df, n1c_variants_df
    print("Loading databases...")
    try:
        clingen_df = pd.read_csv"Clingen-Curation-Activity-Summary-Report-2025-10-15.csv").set_index('gene_symbol')
        goflof_df = pd.read_csv("goflof_HGMD2019_v032021_allfeat.csv").set_index('GENE')
        
        # Load SpliceVarDB from Excel
        splicevar_df = pd.read_excel("splicevardb.xlsx")
        
        # Sanitize SpliceVarDB data (critical for lookups)
        splicevar_df.columns = splicevar_df.columns.str.strip()
        if 'hgvs' in splicevar_df.columns and 'gene' in splicevar_df.columns:
            splicevar_df['hgvs'] = splicevar_df['hgvs'].astype(str).str.strip()
            splicevar_df['gene'] = splicevar_df['gene'].astype(str).str.strip()
        print(f"Loaded and sanitized SpliceVarDB data for {len(splicevar_df)} variants.")

        # Fetch and load N1C N-of-1 Projects data
        print("Fetching N=1 Collaborative Projects registry data...")
        response = requests.get(N1C_API_URL, timeout=30)
        response.raise_for_status() # Will raise an error if the request fails
        n1c_data = response.json()
        n1c_variants_df = pd.DataFrame(n1c_data)
        
        # Sanitize the N1C columns we will search on
        if 'Gene' in n1c_variants_df.columns:
            n1c_variants_df['Gene'] = n1c_variants_df['Gene'].astype(str).str.strip()
        if 'Coding DNA change (c.)' in n1c_variants_df.columns:
            n1c_variants_df['Coding DNA change (c.)'] = n1c_variants_df['Coding DNA change (c.)'].astype(str).str.strip()
        print(f"Successfully loaded {len(n1c_variants_df)} projects from the N1C registry.")

    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"An error occurred during database loading: {e}")
        exit(1)

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
    def vep_hgvs(self, hgvs_string): return self._get(f"/vep/human/hgvs/{hgvs_string.strip()}", params={'variant_class': 1})
    def get_cds_sequence(self, transcript_id):
        data = self._get(f"/sequence/id/{transcript_id}", params={"type": "cds"})
        return data.get("seq") if isinstance(data, dict) else None
    def get_domains(self, protein_id):
        all_features = self._get(f"/overlap/translation/{protein_id}", params={"feature": "protein_feature"})
        if not all_features or not isinstance(all_features, list): return []
        domain_sources = {'CDD','Pfam','SMART','PROSITE profiles','PROSITE patterns','SUPERFAMILY','PRINTS','TIGRFAM','ProDom'}
        preliminary_domains = [f for f in all_features if f.get('type') in domain_sources]
        unique_interpro_domains = {f['interpro']: f for f in preliminary_domains if f.get('interpro')}
        return list(unique_interpro_domains.values())
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
    
# --- Helper & Parsing Functions ---

def _evaluate_splice_variant_position(variant_hgvs: str, vep_data: Dict[str, Any], details: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    """
    Assesses a validated splice-altering variant based on its genomic position.
    This logic is shared by database-found variants and user-validated variants.
    """
    core_hgvs_match = re.search(r'(c\..*)', variant_hgvs, re.IGNORECASE)
    if not core_hgvs_match: 
        return None
    core_canonical_hgvs = core_hgvs_match.group(1).lower()

    result = {"details": details}
    consequence_terms = set(vep_data.get('transcript_consequences', [{}])[0].get('consequence_terms', []))
    
    is_intronic_by_consequence = 'intron_variant' in consequence_terms or 'splice_acceptor_variant' in consequence_terms or 'splice_donor_variant' in consequence_terms
    is_intronic_by_notation = '+' in core_canonical_hgvs or '-' in core_canonical_hgvs
    is_intronic = is_intronic_by_consequence or is_intronic_by_notation

    if is_intronic:
        dist_match = re.search(r'[+-](\d+)', core_canonical_hgvs)
        if not dist_match: 
            return None 
        dist = int(dist_match.group(1))
        
        if '+' in core_canonical_hgvs: # Downstream (e.g., c.123+1G>A)
            if dist <= 5:
                result.update({"classification": "Not Eligible", "reason": "Variant is a validated splice-altering variant located too close to the canonical splice site (<=+5bp)."})
            elif 6 <= dist <= 50:
                result.update({"classification": "Unlikely Eligible", "reason": "Variant is a validated splice-altering variant located near the canonical splice site (+6-+50bp)."})
            else:
                result.update({"classification": "Likely Eligible", "reason": "Variant is a validated splice-altering variant in a favorable deep-intronic position (>+50bp)."})
        
        elif '-' in core_canonical_hgvs: # Upstream (e.g., c.124-2A>G)
            if dist <= 5:
                result.update({"classification": "Not Eligible", "reason": "Variant is a validated splice-altering variant located too close to the canonical splice site (>=-5bp)."})
            elif 6 <= dist <= 100:
                result.update({"classification": "Unlikely Eligible", "reason": "Variant is a validated splice-altering variant located near the canonical splice site (-6-(-100b)p)."})
            else:
                result.update({"classification": "Likely Eligible", "reason": "Variant is a validated splice-altering variant in a favorable deep-intronic position (<-100bp)."})
        
    elif any(c in consequence_terms for c in ['missense_variant', 'synonymous_variant']):
        if 'splice_region_variant' in consequence_terms:
            result.update({"classification": "Not Eligible", "reason": "This validated exonic splice-altering variant is within the canonical splice region, making it high-risk."})
        else:
            result.update({"classification": "Likely Eligible", "reason": "This validated exonic splice-altering variant is outside the immediate splice region, making it a potential candidate for correction."})
    
    else:
        # Fallback if it's splice-altering but not in a recognized position
        result.update({"classification": "Unlikely Eligible", "reason": "Variant is a validated splice-altering variant, but also presumed to cause other effects (e.g. it is a nonsense variant)."})

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
    """
    # Check if the DataFrame was loaded successfully
    if n1c_variants_df is None or n1c_variants_df.empty or not gene_symbol:
        return None

    # Extract the core c. notation from the VEP-formatted HGVS string
    core_hgvs_match = re.search(r'(c\..*)', formatted_hgvs, re.IGNORECASE)
    if not core_hgvs_match:
        return None
    core_hgvs = core_hgvs_match.group(1)

    # Filter the DataFrame for the correct gene (case-insensitive)
    gene_matches = n1c_variants_df[n1c_variants_df['Gene'].str.upper() == gene_symbol.upper()]
    if gene_matches.empty:
        return None
        
    # Search within the gene-specific rows for the variant notation
    for index, row in gene_matches.iterrows():
        # Check if our core HGVS notation is present in the registry's 'Coding DNA change (c.)' field
        registry_c_dot = row.get('Coding DNA change (c.)')
        if registry_c_dot and core_hgvs.lower() in registry_c_dot.lower():
            # Match found! Extract data and return the result.
            status = row.get('Status', 'N/A')
            modality = row.get('Therapeutic Modality', 'N/A')
            therapypublication = row.get('Therapy Publication','N/A')
            n1cid = row.get('ID')
            
            return {
                "classification": "Eligible",
                "reason": (
                    f"<p>A direct match for a variant in <strong>{gene_symbol}</strong> was found in the N=1 Collaborative Projects Registry.</p>"
                    f"<ul>"
                    f"<li><strong>Matched Variant:</strong> {registry_c_dot}</li>"
                    f"<li><strong>Status:</strong> {status}</li>"
                    f"<li><strong>Therapeutic Modality:</strong> {modality}</li>"
                    f"</ul>"
                    f"<a href='https://generegistry.n1collaborative.org/entry.html?id={n1cid}' target='_blank' rel='noopener noreferrer'>Click here to view the N=1 Collaborative registry page.</a>"
                ),
                "link": f"https://generegistry.n1collaborative.org/entry.html?id={n1cid}" 
            }
    return None

def get_gene_characteristics(gene_symbol: str) -> Dict[str, Any]:
    """
    Retrieves MOI, Haploinsufficiency, and MOA from loaded dataframes,
    including source URLs for better explainability.
    """
    characteristics = {
        "moi": [],
        "haploinsufficiency": {"text": "Unknown", "url": None},
        "moa": [],
        "gene_url": None # --- ADDITION: Initialize gene_url ---
    }

    # --- ClinGen Lookup ---
    if clingen_df is not None and gene_symbol in clingen_df.index:
        # Use .loc[gene_symbol] which can return a Series or DataFrame
        gene_data_row = clingen_df.loc[gene_symbol]
        if isinstance(gene_data_row, pd.DataFrame):
            gene_data_row = gene_data_row.iloc[0] # Take the first row if multiple exist

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
            else: haplo_text = "No evidence"
            
            characteristics['haploinsufficiency'] = {
                "text": haplo_text,
                "url": hap_url if pd.notna(hap_url) else None
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

def assess_knockdown(gene_characteristics: Dict[str, Any]) -> Dict[str, Any]:
    """Assesses eligibility for a knockdown strategy."""
    haplo_obj = gene_characteristics.get("haploinsufficiency", {"text": "Unknown"})
    haplo_status_text = haplo_obj.get("text", "Unknown")
    
    moi = gene_characteristics.get("moi", [])
    is_ad_lof = any("Autosomal Dominant" in m for m in moi) and gene_characteristics.get("moa", []) == 'LoF'
    
    checks = {"Gene is not haploinsufficient": haplo_status_text in ["No evidence", "Little evidence"]}
    
    if is_ad_lof or haplo_status_text == "Sufficient evidence":
        return {
            "classification": "Not Eligible", # Changed to be more definitive
            "reason": "Gene is associated with haploinsufficiency or Autosomal Dominant LoF, making knockdown risky.",
            "checks": checks
        }
    else:
        return {
            "classification": "Likely Eligible",
            "reason": "Gene is not known to be sensitive to haploinsufficiency.",
            "checks": checks
        }
    
def get_overlapping_genes(self, gene_id):
        """Fetches all genes that overlap with a given Ensembl Gene ID."""
        data = self._get(f"/overlap/id/{gene_id}", params={"feature": "gene"})
        return data if isinstance(data, list) else []

def assess_wt_upregulation(client: EnsemblClient, gene_id: str, gene_symbol: str) -> Dict[str, Any]:
    """
    Assesses for WT upregulation by checking for overlapping NATs and by
    searching for a conventional antisense gene name ([GENE]-AS1).
    Trusts the '-AS1' naming convention without a strict biotype check.
    """
    if not gene_id or not gene_symbol:
        return {"classification": "Unable to Assess", "reason": "Missing Gene ID or Symbol."}

    found_antisense_genes = {} # Use a dict to store unique NATs by their ID

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

        # --- Evaluate the combined results ---
        if found_antisense_genes:
            nat_list = list(found_antisense_genes.values())
            nat_names = ", ".join([nat.get('external_name', nat['id']) for nat in nat_list])
            nat_ids = [nat['id'] for nat in nat_list]
            
            details = {}
            for nat in nat_list:
                nat_name = nat.get('external_name', nat['id'])
                ensembl_link = f"https://www.ensembl.org/Homo_sapiens/Gene/Summary?g={nat['id']}"
                # The key will be the gene name, and the value will be the link.
                details[nat_name] = ensembl_link
            
            return {
                "classification": "Likely Eligible",
                "reason": f"Found at least one potential natural antisense transcript ({nat_names}) for {gene_symbol}. Targeting this NAT could upregulate the WT allele.",
                "details": details, # The new details dictionary with links
                "antisense_gene_ids": nat_ids,
                "checks": {"Overlapping antisense transcript found": True}
            }
        else:
            return {
                "classification": "Unlikely Eligible",
                "reason": f"No overlapping or conventionally named antisense transcripts were found for {gene_symbol} in Ensembl.",
                "checks": {"Overlapping antisense transcript found": False}
            }

    except Exception as e:
        return {
            "classification": "Unable to Assess",
            "reason": f"An error occurred while searching for antisense transcripts: {e}",
            "checks": {}
        }

def assess_splice_switching(variant_hgvs: str, vep_data: Dict[str, Any], gene_symbol: str) -> Optional[Dict[str, Any]]:
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
        # --- MODIFICATION: Return prompt object if gene not found ---
        return {
            "classification": "Not in Database",
            "reason": "This variant was not found in SpliceVarDB. If there is experimental validation (by qPCR or RNA-seq) that this variant is splice-altering, please confirm below.",
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

            return _evaluate_splice_variant_position(variant_hgvs, vep_data, details)

    # --- Return prompt object if loop finishes with no match ---
    return {
        "classification": "Not in Database",
        "reason": "This variant was not found in SpliceVarDB. If there is experimental validation (by qPCR or RNA-seq) that this variant is splice-altering, please confirm below.",
        "user_validation_prompt": True
    }

def assess_single_exon(client, original_query, transcript, all_exons, target_exon, vep_entry: Dict[str, Any], refseq_id_for_viewer: Optional[str] = None):
    # --- Step 1: Data Gathering and Calculations ---
    gene_id = transcript.get('Parent') 
    transcript_id = transcript.get('id')
    protein_id = transcript.get("Translation", {}).get("id")
    cds_seq = client.get_cds_sequence(transcript_id)
    
    coding_exons = [e for e in all_exons if e['cds_length'] > 0]
    total_coding_exons = len(coding_exons)
    total_cds_len = sum(e['cds_length'] for e in coding_exons)

    if not target_exon.get('coding_exon_number'):
        return {"classification": "Unable to Assess", "reason": f"The variant maps to exon {target_exon['total_exon_number']}, which is non-coding."}
    
    chrom, start, end = target_exon['seq_region_name'], target_exon['start'], target_exon['end']
    variants_in_region = client.overlap_region_variation(chrom, start, end)
    
    clinvar_url = f"https://www.ncbi.nlm.nih.gov/clinvar/?term=GRCh38%3A{chrom}%3A{start}-{end}"

    counts = {'missense': 0, 'inframe_del': 0, 'splice': 0, 'nonsense': 0, 'frameshift': 0, 'benign_splice':0}
    for v in variants_in_region:
        clclass = classify_variant_clinsig(v.get('clinical_significance'))
        conseq = (v.get("consequence_type") or "").lower()
        if clclass == "pathogenic":
            if "missense" in conseq: counts['missense'] += 1
            elif "inframe_deletion" in conseq: counts['inframe_del'] += 1
            elif "splice_donor" in conseq or "splice_acceptor" in conseq: counts['splice'] += 1
            elif "stop_gained" in conseq: counts['nonsense'] += 1
            elif "frameshift" in conseq: counts['frameshift'] += 1
        elif clclass == 'benign' and ("splice_donor" in conseq or "splice_acceptor" in conseq):
            counts['benign_splice'] += 1
    
    exon_cds_len = target_exon['cds_length']
    coding_exon_number = target_exon['coding_exon_number']
    
    # --- Step 2: Condition Checks ---
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
    cond4_small = (exon_cds_len / total_cds_len) < 0.1 if total_cds_len > 0 else False
    
    domains = client.get_domains(protein_id) if protein_id else []
    overlapping_domain_names = []
    if domains:
        cds_pos_start = sum(e['cds_length'] for e in sorted(coding_exons, key=lambda x: x['coding_exon_number']) if e['coding_exon_number'] < coding_exon_number)
        exon_aa_start, exon_aa_end = (cds_pos_start // 3) + 1, ((cds_pos_start + exon_cds_len -1) // 3) + 1
        for d in domains:
            if d.get('start', 0) <= exon_aa_end and d.get('end', 0) >= exon_aa_start:
                overlapping_domain_names.append(d.get('description', d.get('id', 'Unknown Domain')))

    domain_count = len(overlapping_domain_names)
    cond5_no_domain = domain_count == 0
    cond6_missense = counts['missense'] < 3 + counts['frameshift'] + counts['nonsense']
    cond7_splice = counts['splice'] == 0
    cond8_no_inframe_del = counts['inframe_del'] == 0
    cond9_benign_splice = counts['benign_splice'] > 0
    
    # --- Step 3: Classification Logic Chain ---
    classification, reason = "Undetermined", ""
    if cond9_benign_splice:
        classification, reason = "Eligible", "Exon contains benign splice variants, suggesting it may be safely skipped."
    elif not cond3_not_terminal:
        classification, reason = "Not Eligible", "Exon is the first or last coding exon."
    elif not cond1_inframe:
        classification, reason = "Not Eligible", "Exon is out-of-frame, which would disrupt the reading frame."
    elif not cond2_no_stop:
        classification, reason = "Not Eligible", "Skipping this exon is predicted to create a premature stop codon."
    elif not cond7_splice:
        classification, reason = "Not Eligible", f"Exon contains {counts['splice']} pathogenic splice variant(s), indicating exon loss is pathogenic."
    elif not cond4_small and domain_count > 1:
        classification, reason = "Not Eligible", f"Exon overlaps with {domain_count} protein domains and constitutes >10% of the coding region."
    elif not cond8_no_inframe_del:
        classification, reason = "Not Eligible", f"Exon contains {counts['inframe_del']} pathogenic in-frame deletion(s)."
    elif not cond6_missense:
        classification, reason = "Unlikely Eligible", f"Exon is a mutational hotspot with {counts['missense']} pathogenic missense variants."
    elif not cond4_small:
        classification, reason = "Unlikely Eligible", "Exon constitutes >=10% of the protein, risking major functional loss."
    elif not cond5_no_domain:
        classification, reason = "Unlikely Eligible", f"Exon overlaps with {domain_count} protein domain(s)."
    else:
        classification, reason = "Likely Eligible", "Exon meets the primary criteria for a skippable exon."

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
                            "name": domain.get('description', domain.get('id', 'Domain'))
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
        "frac_cds": f"{(exon_cds_len / total_cds_len * 100):.2f}%" if total_cds_len > 0 else "N/A",
        "pathogenic_variant_counts": counts, 
        "domain_count": domain_count,
        "domain_names": overlapping_domain_names,
        "coding_exon_number": coding_exon_number,
        "total_exon_number": target_exon['total_exon_number'], 
        "gene_id": gene_id,
        "transcript_id": transcript_id,
        "clinvar_url": clinvar_url,
        "checks": {
            "Benign splice variant found": cond9_benign_splice, "Is In-Frame": cond1_inframe, 
            "No New Stop Codon": cond2_no_stop, "Not First/Last Exon": cond3_not_terminal,
            "No Pathogenic Splice Variants": cond7_splice, "No Pathogenic In-Frame Deletions": cond8_no_inframe_del,
            "No Domain Overlap": cond5_no_domain, "Low Missense Count": cond6_missense,
            "Is <10% of Protein": cond4_small
        },
        "visualization": visualization_data
    }

def process_single_variant(query: str, client: EnsemblClient, splice_user_input: Optional[str] = None) -> Dict[str, Any]:
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
        final_result = {
            "summary": {"gene": gene_symbol, "transcript_id": definitive_transcript_id, **gene_characteristics},
            "assessments": {}
        }
        
        # --- 4. N1C Registry Check (Exit early if matched) ---
        n1c_result = check_n1c_registry(gene_symbol, query, hgvs_query)
        if n1c_result:
            final_result["assessments"]["N1C_Registry_Check"] = n1c_result
            return final_result

        # --- 5. Gene-Level Strategies (Knockdown, WT Upregulation) ---
        moa = set(gene_characteristics.get("moa", []))
        moi = set(gene_characteristics.get("moi", []))
        haplo_status_text = gene_characteristics.get("haploinsufficiency", {}).get("text", "Unknown")

        if "GoF" in moa:
            final_result["assessments"]["Allele_Specific_Knockdown"] = assess_knockdown(gene_characteristics)
        
        is_lof_ad = "LoF" in moa and any("Autosomal Dominant" in m for m in moi)
        if haplo_status_text == "Sufficient evidence" or is_lof_ad:
            final_result["assessments"]["WT_Upregulation"] = assess_wt_upregulation(client, gene_id, gene_symbol)

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
                details = {"Confirmation Method": "User-provided validation (qPCR/RNA-seq)"}
                splice_assessment = _evaluate_splice_variant_position(variant_identifier_from_vep, vep_entry, details)
            elif splice_user_input == 'no':
                splice_assessment = {"classification": "Not Eligible", "reason": "User confirmed no known splice-altering effect."}
            else:
                splice_assessment = assess_splice_switching(variant_identifier_from_vep, vep_entry, gene_symbol)
            
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
                    exon_skip_result = assess_single_exon(client, query, transcript_data, all_exons, target_exon, vep_entry, refseq_id_for_viewer)
                    if "visualization" in exon_skip_result and exon_skip_result["visualization"]:
                        final_result["visualization"] = exon_skip_result.pop("visualization")
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

@app.route('/')
def index(): return render_template('index.html', title="Tool")
@app.route('/about')
def about(): return render_template('about.html', title="About/Methods")
@app.route('/cite')
def cite(): return render_template('cite.html', title="How to Cite")

@app.route('/api_docs')
def api_docs():
    """Serves the API documentation page."""
    return render_template('api_docs.html', title="API Documentation")

@app.route('/api/v1/assess', methods=['GET'])
def api_assess():
    """
    Handles a single variant assessment via a GET request for programmatic access.
    Returns the full assessment data as JSON.
    """
    query = request.args.get('query')
    if not query:
        return jsonify({"error": "The 'query' parameter is required."}), 400

    client = EnsemblClient()
    result = process_single_variant(query, client)
    
    # Provide more specific HTTP status codes based on the outcome
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
    client = EnsemblClient()
    result = process_single_variant(query, client, splice_user_input=splice_input)
    
    return jsonify(result)
@app.route('/batch_assess', methods=['POST'])
def batch_assess():
    if 'file' not in request.files:
        return jsonify({"error": "No file part"}), 400
    file = request.files['file']
    if file.filename == '':
        return jsonify({"error": "No selected file"}), 400

    try:
        if file.filename.endswith('.xlsx'):
            df = pd.read_excel(file, header=None)
        else: # Handles .csv and .txt
            df = pd.read_csv(file, header=None)
    except Exception as e:
        return jsonify({"error": f"Error reading file: {e}"}), 400

    variants = df[0].dropna().astype(str).tolist()
    client = EnsemblClient()
    output_rows = []

    for variant in variants:
        result = process_single_variant(variant, client)
        
        row = {"Variant": variant}
        summary = result.get("summary", {})
        assessments = result.get("assessments", {})

        row["Gene"] = summary.get("gene", "N/A")
        row["MOI"] = ', '.join(summary.get("moi", []))
        row["MOA"] = ', '.join(summary.get("moa", []))
        
        haplo_info = summary.get("haploinsufficiency", {})
        row["Haploinsufficiency"] = haplo_info.get("text", "N/A")
        row["ClinGen Link"] = haplo_info.get("url", "N/A")
        
        # --- START: NEW COLUMN LOGIC ---

        # 1. Assess if an ASO exists and add the N1C link
        n1c_check = assessments.get("N1C_Registry_Check", {})
        if n1c_check:
            row["Existing ASO (N1C)"] = "Yes"
            row["N1C Registry Link"] = n1c_check.get("link", "N/A")
        else:
            row["Existing ASO (N1C)"] = "No"
            row["N1C Registry Link"] = "N/A"
            
        # 2. Get the Antisense Transcript ID for WT Upregulation
        wt_up = assessments.get("WT_Upregulation", {})
        antisense_ids = wt_up.get("antisense_gene_ids", [])
        row["Antisense Transcript ID"] = ", ".join(antisense_ids) if antisense_ids else "N/A"

        # --- END: NEW COLUMN LOGIC ---

        # Exon Skipping
        skip = assessments.get("Exon_Skipping", {})
        row["Exon Skipping Assessment"] = skip.get("classification", "NA")
        for check, status in skip.get("checks", {}).items():
            row[f"ES Check: {check}"] = status

        # Splice Correction
        splice = assessments.get("Splice_Switching", {})
        row["Splice Correction Assessment"] = splice.get("classification", "Unable to Assess")
        row["Splicing Validation DOI"] = splice.get("details", {}).get("Publication DOI", "NA")

        # WT Upregulation and Knockdown (assessments remain)
        row["WT-Upregulation"] = wt_up.get("classification", "NA")
        row["Knockdown"] = assessments.get("Allele_Specific_Knockdown", {}).get("classification", "NA")
        
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
       
if __name__ == '__main__':
    load_databases()

    app.run(debug=True)
