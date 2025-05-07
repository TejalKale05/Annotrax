import streamlit as st
import pandas as pd
import numpy as np
import base64
from Bio import Entrez, SeqIO
import io
import json
import time
import os

# Set your email address for NCBI Entrez API (mandatory)
Entrez.email = "kaletejal05@mail.com"  # Replace with your actual email address

# Page configuration
st.set_page_config(
    page_title="Annotrax",
    page_icon="🧬",
    layout="wide"
)

# App styling
st.markdown("""
<style>
    /* --- General Theme & Layout --- */
    body {
        font-family: 'Roboto', 'Arial', sans-serif; /* Modern, clean font */
        color: #333; /* Default text color */
        /* Subtle overall background gradient */
        background: linear-gradient(to bottom right, #e9f0f5, #f5f7fa); 
    }
    .main { 
        background-color: transparent; /* Make main container transparent to show body bg */
    }
    /* Overall container for content within tabs to have consistent padding */
    div[data-testid="stVerticalBlock"] div[data-testid="stVerticalBlock"] div[data-testid="stVerticalBlock"] > div[data-testid^="stHorizontalBlock"] {
         /* This targets the main content area inside a tab. May need adjustment if Streamlit structure changes. */
    }


    /* --- Headers --- */
    .main-header {
        font-size: 2.8rem; /* Slightly adjusted */
        color: #1A237E; /* Darker, more authoritative blue */
        text-align: center;
        padding-top: 25px;
        font-weight: 700; /* Bolder */
        letter-spacing: 0.5px;
        text-shadow: 1px 1px 2px rgba(0,0,0,0.05);
    }
    .tagline {
        font-size: 1.2rem; 
        color: #455A64; /* Muted grey-blue */
        text-align: center;
        margin-bottom: 35px; 
        font-style: italic;
        font-weight: 300;
    }
    .sub-header { /* Used for tab main headers like "Welcome", "Single Gene Search" */
        font-size: 2rem;
        color: #0D47A1; /* Strong blue for sub-headers */
        margin-top: 20px; /* Reduced top margin as tab provides some separation */
        margin-bottom: 20px;
        border-bottom: 3px solid #4FC3F7; /* Lighter blue accent line */
        padding-bottom: 10px;
        font-weight: 500;
    }
    .section-header { /* Used for "About Annotrax", "Download Options" */
        font-size: 1.6rem;
        color: #1976D2; /* Medium blue */
        margin-top: 30px;
        margin-bottom: 15px;
        font-weight: 500;
    }

    /* --- Result Box (Gene Search) --- */
    .result-box {
        background-color: #FFFFFF; 
        border-radius: 10px; /* Softer radius */
        padding: 25px 30px; /* More padding */
        margin-bottom: 30px;
        border: 1px solid #CFD8DC; /* Lighter border */
        box-shadow: 0 8px 16px rgba(0,0,0,0.06); /* Slightly more defined shadow */
        transition: box-shadow 0.3s ease-in-out;
    }
    .result-box:hover {
        box-shadow: 0 12px 24px rgba(0,0,0,0.08);
    }
    
    /* --- Text & Links --- */
    .info-text {
        color: #546E7A; 
        font-size: 0.9rem;
    }
    a, a:visited {
        color: #0277BD; /* Link color */
        text-decoration: none;
        transition: color 0.2s ease;
    }
    a:hover, a:active {
        color: #01579B; /* Darker on hover */
        text-decoration: underline;
    }
    .team-member-name {
        font-size: 1.3rem;
        font-weight: 600;
        color: #0D47A1; 
    }
    .team-member-description {
        font-size: 1rem;
        color: #37474F; 
        line-height: 1.6;
    }

    /* --- Tabs --- */
    .stTabs [data-baseweb="tab-list"] {
		gap: 4px; 
        border-bottom: 2px solid #B0BEC5; 
        margin-bottom: 20px; /* Space below tab bar */
	}
	.stTabs [data-baseweb="tab"] {
		height: 50px;
        white-space: pre-wrap;
		background-color: #ECEFF1; /* Lightest grey for inactive tabs */
		border-radius: 8px 8px 0px 0px; 
		padding: 10px 22px; 
        font-weight: 500;
        font-size: 1.05rem; /* Slightly larger tab font */
        color: #455A64; 
        transition: background-color 0.3s ease, color 0.3s ease, border-bottom-color 0.3s ease;
        border-bottom: 3px solid transparent; /* For smooth transition */
    }
    .stTabs [aria-selected="true"] {
  		background-color: #FFFFFF; /* White for active tab, stands out on gradient */
		color: #0D47A1 !important; /* Primary dark blue for active tab text */
        border-bottom: 3px solid #0288D1 !important; /* Bright blue accent on active tab */
        box-shadow: 0 -2px 5px rgba(0,0,0,0.05);
	}
    .stTabs [data-baseweb="tab"]:hover:not([aria-selected="true"]) { /* Hover for inactive tabs */
        background-color: #CFD8DC;
        color: #263238;
    }

    /* --- Sidebar --- */
    [data-testid="stSidebar"] {
        background-color: #fcfdff; /* Very light off-white, almost white */
        padding: 25px;
        border-right: 1px solid #E0E0E0;
        box-shadow: 2px 0px 10px rgba(0,0,0,0.03);
    }
    .stSidebar .sub-header { /* "Search Options" in sidebar */
        font-size: 1.5rem !important;
        color: #0D47A1;
        border-bottom: 2px solid #78909C;
        margin-top: 0px;
        padding-bottom: 8px;
    }
    .stSidebar .stRadio > label, .stSidebar .stSelectbox > label, .stSidebar .stCheckbox > label {
        font-weight: 500;
        color: #37474F;
        font-size: 0.95rem;
    }

    /* --- Buttons --- */
    div[data-testid="stButton"] > button {
        background-color: #0277BD; /* Primary Action Blue */
        color: white;
        border-radius: 6px;
        padding: 10px 20px; 
        border: none;
        font-weight: 500;
        font-size: 0.95rem;
        transition: background-color 0.2s ease-in-out, transform 0.1s ease, box-shadow 0.2s ease;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
    }
    div[data-testid="stButton"] > button:hover {
        background-color: #01579B; /* Darker blue on hover */
        transform: translateY(-2px); 
        box-shadow: 0 4px 8px rgba(0,0,0,0.15);
    }
    div[data-testid="stButton"] > button:active {
        transform: translateY(0px);
        box-shadow: 0 1px 2px rgba(0,0,0,0.1);
    }
    /* Download Button */
    div[data-testid="stDownloadButton"] > button {
        background-color: #4CAF50; /* Green for download - distinct */
        color: white; 
        border: none;
        /* Inherit other button styles like padding, radius etc. if needed */
        border-radius: 6px;
        padding: 10px 20px; 
        font-weight: 500;
        font-size: 0.95rem;
        transition: background-color 0.2s ease-in-out, transform 0.1s ease, box-shadow 0.2s ease;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
    }
    div[data-testid="stDownloadButton"] > button:hover {
        background-color: #388E3C; 
        transform: translateY(-2px); 
        box-shadow: 0 4px 8px rgba(0,0,0,0.15);
    }


    /* --- Input Fields --- */
    .stTextInput input, .stTextArea textarea {
        border: 1px solid #B0BEC5;
        border-radius: 6px;
        padding: 12px;
        background-color: #FAFAFC; /* Slightly off-white input */
        font-size: 0.95rem;
        transition: border-color 0.2s ease, box-shadow 0.2s ease;
    }
    .stTextInput input:focus, .stTextArea textarea:focus {
        border-color: #0288D1;
        background-color: #FFFFFF;
        box-shadow: 0 0 0 0.2rem rgba(2, 136, 209, 0.2); 
    }
    
    /* --- Expander --- */
    .stExpander {
        border: 1px solid #CFD8DC;
        border-radius: 8px;
        margin-bottom: 15px;
        background-color: #fff; /* Give expander a white background */
    }
    .stExpander > summary {
        font-size: 1.1rem;
        font-weight: 500;
        color: #1E88E5; /* Blue for expander header */
        padding: 10px 15px;
        background-color: #f7f9fa; /* Light header background */
        border-radius: 8px 8px 0 0; /* Match container if not expanded */
        border-bottom: 1px solid transparent; /* Prepare for border when open */
    }
    .stExpander[open] > summary {
        border-bottom: 1px solid #CFD8DC;
        border-radius: 8px 8px 0 0;
    }
    .stExpander div[data-testid="stVerticalBlock"] { /* Content inside expander */
        padding: 15px;
        background-color: #fff; /* Ensure content area is white */
        border-radius: 0 0 8px 8px;
    }

    /* --- Specific Link Section in Gene Search --- */
    .crosslink-section-header, .link-section-header { 
        font-size: 1.2rem; 
        font-weight: 600;
        color: #0D47A1; 
        margin-top:20px;
        margin-bottom:8px;
        padding-bottom: 5px;
        border-bottom: 1px dotted #90A4AE; /* Slightly darker dotted line */
    }
    .external-link-item {
        margin-bottom: 5px;
        font-size: 0.95rem;
        line-height: 1.5;
    }
    .external-link-item strong {
        color: #37474F; 
        margin-right: 5px;
    }

    /* --- Home Tab Specific Styling --- */
    .home-tab-content-wrapper { /* Wrapper for content inside home tab */
        background: transparent; 
        padding: 30px 40px; 
        /* min-height: calc(100vh - 200px); */ /* REMOVED to prevent forced height */
        display: flex;
        flex-direction: column;
    }
    .home-tab-content-wrapper .sub-header, 
    .home-tab-content-wrapper .section-header, 
    .home-tab-content-wrapper p,
    .home-tab-content-wrapper li,
    .home-tab-content-wrapper .team-member-name,
    .home-tab-content-wrapper .team-member-description {
        color: #263238; 
    }
    .home-tab-content-wrapper .sub-header { color: #0D47A1; margin-top: 0; } 
    .home-tab-content-wrapper .section-header { color: #1565C0; }
    .home-tab-content-wrapper .info-text { color: #455A64; background-color: transparent !important; }
    .home-tab-content-wrapper ul { padding-left: 25px; } /* Add some indent to lists */
    .home-tab-content-wrapper li { margin-bottom: 8px; } /* Space out list items */


    /* Footer styling */
    .footer-text {
        text-align: center; 
        padding: 25px 0;
        color: #546E7A;
        font-size: 0.85rem;
        border-top: 1px solid #E0E0E0;
        margin-top: 30px;
    }

</style>
""", unsafe_allow_html=True)

# Header
st.markdown("<h1 class='main-header'>🧬 Annotrax: Gene Annotation Explorer</h1>", unsafe_allow_html=True)
st.markdown("<p class='tagline'>Annotating Genes with Computational Excellence</p>", unsafe_allow_html=True)

# --- Helper Functions ---
def extract_id_from_url(url, db_name):
    """Extracts common DB IDs from their URLs. Simplified parser."""
    try:
        if db_name == "uniprot": return url.split('/')[-1].split('?')[0]
        elif db_name == "omim": return url.split('/')[-1].split('?')[0]
        elif db_name == "ensembl":
            if "?g=" in url: return url.split("?g=")[-1].split('&')[0]
            if "/Gene/Summary/" in url: return url.split("/Gene/Summary/")[-1].split("?")[0]
            return url.split('/')[-1].split('?')[0]
        elif db_name == "hgnc":
             # Handle potential trailing slashes or query params after ID
             if "HGNC:" in url: return url.split("HGNC:")[-1].split('/')[0].split('?')[0]
             # Fallback if no HGNC: prefix found (less reliable)
             return url.split('/')[-1].split('?')[0]
        elif db_name == "genecards":
            if "?gene=" in url: return url.split("?gene=")[-1].split('&')[0]
            return url.split('/')[-1].split('?')[0] # Add query param split as fallback
    except Exception: return None
    return None

def get_ucsc_link(chrom, position_str, organism):
    """Generates a UCSC Genome Browser link for a given position."""
    org_db_map = {
        "Homo sapiens": "hg38", "Mus musculus": "mm39", "Drosophila melanogaster": "dm6",
        "Caenorhabditis elegans": "ce11", "Saccharomyces cerevisiae": "sacCer3",
        "Danio rerio": "danRer11", "Arabidopsis thaliana": "araTha1"
        # Add other organisms your app supports and their UCSC codes
    }
    db = org_db_map.get(organism)
    if not db or not chrom or chrom == "N/A" or not position_str or position_str == "N/A" or "-" not in position_str:
        return None
    try:
        start, end = position_str.split('-')
        int(start); int(end) # Validate start/end look numeric
        
        chrom_str = str(chrom)
        # UCSC requires 'chr' prefix for standard chromosomes in many assemblies
        if chrom_str.isdigit() or chrom_str.upper() in ['X', 'Y', 'M', 'MT']:
             chrom_formatted = f"chr{chrom_str}" if not chrom_str.lower().startswith("chr") else chrom_str
        else:
             chrom_formatted = chrom_str # Use original for non-standard names
        
        # Ensure start/end are reasonable (optional basic check)
        if int(start) < 0 or int(end) < int(start):
             return None

        return f"https://genome.ucsc.edu/cgi-bin/hgTracks?db={db}&position={chrom_formatted}:{start}-{end}"
    except Exception: # Catch potential errors during splitting, int conversion, etc.
        return None

# --- NCBI Fetching Functions ---
@st.cache_data(ttl=3600) # Cache NCBI results for 1 hour
def fetch_gene_annotation(gene_name, organism="Homo sapiens"):
    """Fetches gene annotations and cross-links from NCBI."""
    st.session_state['ncbi_error_displayed'] = False # Reset error flag
    try:
        # Search for Gene ID
        time.sleep(0.34) # API courtesy
        search_term = f"{gene_name}[Gene Name] AND {organism}[Organism]"
        handle_search = Entrez.esearch(db="gene", term=search_term, retmax=5) 
        record_search = Entrez.read(handle_search)
        handle_search.close()
        if not record_search["IdList"]: return None
        
        results = []
        # Fetch summary for each ID found
        for gene_id in record_search["IdList"]:
            time.sleep(0.34)
            summary_handle = Entrez.esummary(db="gene", id=gene_id)
            summary = Entrez.read(summary_handle)
            summary_handle.close()
            
            if not summary or 'DocumentSummarySet' not in summary or not summary['DocumentSummarySet']['DocumentSummary']:
                continue 
            docsum = summary['DocumentSummarySet']['DocumentSummary'][0]

            # Initialize data structures
            cross_links = {"pdb": [], "uniprot": [], "omim_general": [], "ensembl": [], "hgnc": [], "genecards": []}
            omim_phenotypes = [] 

            # --- Parse Cross-Links from NCBI Summary ---
            for link_item in docsum.get("OtherLinks", []):
                link_name = link_item.get("Name", "")
                link_url = link_item.get("Url", "")
                link_category = link_item.get("Category","")
                link_subject_type = link_item.get("SubjectType", "") 

                if not link_url: continue

                omim_id_current = extract_id_from_url(link_url, "omim") # Extract potential OMIM ID once

                if "MIM" == link_name and ("Phenotypes" == link_category or "phenotype" in link_subject_type):
                    omim_desc = f"OMIM Entry {omim_id_current}" # Basic description
                    if omim_id_current and not any(d['id'] == omim_id_current for d in omim_phenotypes):
                        omim_phenotypes.append({"id": omim_id_current, "url": link_url, "description": omim_desc})
                elif "MIM" == link_name: # If it's OMIM but not clearly phenotype
                    if omim_id_current and not any(d['id'] == omim_id_current for d in cross_links["omim_general"]) and not any(d['id'] == omim_id_current for d in omim_phenotypes):
                         cross_links["omim_general"].append({"id": omim_id_current, "url": link_url})
                elif "UniProtKB/Swiss-Prot" in link_name or "UniProtKB/TrEMBL" in link_name:
                    uid = extract_id_from_url(link_url, "uniprot")
                    if uid and not any(d['id'] == uid for d in cross_links["uniprot"]): cross_links["uniprot"].append({"id": uid, "url": link_url})
                elif "Ensembl" == link_name:
                    ens_id = extract_id_from_url(link_url, "ensembl")
                    if ens_id and not any(d['id'] == ens_id for d in cross_links["ensembl"]): cross_links["ensembl"].append({"id": ens_id, "url": link_url})
                elif "HGNC" == link_name:
                    hgnc_id = extract_id_from_url(link_url, "hgnc")
                    # Check extracted ID exists and isn't already added with prefix
                    if hgnc_id and not any(d['id'] == f"HGNC:{hgnc_id}" for d in cross_links["hgnc"]): 
                        cross_links["hgnc"].append({"id": f"HGNC:{hgnc_id}", "url": link_url})
                elif "GeneCards" == link_name:
                    gc_id = extract_id_from_url(link_url, "genecards")
                    if gc_id and not any(d['id'] == gc_id for d in cross_links["genecards"]): cross_links["genecards"].append({"id": gc_id, "url": link_url})

            # --- Fetch PDB IDs via Protein links ---
            protein_ids_for_pdb = []
            try:
                time.sleep(0.34)
                handle_prot_link = Entrez.elink(dbfrom="gene", db="protein", id=gene_id, linkname="gene_protein_refseq", retmax=2)
                prot_link_record = Entrez.read(handle_prot_link)
                handle_prot_link.close()
                if prot_link_record and prot_link_record[0].get("LinkSetDb"):
                    for link in prot_link_record[0]["LinkSetDb"][0]["Link"]: protein_ids_for_pdb.append(link["Id"])
            except Exception: pass 

            if protein_ids_for_pdb:
                for prot_id in protein_ids_for_pdb[:1]: 
                    try:
                        time.sleep(0.34)
                        handle_pdb_link = Entrez.elink(dbfrom="protein", db="structure", id=prot_id, retmax=5) 
                        pdb_link_record = Entrez.read(handle_pdb_link)
                        handle_pdb_link.close()
                        if pdb_link_record and pdb_link_record[0].get("LinkSetDb"):
                            for link in pdb_link_record[0]["LinkSetDb"][0]["Link"]:
                                pdb_id = link["Id"]
                                if pdb_id and len(pdb_id) == 4 and not any(d['id'] == pdb_id for d in cross_links["pdb"]):
                                    cross_links["pdb"].append({"id": pdb_id, "url": f"https://www.rcsb.org/structure/{pdb_id}"})
                                if len(cross_links["pdb"]) >= 3: break 
                    except Exception: pass 
                    if len(cross_links["pdb"]) >= 3: break 
            
            # --- Fetch Nucleotide IDs ---
            nuccore_ids = []
            try:
                time.sleep(0.34)
                link_handle_refseq = Entrez.elink(dbfrom="gene", db="nuccore", id=gene_id, linkname="gene_nuccore_refseqrna", retmax=3)
                link_results_refseq = Entrez.read(link_handle_refseq)
                link_handle_refseq.close()
                if link_results_refseq and link_results_refseq[0].get("LinkSetDb"):
                     for link_db_item in link_results_refseq[0]["LinkSetDb"]:
                        for link in link_db_item["Link"]: 
                            if link["Id"] not in nuccore_ids: nuccore_ids.append(link["Id"]) # Avoid duplicates
                # Fallback if no RefSeq RNA found
                if not nuccore_ids:
                    time.sleep(0.34)
                    link_handle_general = Entrez.elink(dbfrom="gene", db="nuccore", id=gene_id, retmax=3)
                    link_results_general = Entrez.read(link_handle_general)
                    link_handle_general.close()
                    if link_results_general and link_results_general[0].get("LinkSetDb"):
                         for link_db_item in link_results_general[0]["LinkSetDb"]:
                            for link in link_db_item["Link"]:
                                if link["Id"] not in nuccore_ids: nuccore_ids.append(link["Id"])
                                if len(nuccore_ids) >=3: break 
                            if len(nuccore_ids) >=3: break
            except Exception: pass 
                
            # --- Parse Genomic Info ---
            genomic_info_list = docsum.get("GenomicInfo", [])
            chromosome, chr_start, chr_stop, exon_count = "N/A", "N/A", "N/A", "N/A"
            if genomic_info_list:
                chrom_info = next((info for info in genomic_info_list if info.get("ChrLoc")), genomic_info_list[0])
                chromosome = chrom_info.get("ChrLoc", "N/A")
                chr_start = chrom_info.get("ChrStart", "N/A")
                chr_stop = chrom_info.get("ChrStop", "N/A")
                exon_count = chrom_info.get("ExonCount", "N/A")
            
            # --- Assemble Final Gene Data ---
            gene_data = {
                "Gene ID": gene_id, 
                "Gene Symbol": docsum.get("Name", gene_name),
                "Official Name": docsum.get("Description", "N/A"),
                "Aliases": ", ".join(docsum.get("OtherAliases", "").split(", ")[:5]),
                "Function": docsum.get("Summary", "No function description available").strip(),
                "Chromosome": chromosome,
                "Position": f"{chr_start}-{chr_stop}" if chr_start != "N/A" and chr_stop != "N/A" else "N/A",
                "Exons": exon_count,
                "Source Organism": docsum.get("Organism", {}).get("ScientificName", organism),
                "Annotation Status": "Complete" if docsum.get("Status", "") == "live" and len(docsum.get("Summary", "")) > 50 else "Partial",
                "Nucleotide IDs": ", ".join(nuccore_ids[:3]) if nuccore_ids else "N/A", # Use processed list
                "CrossLinks": cross_links,
                "OMIM_Phenotypes": omim_phenotypes
            }
            results.append(gene_data)
            
        return results
    except Exception as e:
        st.error(f"Error during NCBI search/summary fetch: {e}")
        st.session_state['ncbi_error_displayed'] = True
        return None

@st.cache_data(ttl=3600)
def fetch_gene_sequence(nuccore_id):
    """Fetches sequence and basic features for a Nucleotide ID."""
    try:
        time.sleep(0.34) # API courtesy
        handle = Entrez.efetch(db="nucleotide", id=nuccore_id, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()
        
        features = []
        for feature in record.features:
            if feature.type in ["CDS", "mRNA", "gene", "exon", "ncRNA", "misc_feature", "regulatory"]:
                qualifiers = {k: str(v[0]) if isinstance(v, list) and v else str(v) 
                              for k, v in feature.qualifiers.items()}
                features.append({
                    "type": feature.type, 
                    "location": str(feature.location), 
                    "qualifiers": qualifiers
                })
        
        sequence_data = {
            "id": record.id, 
            "name": record.name, 
            "description": record.description, 
            "length": len(record.seq),
            "features": features[:20], 
            "sequence": str(record.seq)[:2000] + "..." if len(record.seq) > 2000 else str(record.seq)
        }
        return sequence_data
    except Exception as e:
        st.error(f"Error fetching sequence ({nuccore_id}): {e}")
        return None

def extract_feature_data(gene_data, feature_type):
    """Extracts specific annotation fields for download."""
    result = {}
    if feature_type == "function": result["Function"] = gene_data.get("Function", "")
    elif feature_type == "exons":
        result["Exons"] = gene_data.get("Exons", "")
        result["Position"] = gene_data.get("Position", "")
    elif feature_type == "organism": result["Source Organism"] = gene_data.get("Source Organism", "")
    return result

# --- Sidebar ---
st.sidebar.markdown("<h2 class='sub-header' style='font-size: 1.5rem; margin-bottom:0px;'>Search Options</h2>", unsafe_allow_html=True)
database_option = st.sidebar.radio("Filter by Annotation Status:", ["All Genes", "Fully Annotated Only", "Partially Annotated"], key="db_option_radio")
organism_options = ["Homo sapiens", "Mus musculus", "Drosophila melanogaster", "Caenorhabditis elegans", "Saccharomyces cerevisiae", "Danio rerio", "Arabidopsis thaliana", "Escherichia coli K-12"]
selected_organism = st.sidebar.selectbox("Select Organism:", organism_options, key="org_select")
search_mode = st.sidebar.radio("Search Mode:", ["Single Gene", "Multiple Genes (Batch)"], key="search_mode_radio")
st.sidebar.divider() # Visual separator
with st.sidebar.expander("About Annotrax (Quick View)", expanded=False):
    st.markdown("**Annotrax** simplifies access to gene annotations from NCBI. Use the main 'Home' tab for more details.")
st.sidebar.markdown("<div class='info-text' style='margin-top: 20px;'><small>Using NCBI Entrez. Please respect API usage guidelines. App developed by Tejal Kale.</small></div>", unsafe_allow_html=True)

# --- Main Page Tabs ---
tab_home, tab_search = st.tabs(["🏠 Home", "🔍 Gene Search"])

with tab_home:
    st.markdown("<div class='home-tab-content-wrapper'>", unsafe_allow_html=True) # Wrapper for home content
    st.markdown("<h2 class='sub-header'>Welcome to Annotrax!</h2>", unsafe_allow_html=True)
    st.markdown("""
    <p>Annotrax is a user-friendly web application meticulously designed to streamline the process of fetching and exploring gene annotations. 
    Whether you are a student embarking on your bioinformatics journey or a seasoned researcher in need of quick, reliable gene information, 
    Annotrax aims to be your indispensable assistant.</p>
    """, unsafe_allow_html=True)
    st.markdown("<h3 class='section-header'>About Annotrax</h3>", unsafe_allow_html=True)
    st.markdown("""
    <p><strong>Our Objective:</strong>
    To provide a simple, fast, and efficient gateway to search gene information from the comprehensive NCBI (National Center for Biotechnology Information)
    databases. Annotrax abstracts away the complexities of direct database navigation, presenting key gene annotations
    in an easily digestible format and facilitating convenient data download for offline analysis.</p>
    <p><strong>Key Features:</strong></p>
    <ul>
        <li><strong>Intuitive Gene Search</strong></li>
        <li><strong>Organism Specificity</strong></li>
        <li><strong>Annotation Filtering</strong></li>
        <li><strong>Comprehensive Details & Cross-Database Links (PDB, UniProt, OMIM, etc.)</strong></li>
        <li><strong>Associated Disease Phenotypes (from OMIM)</strong></li>
        <li><strong>Sequence Data Exploration with Features</strong></li>
        <li><strong>Efficient Batch Processing</strong></li>
        <li><strong>Customizable Downloads (CSV, JSON, TXT)</strong></li>
        <li><strong>Modern User Interface & Genome Browser Links</strong></li>
    </ul>
    <p><strong>Target Audience:</strong></p>
    <ul>
        <li>Biology and Bioinformatics students seeking to understand gene structures and functions.</li>
        <li>Researchers requiring quick and easy access to verified gene data.</li>
        <li>Educators looking for effective tools to demonstrate gene annotation concepts in classrooms and labs.</li>
    </ul>
    """, unsafe_allow_html=True)
    st.markdown("<h3 class='section-header'>The Developer</h3>", unsafe_allow_html=True)
    col1_team, col2_team = st.columns([0.8, 3]) 
    with col1_team:
        image_path = "Tejal.jpg" # CONFIRMED: Using Tejal.jpg
        if os.path.exists(image_path): st.image(image_path, width=130, caption="Tejal Kale") 
        else: 
            st.error(f"Error: Developer photo '{image_path}' not found.")
            st.markdown("<small><i>Please ensure the image file is in the same directory as the script.</i></small>", unsafe_allow_html=True)
    with col2_team:
        st.markdown("<p class='team-member-name'>Tejal Kale</p>", unsafe_allow_html=True)
        st.markdown("""
        <p class='team-member-description'>
        A passionate Bioinformatics Student at DES Pune University, Pune, and the developer behind Annotrax. 
        Driven by an interest in making complex biological data more accessible and actionable through computational tools.
        </p>
        """, unsafe_allow_html=True)
    st.markdown("</div>", unsafe_allow_html=True) # Close home-tab-content-wrapper div

# --- GENE SEARCH TAB ---
with tab_search:
    # Initialize session state for search status tracking
    if 'search_performed' not in st.session_state: st.session_state.search_performed = False 
    if 'last_search_input' not in st.session_state: st.session_state.last_search_input = ""
    
    st.session_state['ncbi_error_displayed'] = False # Reset error flag each time tab is viewed
    
    if search_mode == "Single Gene":
        st.markdown("<h2 class='sub-header'>Single Gene Search</h2>", unsafe_allow_html=True)
        
        with st.form(key='single_gene_search_form'):
            gene_name_input = st.text_input(
                "🔍 Enter gene symbol or name:", 
                key="single_gene_input_field",
                value=st.session_state.last_search_input, 
                placeholder="e.g., TP53, BRCA1, insr", 
                help="Enter a standard gene symbol (like TP53) or official name."
            )
            search_button = st.form_submit_button(label="Search Gene")

        # Determine if a search should be triggered
        trigger_search = search_button 

        if trigger_search and gene_name_input: 
            st.session_state.search_performed = True
            st.session_state.last_search_input = gene_name_input 
            with st.spinner(f"Searching for '{gene_name_input}' in {selected_organism}... Fetching details."):
                gene_annotations_list = fetch_gene_annotation(gene_name_input, selected_organism)
        elif not gene_name_input and 'search_performed' in st.session_state and st.session_state.search_performed:
             st.session_state.search_performed = False 
             st.session_state.last_search_input = ""
             gene_annotations_list = None 
        elif 'search_performed' in st.session_state and st.session_state.search_performed:
             with st.spinner(f"Retrieving cached data for '{st.session_state.last_search_input}'..."): 
                gene_annotations_list = fetch_gene_annotation(st.session_state.last_search_input, selected_organism)
        else:
            gene_annotations_list = None

        # Display logic
        if not st.session_state.search_performed:
             st.info("Enter a gene name above and click 'Search Gene'.")
        elif st.session_state.search_performed and gene_name_input: 
            if gene_annotations_list:
                # Filter results based on sidebar option
                if database_option == "Fully Annotated Only": display_list = [g for g in gene_annotations_list if g["Annotation Status"] == "Complete"]
                elif database_option == "Partially Annotated": display_list = [g for g in gene_annotations_list if g["Annotation Status"] == "Partial"]
                else: display_list = gene_annotations_list
                    
                if not display_list:
                    st.warning(f"⚠️ No genes matching '{gene_name_input}' with '{database_option}' criteria were found in {selected_organism}.")
                else:
                    st.success(f"✅ Displaying {len(display_list)} result(s) for '{gene_name_input}' in {selected_organism} matching '{database_option}'.")
                    
                    # --- Display Function ---
                    def display_gene_data_with_crosslinks(gene_data_item):
                        st.markdown("<div class='result-box'>", unsafe_allow_html=True)
                        c1_info, c2_info = st.columns([2.5, 1]) 
                        with c1_info:
                            st.markdown(f"#### {gene_data_item['Gene Symbol']} – {gene_data_item['Official Name']}")
                            st.markdown(f"**Aliases:** {gene_data_item.get('Aliases', 'N/A')}")
                            chrom = gene_data_item.get('Chromosome', 'N/A')
                            position_str = gene_data_item.get('Position', 'N/A')
                            organism = gene_data_item.get('Source Organism', 'N/A')
                            ucsc_url = get_ucsc_link(chrom, position_str, organism) 
                            location_text = f"**Location:** Chr {chrom}, Position {position_str}"
                            if ucsc_url: location_text += f" [<a href='{ucsc_url}' target='_blank' title='View on UCSC Genome Browser'>🧬</a>]"
                            st.markdown(location_text, unsafe_allow_html=True)
                            st.markdown(f"**Exons:** {gene_data_item.get('Exons', 'N/A')}")
                        with c2_info:
                            st.markdown(f"**Organism:** {organism}") 
                            st.markdown(f"**Gene ID:** [{gene_data_item['Gene ID']}](https://www.ncbi.nlm.nih.gov/gene/{gene_data_item['Gene ID']})")
                            st.markdown(f"**Annotation:** {gene_data_item['Annotation Status']}")
                        st.markdown(f"**Function Summary:**")
                        st.markdown(f"<small style='display:block; max-height: 100px; overflow-y:auto; padding:5px; border: 1px solid #eee; border-radius:4px; background-color:#fdfdfd;'>{gene_data_item.get('Function', 'No function description available.')}</small>", unsafe_allow_html=True)

                        omim_phenotypes_data = gene_data_item.get("OMIM_Phenotypes", [])
                        if omim_phenotypes_data:
                            st.markdown("<div class='link-section-header'>Associated Disease Phenotypes (OMIM):</div>", unsafe_allow_html=True)
                            for pheno in omim_phenotypes_data[:3]: st.markdown(f"<div class='external-link-item'>• <a href='{pheno['url']}' target='_blank'>OMIM:{pheno['id']}</a> - {pheno.get('description','View details')}</div>", unsafe_allow_html=True)
                        
                        cross_links_data = gene_data_item.get("CrossLinks", {})
                        other_db_links = {k: v for k, v in cross_links_data.items() if k != "omim_general" and v} 
                        omim_general_links = cross_links_data.get("omim_general", [])
                        if any(other_db_links.values()) or omim_general_links: 
                            st.markdown("<div class='crosslink-section-header'>Access to Protein Data Bank & Other Databases:</div>", unsafe_allow_html=True)
                            link_cols = st.columns(3)
                            link_col_map = { "PDB": other_db_links.get("pdb", []), "UniProt": other_db_links.get("uniprot", []), 
                                             "OMIM (General)": omim_general_links, "Ensembl": other_db_links.get("ensembl", []), 
                                             "HGNC": other_db_links.get("hgnc", []), "GeneCards": other_db_links.get("genecards", []) }
                            current_col_idx = 0
                            for db_name, links in link_col_map.items():
                                if links:
                                    with link_cols[current_col_idx % 3]:
                                        link_html = f"<strong>{db_name}:</strong> " + ", ".join([f"<a href='{link_info['url']}' target='_blank'>{link_info['id']}</a>" for link_info in links[:3]])
                                        st.markdown(f"<div class='external-link-item'>{link_html}</div>", unsafe_allow_html=True)
                                    current_col_idx +=1
                        
                        st.divider() 

                        if gene_data_item.get('Nucleotide IDs', 'N/A') != "N/A":
                            first_nuccore_id_val = gene_data_item['Nucleotide IDs'].split(", ")[0]
                            seq_button_key = f"seq_btn_item_{gene_data_item['Gene ID']}_{first_nuccore_id_val}"
                            if st.button(f"🧬 View Sequence Data for {first_nuccore_id_val}", key=seq_button_key):
                                with st.spinner(f"Fetching sequence {first_nuccore_id_val}..."): 
                                    seq_data_item = fetch_gene_sequence(first_nuccore_id_val)
                                if seq_data_item:
                                    with st.expander(f"Sequence Details for {seq_data_item['id']}", expanded=False):
                                        st.info(f"Description: {seq_data_item['description']}")
                                        st.markdown(f"**Length:** {seq_data_item['length']} bp")
                                        st.markdown("**Sequence (first 2000 bp):**")
                                        st.code(seq_data_item['sequence'], language=None)
                                        st.markdown("**Features (Top 15):**")
                                        if seq_data_item['features']:
                                            for feat in seq_data_item['features']:
                                                st.markdown(f"**{feat['type']}** at `{feat['location']}`")
                                                quals = []
                                                for k_qual, v_qual in feat['qualifiers'].items():
                                                    if k_qual in ['product', 'gene', 'note', 'protein_id', 'transcript_id', 'db_xref', 'standard_name']: 
                                                        quals.append(f"_{k_qual}_: {v_qual}")
                                                if quals: st.markdown(f"  > {'; '.join(quals)}")
                                        else: st.write("No relevant features extracted.")
                                else:
                                     st.error(f"Could not fetch sequence data for {first_nuccore_id_val}.") # Show error if fetch fails
                        st.markdown("</div>", unsafe_allow_html=True) # Close result-box

                    # --- Display Trigger ---
                    if len(display_list) > 1:
                        # Ensure unique keys for tabs if multiple searches happen without full reload
                        tab_keys = [f"tab_{g['Gene ID']}_{st.session_state.last_search_input}" for g in display_list] 
                        result_display_tabs = st.tabs([f"{g['Gene Symbol']} (ID: {g['Gene ID']})" for g in display_list])
                        for i, tab_content_item in enumerate(result_display_tabs):
                            with tab_content_item: display_gene_data_with_crosslinks(display_list[i])
                    elif len(display_list) == 1: 
                        display_gene_data_with_crosslinks(display_list[0])
                    
                    # --- Download Options ---
                    if display_list: 
                        st.divider() 
                        st.markdown("<h3 class='section-header'>📥 Download Annotation Data</h3>", unsafe_allow_html=True)
                        genes_to_prepare_for_dl = display_list 
                        col_dl_f, col_dl_fmt = st.columns(2)
                        # Add unique keys based on search term to avoid widget state issues on re-search
                        search_key_suffix = f"_{st.session_state.last_search_input}" 
                        with col_dl_f:
                            st.markdown("##### Select fields:")
                            dl_all_key = f"s_dl_all_fields{search_key_suffix}" 
                            dl_all_fields = st.checkbox("Download all available fields (excluding links/phenotypes)", value=True, key=dl_all_key, help="Downloads core annotations like ID, names, function, location, exons.")
                            dl_func, dl_ex, dl_org = False, False, False
                            if not dl_all_fields:
                                dl_func = st.checkbox("Gene function", value=True, key=f"s_dl_func{search_key_suffix}")
                                dl_ex = st.checkbox("Exon info & Position", value=False, key=f"s_dl_ex{search_key_suffix}")
                                dl_org = st.checkbox("Source organism", value=False, key=f"s_dl_org{search_key_suffix}")
                        with col_dl_fmt:
                            st.markdown("##### Select format:")
                            dl_format = st.radio("", ["CSV", "JSON", "TXT"], key=f"s_dl_format_radio{search_key_suffix}", horizontal=True)
                        
                        if st.button("📄 Download Data", key=f"s_dl_btn{search_key_suffix}"):
                            export_list = []
                            for g_data in genes_to_prepare_for_dl:
                                data_for_export = {k: v for k, v in g_data.items() if k not in ['CrossLinks', 'OMIM_Phenotypes']}
                                if dl_all_fields: export_list.append(data_for_export)
                                else:
                                    item_export = {"Gene Symbol": g_data.get("Gene Symbol"), "Gene ID": g_data.get("Gene ID")}
                                    if dl_func: item_export.update(extract_feature_data(g_data, "function"))
                                    if dl_ex: item_export.update(extract_feature_data(g_data, "exons"))
                                    if dl_org: item_export.update(extract_feature_data(g_data, "organism"))
                                    export_list.append(item_export)
                            df_export = pd.DataFrame(export_list)
                            dl_fname_prefix = st.session_state.last_search_input.replace(" ", "_").replace("/","_") if len(genes_to_prepare_for_dl) == 1 else "gene_selection" # Sanitize filename
                            if not df_export.empty:
                                dl_btn_key_base = f"dl_btn_{dl_format.lower()}{search_key_suffix}"
                                if dl_format == "CSV": st.download_button("⬇️ Download CSV", df_export.to_csv(index=False).encode("utf-8"), f"{dl_fname_prefix}_annotations.csv", "text/csv", key=f"{dl_btn_key_base}_csv")
                                elif dl_format == "JSON": st.download_button("⬇️ Download JSON", df_export.to_json(orient="records", indent=2).encode("utf-8"), f"{dl_fname_prefix}_annotations.json", "application/json", key=f"{dl_btn_key_base}_json")
                                else: st.download_button("⬇️ Download TXT", df_export.to_string(index=False).encode("utf-8"), f"{dl_fname_prefix}_annotations.txt", "text/plain", key=f"{dl_btn_key_base}_txt")
                            else: st.error("No data selected or available for download.")
            
            # Handle case where search triggered but fetch returned None
            elif trigger_search and not gene_annotations_list and not st.session_state.get('ncbi_error_displayed', False):
                 st.warning(f"⚠️ No gene data found for '{gene_name_input}' in {selected_organism}. Please check spelling or try another gene/organism combination.")


    else: # Batch search mode
        st.markdown("<h2 class='sub-header'>Batch Gene Search</h2>", unsafe_allow_html=True)
        st.markdown("Enter multiple gene symbols separated by commas, spaces, or new lines. (Max 50 genes per batch for API courtesy)")
        batch_gene_input = st.text_area("Enter gene names:", height=150, key="batch_gene_input_area", placeholder="e.g.\nTP53\nBRCA1, INS\nMYC EGFR")
        
        submit_batch = st.button("🚀 Search Batch", key="batch_search_exe_button")

        # Store results in session state to persist across reruns if GO button is clicked etc.
        if 'batch_results' not in st.session_state: st.session_state.batch_results = None
        if 'batch_gene_ids' not in st.session_state: st.session_state.batch_gene_ids = []

        if submit_batch and batch_gene_input: 
            raw_gene_names = [name.strip() for name in batch_gene_input.replace(",", " ").replace("\n", " ").split()]
            unique_gene_names = list(dict.fromkeys([name for name in raw_gene_names if name])) 
            if unique_gene_names:
                st.markdown(f"Attempting to search for **{len(unique_gene_names)}** unique gene names in **{selected_organism}**.")
                
                # Reset previous results before new search
                st.session_state.batch_results = [] 
                st.session_state.batch_gene_ids = []

                prog_bar, stat_text, MAX_BATCH_GENES = st.progress(0), st.empty(), 50
                genes_to_process = unique_gene_names[:MAX_BATCH_GENES]
                if len(unique_gene_names) > MAX_BATCH_GENES: st.warning(f"Processing the first {MAX_BATCH_GENES} genes due to API limits.")
                
                found_count = 0
                processed_count = 0
                
                progress_container = st.container()
                with progress_container:
                    prog_bar = st.progress(0)
                    stat_text = st.empty()

                for idx, g_name in enumerate(genes_to_process):
                    processed_count += 1
                    stat_text.text(f"Fetching data for '{g_name}'... ({processed_count}/{len(genes_to_process)})")
                    annotations = fetch_gene_annotation(g_name, selected_organism) 
                    
                    if annotations:
                        filtered_for_criteria = []
                        if database_option == "Fully Annotated Only": filtered_for_criteria = [g for g in annotations if g["Annotation Status"] == "Complete"]
                        elif database_option == "Partially Annotated": filtered_for_criteria = [g for g in annotations if g["Annotation Status"] == "Partial"]
                        else: filtered_for_criteria = annotations
                        
                        if filtered_for_criteria:
                            st.session_state.batch_results.extend(filtered_for_criteria)
                            found_count += len(filtered_for_criteria)
                            for gene_item in filtered_for_criteria: 
                                st.session_state.batch_gene_ids.append(str(gene_item["Gene ID"]))
                    
                    prog_bar.progress((idx + 1) / len(genes_to_process))
                
                # Final status update outside the loop
                if st.session_state.batch_results:
                    stat_text.success(f"Batch search complete! Found {found_count} gene record(s) matching criteria.")
                elif genes_to_process: 
                     stat_text.warning(f"⚠️ Search complete. No genes found matching your batch list and the filter criteria ('{database_option}') for {selected_organism}.")
            
            elif not unique_gene_names and batch_gene_input: 
                st.warning("Please enter valid gene names for batch search.")
            elif not unique_gene_names and not batch_gene_input:
                 st.warning("Please enter gene names in the text area before searching.")
        
        # Display results if they exist in session state
        if st.session_state.batch_results:
            df_batch_display_data = [{"Symbol": r.get("Gene Symbol"), "ID": r.get("Gene ID"), "Name": r.get("Official Name"), "Organism": r.get("Source Organism"), "Chr": r.get("Chromosome"), "Exons": r.get("Exons"), "Status": r.get("Annotation Status"), "PDBs": len(r.get("CrossLinks",{}).get("pdb",[])), "UniProt": len(r.get("CrossLinks",{}).get("uniprot",[])), "OMIM Pheno.": len(r.get("OMIM_Phenotypes",[]))} for r in st.session_state.batch_results]
            st.dataframe(pd.DataFrame(df_batch_display_data), height=300, use_container_width=True) 

            # --- Placeholder for GO Enrichment ---
            unique_go_ids = list(set(st.session_state.batch_gene_ids))
            if unique_go_ids and selected_organism == "Homo sapiens": 
                st.divider()
                st.markdown("<h3 class='section-header'>📊 Gene Ontology Enrichment</h3>", unsafe_allow_html=True)
                if st.button("Run GO Enrichment Analysis (Placeholder)", key="go_enrich_btn"):
                    st.info("GO Enrichment feature is under development.")
                    st.write(f"Genes for GO analysis (NCBI IDs): {', '.join(unique_go_ids)}") 
                    pass # Placeholder for actual analysis
                st.divider()
            
            # --- Batch Download Options ---
            st.markdown("<h3 class='section-header'>📥 Batch Download Options</h3>", unsafe_allow_html=True)
            col_b_dl_f, col_b_dl_fmt = st.columns(2)
            with col_b_dl_f:
                st.markdown("##### Select fields:")
                b_dl_all = st.checkbox("Download all fields (excluding links/phenotypes)", value=True, key="b_dl_all_f")
                b_dl_func, b_dl_ex, b_dl_org = False, False, False
                if not b_dl_all:
                    b_dl_func, b_dl_ex, b_dl_org = st.checkbox("Function", value=True, key="b_dl_fn"), st.checkbox("Exon/Position", value=False, key="b_dl_exn"), st.checkbox("Organism", value=False, key="b_dl_orgn")
            with col_b_dl_fmt:
                st.markdown("##### Select format:")
                b_dl_fmt = st.radio("", ["CSV", "JSON", "TXT"], key="b_dl_fmt_rad", horizontal=True)
            if st.button("📄 Download Batch Results", key="b_dl_res_btn"):
                b_export_list = []
                for g_res in st.session_state.batch_results: # Use results from session state
                    data_for_export_batch = {k: v for k, v in g_res.items() if k not in ['CrossLinks', 'OMIM_Phenotypes']}
                    if b_dl_all: b_export_list.append(data_for_export_batch)
                    else:
                        b_item = {"Gene Symbol": g_res.get("Gene Symbol"), "Gene ID": g_res.get("Gene ID")}
                        if b_dl_func: b_item.update(extract_feature_data(g_res, "function"))
                        if b_dl_ex: b_item.update(extract_feature_data(g_res, "exons"))
                        if b_dl_org: b_item.update(extract_feature_data(g_res, "organism"))
                        b_export_list.append(b_item)
                df_b_export = pd.DataFrame(b_export_list)
                if not df_b_export.empty:
                    if b_dl_fmt == "CSV": st.download_button("⬇️ Download CSV", df_b_export.to_csv(index=False).encode("utf-8"), "batch_annotations.csv", "text/csv", key="dl_batch_csv")
                    elif b_dl_fmt == "JSON": st.download_button("⬇️ Download JSON", df_b_export.to_json(orient="records", indent=2).encode("utf-8"), "batch_annotations.json", "application/json", key="dl_batch_json")
                    else: st.download_button("⬇️ Download TXT", df_b_export.to_string(index=False).encode("utf-8"), "batch_annotations.txt", "text/plain", key="dl_batch_txt")
                else: st.error("No data to download based on current selections.")
        
        # Initial state message for batch mode
        elif not submit_batch and not batch_gene_input:
             st.info("Enter one or more gene names/symbols separated by spaces, commas, or new lines, then click 'Search Batch'.")


# Footer
st.divider()
st.markdown("<div class='footer-text'><small>Annotrax retrieves data from NCBI's Entrez API. Please use responsibly and respect NCBI's usage guidelines.<br>This tool is for educational and research purposes. Data accuracy depends on NCBI records. Developed by Tejal Kale.</small></div>", unsafe_allow_html=True)