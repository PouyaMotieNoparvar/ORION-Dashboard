import streamlit as st
import json
import os
import pandas as pd
import altair as alt
from Bio.Seq import Seq

# Set page configuration
st.set_page_config(page_title="ORION Evolution Dashboard", layout="wide", initial_sidebar_state="expanded")

# --- Custom CSS for Styling ---
st.markdown("""
    <style>
    .main {
        background-color: #0e1117;
    }
    .stSidebar {
        background-color: #1a1c24;
    }
    h1 {
        color: #00d4ff;
        font-family: 'Inter', sans-serif;
    }
    </style>
""", unsafe_allow_html=True)

def load_data():
    viz_data_path = 'data/viz_data.json'
    if not os.path.exists(viz_data_path):
        st.error(f"Missing data: {viz_data_path}. Please run scripts/07_generate_viz_data.py first.")
        return None
    with open(viz_data_path, 'r') as f:
        data = json.load(f)
    
    # Convert list of MSA dicts to a single dict for easier access
    if isinstance(data.get('msa'), list):
        data['msa'] = {item['name']: item['sequence'] for item in data['msa']}
        
    return data

def load_disruption_data():
    path = 'data/non_coding_quantification.json'
    if os.path.exists(path):
        with open(path, 'r') as f:
            return json.load(f)
    return {}

def load_gff():
    # Use the fixed GFF with correct exon coordinates
    gff_path = "data/orion_msa_fixed.gff"
    if not os.path.exists(gff_path): 
        # Fallback to old one if fixed not found
        gff_path = "data/orion_msa_ugene_final.gff"
        if not os.path.exists(gff_path): return {}
    
    features = {}
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.strip().split("\t")
            if len(parts) < 9: continue
            seq_id, source, f_type, start, end, score, strand, phase, attr = parts
            if seq_id not in features: features[seq_id] = []
            features[seq_id].append({
                "type": f_type,
                "start": int(start) - 1,
                "end": int(end),
                "strand": strand
            })
    return features

def main():
    st.title("🧬 ORION Gene Evolution Dashboard")
    st.markdown("---")

    data = load_data()
    if not data: return
    disruption = load_disruption_data()
    gff_features = load_gff()

    # --- 1. Root Cause Resolution: Robust Metadata Initialization ---
    # Normalize MSA keys for robust matching
    msa_norm = {k.strip(): v for k, v in data['msa'].items()}
    alignment_len = len(next(iter(msa_norm.values())))
    
    # Robust Tree Traversal (Handles all node formats)
    species_order = []
    def traverse(node):
        if 'name' in node and node['name'] in msa_norm:
            if node['name'] not in species_order:
                species_order.append(node['name'])
        if 'children' in node:
            for child in node['children']:
                traverse(child)
    traverse(data['tree'])
    
    # Fill missing species
    for s in msa_norm:
        if s not in species_order:
            species_order.append(s)

    # --- 2. Strand Handling ---
    st.sidebar.subheader("DNA Strand Orientation")
    if 'strand' not in st.session_state:
        st.session_state.strand = "Forward"
    
    strand = st.sidebar.radio("Sense Strand", ["Forward", "Reverse Complement"], 
                             index=0 if st.session_state.strand == "Forward" else 1)
    
    if strand != st.session_state.strand:
        st.session_state.strand = strand
        st.rerun()

    # Pre-process MSA based on strand
    if strand == "Reverse Complement":
        msa_viz = {k: str(Seq(v).reverse_complement()) for k, v in msa_norm.items()}
    else:
        msa_viz = msa_norm

    # Locate Arabidopsis Key
    ath_key = next((k for k in msa_norm if 'thaliana' in k.lower()), None)

    # --- 3. Navigation ---
    st.sidebar.subheader("Navigation")
    if st.sidebar.button("🚀 Jump to Start Codon"):
        full_seq = msa_viz[ath_key].upper() if ath_key else ""
        clean_dna_str = full_seq.replace('-', '')
        ath_cds_seed = "ATGGCAACATTGGCACCAGTTCATGCTGTT"
        seed_idx = clean_dna_str.find(ath_cds_seed)
        
        if seed_idx != -1:
            curr_clean = 0
            for idx, char in enumerate(full_seq):
                if char != '-':
                    if curr_clean == seed_idx:
                        st.session_state.start_pos = max(0, idx - 10)
                        break
                    curr_clean += 1
        st.rerun()

    if 'start_pos' not in st.session_state:
        st.session_state.start_pos = 0

    window_size = st.sidebar.slider("Window Size (bp)", 20, 300, 100)
    start_pos = st.sidebar.slider("Genomic Position", 0, alignment_len - window_size, st.session_state.start_pos)
    st.session_state.start_pos = start_pos
    end_pos = start_pos + window_size

    show_fs = st.sidebar.toggle("Show Frameshifts", value=True)
    show_stops = st.sidebar.toggle("Show Stop Codons", value=True)
    
    # --- Species Filtering ---
    st.sidebar.subheader("Filter Species")
    categories = {
        "All Species": species_order,
        "Arabidopsis Only": [ath_key] if ath_key else [],
        "Coding Orthologs": [s for s in species_order if not disruption.get(s, {}).get('frameshift', True)],
        "Non-coding Precursors": [s for s in species_order if disruption.get(s, {}).get('frameshift', False) and 'aethionema' in s.lower()],
        "Outgroups": [s for s in species_order if 'aethionema' in s.lower() or 'schrenkiella' in s.lower()]
    }
    selected_cat = st.sidebar.multiselect("View Groups", list(categories.keys()), default=["All Species"])
    
    filtered_species = []
    if "All Species" in selected_cat:
        filtered_species = species_order
    else:
        for cat in selected_cat:
            filtered_species.extend(categories[cat])
    filtered_species = list(dict.fromkeys(filtered_species)) # unique while preserving order

    # --- Optional Analysis ---
    st.sidebar.subheader("Analysis Options")
    show_table = st.sidebar.checkbox("📊 Show Species Coding Potential", value=False)
    
    st.sidebar.markdown("---")

    # --- 4. Data Processing (Gap-Aware Translation) ---
    if 'msa_viz' not in locals():
        msa_viz = msa_norm # Fallback

    # --- 3. Ground-Truth Data Integration ---
    # Load official sequences
    # Load official sequences
    ath_protein = "MATLAPVHAVSVLHSRKRGGCNGVPLPYLIRPRRKIPPVDTFDFYNDSDDESGEGDKGHSDDEYVVVDKAEAIMPNVKTEKPGAIQMMELSKDAVGDDESKAKDAGYSSDEWFVV"
    ath_cds = "ATGGCAACATTGGCACCAGTTCATGCTGTTTCTGTTCTGCATTCTCGGAAGAGGGGAGGCTGCAACGGAGTTCCTTTACCATATCTGATCCGACCTCGTAGGAAAATTCCCCCTGTTGATACGTTTGATTTCTACAACGATTCTGATGATGAGAGTGGTGAAGGTGATAAAGGTCATTCAGATGACGAATATGTAGTCGTAGACAAAGCTGAAGCTATTATGCCCAATGTTAAAACAGAGAAACCAGGAGCTATTCAAATGATGGAACTCTCAAAAGATGCAGTCGGTGATGATGAGTCTAAAGCCAAAGATGCTGGATACAGCAGTGATGAGTGGTTTGTTGTGTAA"
    
    rows = []
    aa_rows = []
    codon_data = [] 
    ath_atg_indices = set()
    ath_stop_indices = set()
    first_codon_pos = 0

    if ath_key:
        full_seq = msa_viz[ath_key].upper()
        clean_dna = []
        clean_to_msa = []
        for idx, char in enumerate(full_seq):
            if char != '-':
                clean_dna.append(char)
                clean_to_msa.append(idx)
        clean_dna_str = "".join(clean_dna)
        
        # Search for Exons individually to handle introns
        # We use the first 30bp as a seed for the Start Codon
        exon_start_seed = ath_cds[:30]
        cds_start_in_clean = clean_dna_str.find(exon_start_seed)
        
        if cds_start_in_clean != -1:
            first_codon_pos = clean_to_msa[cds_start_in_clean]
            
            # For a simpler approach that handles introns, we map the CDS characters one by one
            # by looking for matches in the genomic sequence starting from the seed
            dna_map = []
            curr_pos = cds_start_in_clean
            for char in ath_cds:
                # Find this character in clean_dna_str starting from curr_pos
                # (Simple greedy approach for high-identity match)
                match_pos = clean_dna_str.find(char, curr_pos)
                if match_pos != -1:
                    dna_map.append(clean_to_msa[match_pos])
                    curr_pos = match_pos + 1
                else:
                    # If we hit a gap/mismatch, just skip
                    pass
            
            # Map protein to dna_map
            for c_id, aa in enumerate(ath_protein):
                c_idx = c_id * 3
                if c_idx + 2 < len(dna_map):
                    m1, m2, m3 = dna_map[c_idx], dna_map[c_idx+1], dna_map[c_idx+2]
                    if aa == 'M' and c_id == 0: ath_atg_indices.update([m1, m2, m3])
                    
                    # Color all 3, but only text for the middle one
                    codon_data.append({"Species": "Codon Frame", "Position": m1, "AA": aa, "DisplayAA": "", "CodonID": c_id})
                    codon_data.append({"Species": "Codon Frame", "Position": m2, "AA": aa, "DisplayAA": aa, "CodonID": c_id})
                    codon_data.append({"Species": "Codon Frame", "Position": m3, "AA": aa, "DisplayAA": "", "CodonID": c_id})
                    
                    mid_pos = m2
                    if mid_pos >= start_pos and mid_pos < end_pos:
                        aa_rows.append({"Position": mid_pos, "AA": aa})
            
            # Add Stop Codon
            c_idx = len(ath_protein) * 3
            if c_idx + 2 < len(dna_map):
                m1, m2, m3 = dna_map[c_idx], dna_map[c_idx+1], dna_map[c_idx+2]
                ath_stop_indices.update([m1, m2, m3])
                codon_data.append({"Species": "Codon Frame", "Position": m1, "AA": "*", "DisplayAA": "", "CodonID": 999})
                codon_data.append({"Species": "Codon Frame", "Position": m2, "AA": "*", "DisplayAA": "*", "CodonID": 999})
                codon_data.append({"Species": "Codon Frame", "Position": m3, "AA": "*", "DisplayAA": "", "CodonID": 999})

    # Update Jump to Start logic to use the found first_codon_pos
    # (I'll do this in a second edit to avoid target overlap)

    # Prepare Rows for species
    for species in filtered_species:
        seq = msa_viz[species].upper()
        sub_seq = seq[start_pos:end_pos]
        meta = disruption.get(species, {})
        
        for i, char in enumerate(sub_seq):
            abs_pos = start_pos + i
            feature = "Normal"
            if show_fs and meta.get('frameshift', False): feature = "Frameshift"
            
            if species == ath_key:
                if abs_pos in ath_atg_indices:
                    feature = "Start Codon"
                    if 'target_atg' in locals() and abs_pos >= target_atg and abs_pos < target_atg + 3:
                        feature = "Primary Start"
                elif abs_pos in ath_stop_indices:
                    feature = "Stop Codon"

            # Stroke Logic
            stroke_color = 'rgba(255,255,255,0.1)'
            stroke_width = 1
            if feature == "Primary Start": stroke_color = '#00d4ff'; stroke_width = 4
            elif feature == "Start Codon": stroke_color = '#60a5fa'; stroke_width = 2
            elif feature == "Stop Codon": stroke_color = 'red'; stroke_width = 2
            elif feature == "Frameshift": stroke_color = '#ff00ff'; stroke_width = 1

            rows.append({
                "Species": species.replace('_', ' ').title(),
                "Position": abs_pos,
                "Nucleotide": char,
                "Feature": feature,
                "StrokeColor": stroke_color,
                "StrokeWidth": stroke_width
            })
    
    # Add AA Sequence Row Data
    df_aa_row = pd.DataFrame([{"Species": "AA Sequence", "Position": d['Position'], "AA": d['AA']} for d in aa_rows])
    df_codon = pd.DataFrame([d for d in codon_data if d['Position'] >= start_pos and d['Position'] < end_pos])
    df = pd.DataFrame(rows)
    
    display_order = ["AA Sequence", "Codon Frame"] + [s.replace('_', ' ').title() for s in filtered_species]

    # --- 4. Chart Construction ---
    color_scale = alt.Scale(
        domain=['A', 'T', 'C', 'G', '-'],
        range=['#f87171', '#60a5fa', '#34d399', '#fbbf24', '#333333']
    )
    
    # Full 20 AA Domain for consistent coloring across all residues
    aa_list = list("ACDEFGHIKLMNPQRSTVWY*")
    aa_color_scale = alt.Scale(
        domain=aa_list,
        scheme='turbo'
    )

    # Custom CSS for table legibility and full-width charts
    st.markdown("""
        <style>
        /* Force table text to be white for readability in dark mode */
        [data-testid="stTable"] table {
            color: white !important;
        }
        [data-testid="stTable"] th, [data-testid="stTable"] td {
            color: white !important;
        }
        /* Ensure the main container uses all available space */
        .main .block-container {
            max-width: 95% !important;
            padding-top: 1rem;
        }
        </style>
    """, unsafe_allow_html=True)

    # --- TOP CHART: Unified Functional Frame ---
    codon_base = alt.Chart(df_codon).encode(
        x=alt.X('Position:O', axis=alt.Axis(grid=False), title=None, type='ordinal'),
        y=alt.Y('Species:N', sort=display_order[:2], axis=alt.Axis(grid=False), title=None)
    )

    codon_rects = codon_base.mark_rect(strokeWidth=0).encode(
        color=alt.Color('AA:N', scale=aa_color_scale, legend=None, type='nominal'),
        tooltip=[alt.Tooltip('AA:N'), alt.Tooltip('CodonID:Q')]
    ).properties(height=50)

    aa_text = codon_base.mark_text(fontWeight='bold', color='black').encode(
        text=alt.Text('DisplayAA:N')
    )

    top_chart = (codon_rects + aa_text).properties(
        width=1600,
        title="Functional Coding Frame (Arabidopsis Reference)"
    )

    # --- BOTTOM CHART: Species Alignment ---
    base = alt.Chart(df).encode(
        x=alt.X('Position:O', axis=alt.Axis(labelAngle=0, grid=False), title="MSA Position"),
        y=alt.Y('Species:N', sort=display_order[2:], title=None)
    )

    nt_rects = base.mark_rect().encode(
        color=alt.Color('Nucleotide:N', scale=color_scale, legend=None),
        stroke=alt.Stroke('StrokeColor:N', scale=None),
        strokeWidth=alt.StrokeWidth('StrokeWidth:Q', scale=None),
        tooltip=['Species', 'Position', 'Nucleotide', 'Feature']
    )

    nt_text = base.mark_text(baseline='middle').encode(
        text='Nucleotide:N',
        color=alt.condition(alt.datum.Nucleotide == '-', alt.value('transparent'), alt.value('white'))
    )

    bottom_chart = (nt_rects + nt_text).properties(
        width=1600,
        height=len(filtered_species) * 25
    )

    # Combine charts - Sync X-axis
    combined_chart = alt.vconcat(top_chart, bottom_chart, spacing=15).resolve_scale(x='shared')

    st.altair_chart(combined_chart, use_container_width=True)

    # --- 5. Optional Species Table ---
    if show_table:
        st.markdown("### 📊 Species Analysis (Current Window)")
        table_data = []
        for sp in filtered_species:
            # Calculate mean coding potential (percentage of non-gap characters in this window)
            sp_seq = msa_viz[sp].upper()
            win_seq = sp_seq[start_pos:end_pos]
            non_gaps = len([c for c in win_seq if c != '-'])
            potential = non_gaps / len(win_seq) if len(win_seq) > 0 else 0
            
            # Frameshift status
            meta = disruption.get(sp, {})
            has_fs = "⚠️ Yes" if meta.get('frameshift', False) else "✅ No"
            
            table_data.append({
                "Species": sp.replace('_', ' ').title(),
                "Mean Coding Potential": f"{potential:.2%}",
                "Frameshift": has_fs
            })
        
        st.table(pd.DataFrame(table_data))

    # --- Evolutionary Context Panel ---

if __name__ == "__main__":
    main()
