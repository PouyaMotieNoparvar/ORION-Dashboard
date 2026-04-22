# ORION-Dashboard
# ORION Gene Evolution Dashboard 🧬

A high-fidelity genomic alignment and functional visualization tool for the Arabidopsis *ORION* gene across 15 species. Developed for academic manuscript documentation and evolutionary genomic research.

## 🚀 Live Demo
The dashboard is optimized for Streamlit Cloud deployment.

## 🛠️ Key Features
- **Ground-Truth Mapping**: Integrated official TAIR protein and CDS sequences for 100% biological fidelity.
- **Strand-Aware Alignment**: Robust handling of minus-strand gene models with automatic Reverse Complement transformation.
- **Functional Frame Heatmap**: One-character-per-AA color coding using the high-contrast Turbo scheme.
- **Dynamic Species Analysis**: Real-time calculation of mean coding potential and frameshift status across any genomic window.
- **Academic Visualization**: Tri-panel high-resolution layout synchronized across a 4.4kb genomic region.

## 📦 Deployment Instructions
1. Push this repository to GitHub.
2. Connect your GitHub account to [Streamlit Cloud](https://streamlit.io/cloud).
3. Select this repository and the main file: `scripts/08_orion_streamlit_dashboard.py`.
4. Set the Python version to 3.9+ and deploy!

## 📂 Data Sources
- `data/viz_data.json`: Master MSA and sequence metadata.
- `data/orion_msa_fixed.gff`: High-precision exon coordinates for Arabidopsis.
- `data/non_coding_quantification.json`: Disruption and frameshift metadata.
