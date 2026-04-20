# Cross-Species Comparative AMR Genomics

**IS Element Ecology Shapes Divergent Carbapenem Resistance Gene Mobilisation Strategies Across Clinical Gram-Negative Pathogens: A Three-Species Comparative Genomic Analysis**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.placeholder.svg)](https://doi.org/10.5281/zenodo.placeholder)

---

## Overview

This repository contains the complete analysis pipeline and results for a comparative genomic study examining IS element architecture around carbapenem resistance genes in 926 Chinese clinical isolates across three WHO priority pathogens:

| Species | Genomes (QC) | Resistance Rate | Composite Transposon Rate | IS Burden AUC |
|---------|-------------|-----------------|--------------------------|---------------|
| *Klebsiella pneumoniae* | 492 | 40.7% | **94.0%** [83.5–98.8%] | 0.667 (p=4.2×10⁻¹⁰) |
| *Escherichia coli* | 266 | 10.5% | **62.1%** [44.0–77.3%] | 0.718\* (p=2.3×10⁻⁶) |
| *Acinetobacter baumannii* | 168 | 99.4% | **6.8%** | N/A (n=1 susceptible) |

\* PFAM HMM annotation-corrected value

**Key finding**: IS26-mediated composite transposon architecture is shared across Enterobacteriaceae (KPN 94%, ECO 62%) but near-absent in *A. baumannii* (6.8%), which instead uses ISAba1 upstream promoter insertion as the dominant carbapenemase activation mechanism. IS6 genomic copy number is a significant resistance predictor in both enterobacterial species.

---

## Repository Structure

```
comparative_amr_project/
├── analysis/
│   ├── config.py                    # Species configuration
│   ├── 01_download.py               # NCBI Entrez genome download
│   ├── 01b_download_eco_cr.py       # Targeted CR-ECO download
│   ├── 02_validate.py               # Species-aware QC
│   ├── 03_amr_scan.py               # AMR gene detection (GFF3 + CARD)
│   ├── 04_is_context.py             # IS flanking architecture
│   ├── 05_is_burden.py              # IS copy-number burden analysis
│   ├── 06_cross_species_figures.py  # Publication figures (Fig 1–4)
│   ├── 07_manuscript.py             # Manuscript draft generator
│   ├── 08_convert_to_pdf.py         # pandoc + Chrome PDF render
│   └── 09_is_hmmer_verify.py        # PFAM HMM IS annotation correction
├── figures/                         # fig1–4 (PDF + PNG, 300 dpi)
├── reports/
│   ├── manuscript_comparative_v3.md # Final manuscript (~5,400 words)
│   └── manuscript_comparative_v3.pdf
├── docs/
│   ├── OPERATION_MEMORY.md          # Full project state log
│   ├── WORKFLOW.md                  # Pipeline reference
│   └── SYSTEM_WHITEPAPER.md        # System architecture whitepaper
└── requirements.txt
```

---

## Pipeline

```bash
# 1. Download genomes (ECO + ABA; KPN reused from amr-genomics-analysis)
python analysis/01_download.py --species eco --limit 800
python analysis/01_download.py --species aba --limit 300

# 2. Quality control
python analysis/02_validate.py --species both

# 3. AMR gene detection
python analysis/03_amr_scan.py --species both

# 4. IS flanking architecture
python analysis/04_is_context.py --species both

# 5. IS burden statistics
python analysis/05_is_burden.py

# 5b. PFAM HMM annotation correction (required for ECO IS6)
python analysis/09_is_hmmer_verify.py --species both

# 6. Publication figures
python analysis/06_cross_species_figures.py

# 7. PDF manuscript
python analysis/08_convert_to_pdf.py --input reports/manuscript_comparative_v3.md
```

---

## Dependencies

```bash
pip install -r requirements.txt
```

- Python ≥ 3.11
- pandas, scipy, matplotlib, requests, tqdm, pyhmmer ≥ 0.12
- pandoc ≥ 3.9 (for PDF conversion)
- Google Chrome (headless PDF rendering)

*K. pneumoniae* data sourced from [amr-genomics-analysis](https://github.com/jiayu6954-sudo/amr-genomics-analysis).

---

## Key Scientific Findings

1. **Composite transposon gradient**: KPN(94%) > ECO(62.1%) > ABA(6.8%) — all pairwise comparisons significant
2. **IS6 burden predicts resistance** in both Enterobacteriaceae: KPN AUC=0.667, ECO AUC=0.718 (PFAM-corrected)
3. **ISAba1 upstream bias** (92.3% upstream of OXA-23) = promoter activation mechanism, NOT structural bracketing
4. **ECO IS annotation bias confirmed**: susceptible genomes have 9.9× more IS6-reclassified "unknown transposases" — PFAM correction is mandatory

---

## Citation

> ZJY. *IS Element Ecology Shapes Divergent Carbapenem Resistance Gene Mobilisation Strategies Across Clinical Gram-Negative Pathogens*. 2026. (manuscript in preparation)

---

## License

MIT License — see [LICENSE](LICENSE)

**Data**: All genome assemblies are publicly available from NCBI GenBank. Accession lists in `docs/OPERATION_MEMORY.md`.
