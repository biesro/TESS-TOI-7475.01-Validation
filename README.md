# Statistical Validation and Photometric Characterization of TOI 7475.01

⚠️ IMPORTANT: This repository is old and may be inaccurate because of my low level at the moment, for a real corrected and robust pipeline go see the https://github.com/biesro/TESS-TOI-7701.01-Validation

**Code repository for the paper:**
> *Statistical Validation and Photometric Characterization of the Hot Jupiter Candidate TOI 7475.01*

![Python Version](https://img.shields.io/badge/python-3.10-blue)
![License](https://img.shields.io/badge/license-MIT-green)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18641713.svg)](https://doi.org/10.5281/zenodo.18641713)
[![Status](https://img.shields.io/badge/Status-Validated%20Planet-success.svg)]()
[![arXiv](https://img.shields.io/badge/arXiv-2602.14840-b31b1b.svg)](https://arxiv.org/abs/2602.14840)

**Principal Investigator:** Biel Escolà-Rodrigo  
**Date:** March 2026 

---

## 🔭 Abstract

This repository contains the complete analysis pipeline used to statistically validate and characterize the exoplanet **TOI 7475.01** (TIC 376866659).

Using a custom pipeline and the BLS algorithm, we identified a transit signal with a period of **3.2538 days** and a depth of **~4600 ppm**. Our analysis confirms a clean spatial environment and a False Positive Probability (FPP) of **~0**, validating the signal as a planetary companion. We further performed a Bayesian transit fit to derive precise physical parameters.

Key findings supporting this validation include:

* **High Signal-to-Noise:** A detection SNR of **294.13**.
* **Clean Environment:** Gaia DR3 analysis shows no contaminating sources within the critical aperture (nearest neighbor at 28.3", $\Delta G \approx 6$).
* **Stable Centroids:** No photocenter shift during transit, confirming the signal originates from the target star.
* **Statistical Validation:** A rigorous `triceratops` analysis yields an FPP of **0.000000**.

## 🌍 Key Findings & Parameters

Based on the analysis of **TESS Sector 91**:

| Parameter | Value | Note |
| :--- | :--- | :--- |
| **Target** | TIC 376866659 |
| **Signal Period** | **3.253773 days** | Ultra-stable periodicity |
| **Transit Depth** | **4601 ppm** | Consistent U-shape |
| **Planet Radius** | **1.18 R_Jup** | Jupiter-sized body |
| **Equilibrium Temp** | **1455 K** | Hot Jupiter classification |
| **Planet Mass (MAP)** | **~2.2 M_Jup** | Estimated via Chen & Kipping |
| **Validation Status** | **FPP ~ 0** | Statistically Validated |

## ⚙️ Methodology & Software Stack
The analysis follows a rigorous forensic protocol implemented in Python.

### 1. Detection & Signal Recovery 
* **Tool:** `Lightkurve` (v2.x) + Custom Pipeline
* **Finding:** Recovers a periodic signal with a distinct U-shaped transit and consistent odd/even depths.
* **Script:** `01_detection_BLS.py`

### 2. Spatial Vetting & Centroid Analysis 
* **Tool:** `astroquery` (Gaia DR3) + `Lightkurve` (Moments Method) 
* **Finding:** The field is clean with a single-star RUWE of 1.02. Centroids remain stable during the transit window.
* **Script:** `02_centroid_FINAL.py`

### 3. Statistical Validation (TRICERATOPS) 
* **Tool:** `triceratops` (v1.0.18) 
* **Finding:** Mean FPP of **0.000000** across 20 runs, ruling out eclipsing binaries.
* **Script:** `03_triceratops_vetting.py`

### 4. Bayesian Characterization 
* **Tool:** `juliet` + `Dynesty` 
* **Finding:** Derived physical parameters via Monte Carlo error propagation. Note: The impact parameter (b = 0.46 +-0.35) remains unconstrained.
* **Script:** `04_characterization.py`

## 📂 Repository Structure

```text
TOI-7475.01-Validation/
│
├── LICENSE                     # MIT License
├── README.md                   # Project documentation
├── requirements.txt            # Python dependencies
├── TOI7475.01_Folded_Raw_GrayPoints.csv
├── characterization_summary.txt
│
├── code/                    # Analysis scripts
│   ├── 01_detection_BLS.py
│   ├── 01.2_spatial_contamination_check.py
│   ├── 02_centroid_FINAL.py
│   ├── 03_triceratops_vetting.ipynb
│   └── 04_characterization.py
│
└── figures/                    # Generated plots
    ├── Figure_01_BLS+Transit+OddEven.png
    ├── Figure_02_spatial_contam_check.png
    ├── Figure_03_Centroid.png
    ├── TOI7475_triceratops_20runs_histograms.png
    ├── corner_characterization.png
    ├── diagnostic_full_lc.png
    ├── phased_check.png
    ├── population_comparison.png
    └── posterior_distributions.png

```
## 🚀 Usage & Reproducibility
To reproduce the analysis, please note that the ipynb (03) requires a specific environment configuration to support triceratops.

1. Clone the repository
Bash
```text
git clone https://github.com/biesro/TESS-TOI-7475.01-Validation.git
cd TESS-TOI-7475.01-Validation
```
2. Running the Python Scripts (Modeling + Characterization)
For general scripts:
Bash
```text
pip install -r requirements.txt
python code/01_detection_BLS.py
```
### 3. Running the Validation Notebook (Triceratops)
**⚠️ IMPORTANT:** To use Triceratops I recommend following steps in the readme.md of https://github.com/JGB276/TRICERATOPS-plus/tree/main and using jupyter lab (in an isolated python 3.10 environment).

## 📄 Citation
If you use this data or methodology, please cite the arXiv paper. If you use the specific code pipeline, you may also cite the software record.

**Paper (BibTeX):**
```bibtex
@article{escola2026toi7475,
  title={Statistical Validation and Photometric Characterization of the Hot Jupiter Candidate TOI 7475.01},
  author={Escolà Rodrigo, Biel},
  journal={arXiv preprint arXiv:2602.14840 },
  year={2026},
  url={https://arxiv.org/abs/2602.14840 }
}
```
**Software (BibTeX):**
```bibtex

@software{escola2026code,
  author       = {Biel Escolà-Rodrigo},
  title        = {TESS-TOI-7475.01-Validation: Analysis Pipeline},
  year         = 2026,
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.18641713},
  url          = {https://doi.org/10.5281/zenodo.18641713}
}
```

This research made use of the NASA Exoplanet Archive and TESS mission data.

**⚠️ IMPORTANT:** There was some corrections on the results (Temperature and derivations) that must be seen at corrected_characterization_summary.txt for a full comprehension of the case and planet due to a metodological error.
