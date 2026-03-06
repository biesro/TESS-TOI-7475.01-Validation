# Statistical Validation and Photometric Characterization of TOI 7475.01

**Code repository for the paper:**
> *Statistical Validation and Photometric Characterization of the Hot Jupiter Candidate TOI 7475.01*

![Python Version](https://img.shields.io/badge/python-3.10-blue)
![License](https://img.shields.io/badge/license-MIT-green)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18641713.svg)](https://doi.org/10.5281/zenodo.18641713)
[![Status](https://img.shields.io/badge/Status-Validated%20Planet-success.svg)]()
[![arXiv](https://img.shields.io/badge/arXiv-2602.14840-b31b1b.svg)](https://arxiv.org/abs/2602.14840)

**Principal Investigator:** Biel Escolà-Rodrigo  
**Date:** February 2026  

---

## 🔭 Abstract

This repository contains the complete analysis pipeline used to statistically validate and physically characterize the exoplanet candidate **TOI 7475.01** (TIC 376866659).

Combining a custom photometric pipeline with rigorous vetting steps, we identified a robust transit signal with a period of **3.2538 days** and a depth of **~4600 ppm**. Unlike many candidates that fail due to background contamination, our analysis confirms a clean spatial environment and a False Positive Probability (FPP) of **$\approx 0$**, statistically validating the signal as a bona fide planetary companion.

Key findings supporting this validation include:

* **High Signal-to-Noise:** A detection SNR of **294.13**, indicating an exceptionally clear signal.
* **Clean Environment:** Gaia DR3 analysis shows no contaminating sources within the critical aperture (nearest neighbor at 28.3", $\Delta G \approx 6$).
* **Stable Centroids:** No photocenter shift during transit, confirming the signal originates from the target star.
* **Statistical Validation:** A rigorous `triceratops` analysis (20 independent runs) yields an FPP $< 10^{-6}$, ruling out astrophysical false positives.

This repository provides the code to reproduce the detection, the spatial vetting, and the statistical validation analysis.

## 🌍 Key Findings & Parameters

Based on the analysis of **TESS Sector 91**:

| Parameter | Value | Note |
| :--- | :--- | :--- |
| **Target** | TIC 376866659 | V = 8.56 |
| **Signal Period** | **3.253773 days** | Ultra-stable periodicity |
| **Transit Depth** | ~4601 ppm | Consistent U-shape |
| **Signal-to-Noise** | **294.13** | High-confidence detection |
| **Contamination** | None | Clean aperture (Gaia DR3) |
| **Validation Status** | **FPP $\approx 0$** | Statistically Validated |
| **Classification** | **Confirmed Planet** | Validated via Triceratops |

## ⚙️ Methodology & Software Stack
The analysis follows a rigorous forensic protocol implemented in Python.

### 1. Detection & Signal Recovery
* **Tool:** `Lightkurve` (v2.x) + Custom "Fusion-Elite" Pipeline
* **Finding:** Recovers a periodic signal with a distinct U-shaped transit and no odd/even depth variations (ruling out eclipsing binaries).
* **Script:** `01_detection_BLS.py`

### 2. Spatial Vetting & Contamination
* **Tool:** `astroquery` (Gaia DR3) + `matplotlib`
* **Finding:** The target star is isolated within the TESS pixel scale (~21"). The nearest neighbor is 28.3" away and too faint to mimic the signal.
* **Script:** `01.2_spatial_contamination_check.py`

### 3. Centroid Analysis
* **Tool:** `Lightkurve` (Moments Method)
* **Finding:** Flux-weighted centroids show no significant shift during the transit window, confirming the signal is on-target.
* **Script:** `02_centroid_FINAL.py`

### 4. Statistical Validation (TRICERATOPS)
* **Tool:** `triceratops` (v1.0.18)
* **Finding:** The analysis yields a mean False Positive Probability (FPP) of **0.000000** across 20 independent runs. The probability mass is concentrated in scenarios where the planet orbits the target star (TP) or a bound companion (PTP).
* **Script:** `03_triceratops_vetting.py` (and robustness test)

## 📂 Repository Structure

```text
TOI-7475.01-Validation/
│
├── LICENSE                     # MIT License
├── README.md                   # Project documentation
├── requirements.txt            # Python dependencies
├── TOI7475.01_Folded_Raw_GrayPoints.csv
│
├── code/                    # Analysis scripts
│   ├── 01_detection_BLS.py
│   ├── 01.2_spatial_contamination_check.py
│   ├── 02_centroid_FINAL.py
│   └── 03_triceratops_vetting.py
│
└── figures/                    # Generated plots
    ├── Figure_01_BLS+Transit+OddEven.png
    ├── Figure_02_spatial_contam_check.png
    ├── Figure_03_Centroid.png
    └── TOI7475_triceratops_20runs_histograms.png

```
## 🚀 Usage & Reproducibility
To reproduce the analysis, please note that the ipynb (03) requires a specific environment configuration to support triceratops.

1. Clone the repository
Bash
```text
git clone https://github.com/biesro/TESS-TOI-864.01-Validation.git
cd TESS-TOI-7475.01-Validation
```
2. Running the Python Scripts (Modeling)
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
