# astrometrica-log-analyzer

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21342573.svg)](https://doi.org/10.5281/zenodo.21342573)
![Python](https://img.shields.io/badge/Python-3.10+-blue)
![License](https://img.shields.io/badge/License-MIT-green)

A quantitative framework for evaluating Astrometrica asteroid observations through automated extraction and visualization of astrometric, photometric, calibration, and kinematic quality metrics.

---

## Overview

> **A Quantitative Framework for Evaluating Asteroid Observations from Astrometrica Log Files**

Astrometrica produces extensive reduction logs containing valuable information about the quality of asteroid observations. However, these outputs are primarily intended for manual inspection and provide no integrated mechanism for quantitatively assessing observation quality.

This project implements a reproducible analysis pipeline that automatically parses Astrometrica's standard output files and extracts observational quality metrics for every detected object.

The framework operates directly on

- `MPCReport.txt`
- `PhotReport.txt`
- `Astrometrica.log`

and generates publication-quality figures, summary statistics, and quality diagnostics.

The framework is independent of any particular observing session and can be applied directly to any IASC Astrometrica reduction.

---

## Features

The single Python source file performs the complete analysis presented in the accompanying paper.

It automatically

- parses all three Astrometrica output files
- synchronizes observations using timestamps
- extracts ten astrometric, photometric, calibration, and kinematic metrics
- computes derived angular motion
- generates publication-quality plots
- produces summary statistics
- evaluates plate calibration stability
- evaluates astrometric residuals
- quantifies signal-to-noise evolution
- analyzes observing conditions throughout an observing session

Every figure appearing in the paper is generated automatically.

---

## Repository Structure

```
astrometrica-log-analyzer/
│
├── astrometrica_report_parser.py
├── README.md
└── LICENSE
```

---

## Requirements

Python 3.10+

Required packages

```
numpy
pandas
matplotlib
scipy
```

Install using

```bash
pip install numpy pandas matplotlib scipy
```

---

## Running the Code

Simply execute

```bash
python astrometrica_report_parser.py
```

The script parses Astrometrica's output files, computes all derived quality metrics, and reproduces every figure and table presented in the accompanying publication.

---

## Generated Figures

Running the script generates the figures presented in the paper, including

- Figure 1 — Apparent Sky Motion and Linear Regression Model
- Figure 2 — Photometric, Image Quality, Astrometric, and Kinematic Metrics
- Figure 3 — Astrometric Plate Solution Residuals
- Figure 4 — Plate Calibration Stability

---

## Extracted Metrics

The parser automatically extracts and analyzes ten quantitative quality metrics from Astrometrica's standard output files.

| Metric | Category | Purpose |
|---------|----------|---------|
| Magnitude | Photometric | Apparent brightness of the detected object. |
| Flux | Photometric | Measured signal from the object. |
| SNR | Inmage Quality | Strength of the detected signal relative to background noise. |
| FWHM | Image Quality | Apparent width of the object's image, indicating focus and seeing. |
| Fit RMS | Astrometric Accuracy | Accuracy of the plate solution. |
| dRA / dDec Residuals | Astrometric Accuracy | Residual errors of the plate solution in right ascension and declination. |
| Angular Velocity | Kinematic | Apparent rate of motion across the sky. |
| Plate Zero Point | Plate Calibration Stability | Calibration constant relating instrumental and standard magnitudes. |
| Calibration Scatter (dmag) | Plate Calibration Stability | Consistency of the reference-star calibration. |
| Reference Stars Used | Plate Calibration Stability | Number of stars contributing to the plate solution. |

## Demonstrated Results

Using an object from an IASC practice image set, the framework demonstrates

- successful extraction of ten observational quality metrics
- automatic generation of publication-ready figures
- trajectory fitting with **R² = 0.9998**
- derived angular velocity of approximately **20.5 arcsec/hour**
- continuous monitoring of photometric, astrometric, and calibration quality throughout the observing session

The framework itself, rather than the demonstration dataset, is the primary contribution.

---

## Supported Input Files

The parser operates directly on Astrometrica's standard output files.

```
MPCReport.txt
PhotReport.txt
Astrometrica.log
```

No modification of these files is required.

---

## License

This project is released under the MIT License.

---

## Citation

If you use this software in research, please cite both the software and the accompanying technical report.

**Software (GitHub):**

```text
Ratnam, A. R. (2026). *astrometrica-log-analyzer* (Version 1.0.0) [Computer software]. GitHub. https://github.com/adityaratnam09/astrometrica-log-analyzer
```

**Technical report (Zenodo):**

```text
Ratnam, A. R. (2026). *A Quantitative Framework for Evaluating Asteroid Observations from Astrometrica Log Files*. Zenodo. https://doi.org/10.5281/zenodo.21342573
```
