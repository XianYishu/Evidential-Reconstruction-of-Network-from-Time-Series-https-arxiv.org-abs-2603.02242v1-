# Evidential Network Reconstruction from SIS Epidemic Time Series

[![MATLAB](https://img.shields.io/badge/MATLAB-R2020b%2B-blue.svg)](https://www.mathworks.com/products/matlab.html)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

A robust network reconstruction framework that infers hidden network structures from multiple independent SIS (Susceptible-Infected-Susceptible) epidemic time series using **Dempster‑Shafer evidence theory**.

## 📌 Overview

Given only binary time series of node states (infected/susceptible) from several independent epidemic outbreaks on an unknown network, this code reconstructs the underlying adjacency matrix. No prior knowledge of the network structure or transmission parameters is required.

**Key features:**
- Multi‑source evidence fusion via Dempster’s rule
- Automatic estimation of infection and recovery rates from data
- Two decision rules for optimal threshold selection:
  - **DR‑MR** (Minimum Robustness) – ensures a connected backbone
  - **DR‑MS** (Maximum Similarity) – maximises Jaccard similarity of new infections
- KL‑divergence based validation of reconstruction quality
- Parallel processing (`parfor`) for multi‑source simulations

## 📁 Repository contents

| File | Description |
|------|-------------|
| `main_Code.m` | Main reconstruction script (run this) |
| `sis_link_new.m` | SIS simulation + first‑level BPA fusion |
| `similarity.m` | Compares reconstructed network with ground truth |
| `compute_infection_jaccard.m` | Jaccard similarity for DR‑MS |
| `estimate_gamma_from_time_series.m` | Estimates recovery rate γ from time series |
| `KL.m` | Kullback‑Leibler divergence between two time series |
| `simulate_sis_with_seeds.m` | SIS simulator with fixed initial infected nodes |
| `network-Karate_adj_matrix.mat` | Zachary’s Karate Club (ground truth for testing) |
| `estimated_gamma.m` | Script to display estimated gamma |
| `KL_Validation.m` | Script to compute KL validation |

## 🚀 Getting started
### Citation
``bib
@article{xian2026evidential,
  title={Evidential Reconstruction of Network from Time Series},
  author={Xian, Yishu and Zhang, Zhaobo and Zhang, Cai and Li, Meizhu and Zhang, Qi},
  journal={arXiv preprint arXiv:2603.02242},
  year={2026}
}

### Prerequisites
- **MATLAB R2020b** or later
- **Parallel Computing Toolbox** (recommended for `parfor`)
- **Statistics and Machine Learning Toolbox** (for `randperm` etc.)

### Installation
1. Clone or download this repository.
2. Ensure all `.m` files and the `.mat` file are in the same folder (or in MATLAB’s path).
3. Open MATLAB, navigate to the folder, and run:
   ```matlab
   main_Code
