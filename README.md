# Evidential-Reconstruction-of-Network-from-Time-Series-https-arxiv.org-abs-2603.02242v1-
# Evidential Network Reconstruction from SIS Time Series

This repository provides a MATLAB implementation of the evidential network reconstruction framework proposed in:

> *Yishu Xian, Zhaobo Zhang, Cai Zhang, Meizhu Li, Qi Zhang*  
> *“Reconstructing complex networks from time series using Dempster‑Shafer evidence theory”*

The method infers the topology of an unknown network solely from multiple SIS (Susceptible‑Infected‑Susceptible) epidemic time series. It uses Dempster‑Shafer evidence theory to fuse multi‑source uncertain information and reconstructs the network structure without any prior knowledge of the connectivity.

## Features

- **SIS dynamics simulation** with configurable infection rate `β` and recovery rate `γ`
- **Basic Probability Assignment (BPA)** generation from association (`M_A`) and non‑association (`M_N`) matrices
- **Two‑level Dempster fusion**:
  - Level 1: combine BPA from `M_A` and `M_N` for each single time series
  - Level 2: fuse BPAs from multiple independent infection sources
- **Two decision rules** to obtain the final adjacency matrix:
  - **DR‑MR** (Decision Rule based on Minimum Robustness): ensures a connected backbone
  - **DR‑MS** (Decision Rule based on Maximum Similarity): optimises infection‑prediction Jaccard similarity
- **Validation using KL divergence** (relative entropy) to assess reconstruction credibility
- **Parallel computing** support (`parfor`) for efficient processing

## Requirements

- MATLAB R2024b or later
- Parallel Computing Toolbox (optional, for `parfor` loops – code falls back to normal loops if not available)

No additional toolboxes are required.

## File Structure

**Note:** The main script `main_reconstruction.m` includes all helper functions at the bottom. You can use it as a single file, or split it into separate files as described above.

## Usage

### 1. Prepare your own network

Load your ground‑truth adjacency matrix (binary, no self‑loops) as an `n × n` matrix named e.g. `G_graph`.

```matlab
load("your_network.mat")   % contains variable 'A' or rename to 'G_graph'
