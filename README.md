# Module Gene Overlap Analysis

A web application for testing gene overlap between user-defined modules and reference modules across datasets. It helps researchers identify statistically significant overlaps and explore shared genes.

**Live app:** https://overlapmodule.onrender.com

![App screenshot](https://github.com/user-attachments/assets/379064a4-0676-46cf-a9cf-90c90fd729c5)

## About

This tool was developed as part of:

**Li H, Xie Z, Tian Y, et al.**  
*Genome-wide consensus transcriptional signatures identify synaptic pruning linking Alzheimer’s disease and epilepsy.*  
**Molecular Psychiatry** (2026) 31:1774–1784  
DOI: 10.1038/s41380-025-03318-0

It enables researchers to validate and compare gene modules against the consensus modules described in the study.

## Features

- Upload gene-module data (**CSV** or **Excel**)
- Compare with built-in or custom reference datasets
- Statistical testing:
  - Fisher’s Exact Test
  - Hypergeometric Test
- Multiple testing correction (Benjamini–Hochberg)
- Interactive heatmap visualization
- Adjustable thresholds and filters
- Download overlapping gene lists

## Installation

```bash
git clone https://github.com/HuihongLi/OverlapModule.git
cd OverlapModule
python -m venv venv
source venv/bin/activate   # Windows: venv\Scripts\activate
pip install -r requirements.txt
mkdir -p data
