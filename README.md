# Gut Microbiome and Parkinson's Disease

## Integrating Machine Learning and Genome-Scale Metabolic Modelling for Microbial and Metabolic Biomarker Discovery

### Project Overview

This project investigates gut microbiome alterations associated with Parkinson's disease (PD) using publicly available microbiome datasets.

The study integrates:

- microbial association analysis
- machine learning
- independent cohort validation
- genome-scale metabolic modelling
- biological interpretation

The primary focus is not simply PD classification, but identification of reproducible microbial signatures and investigation of their potential metabolic consequences.

### Research Question

Which gut microbial features are reproducibly associated with Parkinson's disease across independent cohorts, and what metabolic pathways or metabolites may be associated with these microbial changes?

### Objectives

1. Identify gut microbial features associated with Parkinson's disease using public microbiome datasets.
2. Identify important microbial features using machine learning and evaluate their reproducibility across independent cohorts.
3. Investigate the metabolic potential of key microbial species using available genome-scale metabolic models.
4. Biologically interpret and validate the microbial and metabolic signatures using independent datasets and published literature.

### Core Workflow

Public PD microbiome datasets
→ dataset curation and metadata harmonization
→ quality control / preprocessing
→ taxonomic and functional profiling
→ statistical analysis
→ machine learning
→ independent cohort validation
→ key microbial features
→ genome-scale metabolic modelling
→ predicted pathways/metabolites
→ biological interpretation and validation

### Stretch Goal

A graph neural network (GNN) may be explored to prioritize microbe–metabolite relationships predicted from metabolic modelling.

The GNN is optional and is not required for completion of the core study.

### Scientific Principles

This project will:

- distinguish association from causation
- distinguish predicted metabolic effects from experimentally validated metabolites
- avoid data leakage during machine learning
- keep discovery and validation cohorts separate
- account for study and technical differences
- avoid blindly merging incompatible 16S datasets
- document dataset provenance, software versions and analysis parameters

### Reproducibility

The project will maintain:

- Conda environment specification
- dataset accession and provenance information
- metadata
- analysis scripts
- configuration files
- documented software versions
- fixed random seeds where appropriate
- reproducible analysis steps

Raw sequencing data and other large files will remain outside GitHub.

### Repository Structure

```text
PD_Gut_Microbiome/
├── data/
│   ├── raw/
│   ├── processed/
│   └── metadata/
├── results/
│   ├── qc/
│   ├── taxonomy/
│   ├── statistics/
│   ├── machine_learning/
│   ├── metabolic_modeling/
│   └── gnn/
├── scripts/
├── notebooks/
├── references/
├── docs/
├── environment.yml
├── README.md
├── .gitignore
└── LICENSE

