# Genomic selection and demography in dairy dromedaries

**R code for modeling demography, genomic selection, and inbreeding control in *Camelus dromedarius***

Livestock and Wiledlife laboratory, Arid Regions Institute, Medenine, Tunisia  
Contact: Samia Kdidi, kdidi_samia@yahoo.fr

---

## About this project

Dairy dromedary camels reproduce slowly, which creates a major problem for traditional genetic progress. Welcome to the code repository for our modeling study, where we explore how genomic selection might help overcome these long generation intervals.

This repository contains the R scripts needed to reproduce our entire analytical workflow. The code is structured to perform four main tasks:
1. **Demography:** It builds an age-structured female population model, solving the Euler–Lotka equation for the finite rate of increase (λ) and generation time (T), and verifies the results against a Leslie matrix.
2. **Genomic vs. Progeny Testing:** It uses deterministic genomic prediction and four-path selection theory to calculate the exact reference population size (N*) required for genomic selection to match the annual genetic gain of traditional progeny testing.
3. **Inbreeding Control:** It runs Optimal Contribution Selection (OCS) via quadratic programming on a simulated pedigree to map the trade-off between genetic gain and inbreeding.
4. **Sensitivity Analysis:** It merges the demographic and genetic components into a joint model and performs a variance-based Sobol sensitivity analysis across 15 parameters.

**A quick note on the data:** All parameter values used in this code are illustrative. While they were carefully chosen within ranges reported in the literature for dromedaries, they do not describe one specific real-world population. The pedigree is entirely simulated and no *in vivo* animal data were used. This code is meant to demonstrate how the mathematical model behaves, not to provide immediate breeding recommendations.

## Repository structure

```text
.
├── README.md
├── LICENSE
├── .gitignore
├── data/
│   └── pedigree_demo.csv      # Simulated pedigree used in Module 3
└── src/
    ├── module 12.R           # Modules 1 (Demography) and 2 (Genomics)
    ├── module 3.R             # Module 3 (Optimal Contribution Selection)
    └── module 4.R             # Module 4 (Sobol sensitivity analysis)