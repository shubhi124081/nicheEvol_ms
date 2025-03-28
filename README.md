<<<<<<< HEAD
# Measuring the Evolution of n-Dimensional Environmental Niches

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15090785.svg)](https://doi.org/10.5281/zenodo.15090785)
[![Published in Ecography](https://img.shields.io/badge/Published%20In-Ecography-4b8bbe)](https://doi.org/10.1111/ecog.07285)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![GitHub Repo stars](https://img.shields.io/github/stars/shubhi124081/nicheEvol_ms?style=social)](https://github.com/shubhi124081/nicheEvol_ms/stargazers)
[![GitHub forks](https://img.shields.io/github/forks/shubhi124081/nicheEvol_ms?style=social)](https://github.com/shubhi124081/nicheEvol_ms/network/members)
[![GitHub issues](https://img.shields.io/github/issues/shubhi124081/nicheEvol_ms)](https://github.com/shubhi124081/nicheEvol_ms/issues)
[![Last Commit](https://img.shields.io/github/last-commit/shubhi124081/nicheEvol_ms)](https://github.com/shubhi124081/nicheEvol_ms/commits/main)

_This repository accompanies the manuscript:_  
**"Measuring the evolution of n-dimensional environmental niches"**  
Published in **Ecography (2024)**  
[DOI: 10.1111/ecog.07285](https://doi.org/10.1111/ecog.07285)

---

## 📄 Citation

If you use this code or build upon this work, please cite:

```bibtex
@article{sharma2024nicheevol,
  title = {Measuring the evolution of n-dimensional environmental niches},
  author = {Sharma, Shubhi, Winner, Kevin, Mäkinen, Jussi and Jetz, Walter},
  journal = {Ecography},
  year = {2024},
  doi = {10.1111/ecog.07285}
}
```

---

## Data Availability

All data and results have been archived and are available via Zenodo:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15090785.svg)](https://doi.org/10.5281/zenodo.15090785)

---

## Summary

This manuscript presents a novel species distribution modeling framework that connects spatial and phylogenetic data to jointly estimate both environmental niches and their underlying evolutionary processes.

- The model is a **Bayesian hierarchical framework** that integrates species distribution modeling with phylogenetics.
- Environmental niches are estimated across species using a **latent Gaussian process**, linking **niche similarity to phylogenetic distance**.
- We connect the estimated model parameters to classical **evolutionary models** such as **Ornstein–Uhlenbeck** and **Brownian Motion**, enabling specific hypothesis testing about niche evolution.
- Simulations demonstrate that the framework improves niche estimation across a clade, even when data is sparse or noisy.

---

## Installation

This analysis is lightweight in terms of package dependencies, as most functions required are written in base R (see `scripts/00-functions.R`).

You will need the following R packages:
```r
install.packages(c("rstan", "ape", "phytools", "terra", "MASS", "grDevices", "ggplot2"))
```

We also recommend using [`renv`](https://rstudio.github.io/renv/) for reproducible environments. If available:
```r
renv::restore()
```

---

## 📁 Repository Structure

```text
nicheEvol_ms/
├── data/                 # Simulated data for experiments
│   └── trees/            # Simulated phylogenetic trees
├── raw_data/             # (To be populated) Real-world occurrence data
│   └── trees/            # (To be populated) Real phylogenetic trees
├── res/                  # Output directory for results
├── scripts/              # Core scripts for model and analysis
│   ├── 00-functions.R         # Core functions for modeling and plotting
│   ├── 01-gen_tree.R          # Simulated tree generation (skip if replicating paper)
│   ├── 02-gen_data.R          # Simulated data generation
│   ├── 03-stan_model.R        # Stan model specification
│   ├── 04-run_file.R          # Model fitting and result saving
│   ├── niche_space.R          # Script for generating Figure 4 (part 1)
│   ├── params_fig.R           # Script for generating Figure 2
│   ├── posterior_analysis.R   # Script for generating Figure 3 & Figure 4 (part 2)
├── Appendix_2.Rmd        # Simulation framework vignette
├── .lintr                # Linting configuration
├── .gitignore
└── README.md             # You're here!
```

---

## Reproducing Results

### 1. Clone the repository
```bash
git clone https://github.com/shubhi124081/nicheEvol_ms.git
cd nicheEvol_ms
```

### 2. Install R dependencies

As noted above, either manually install required packages or use `renv`.

### 3. Run simulations

You can replicate the simulations and model fitting pipeline using:

```r
source("scripts/00-functions.R")
source("scripts/02-gen_data.R")
source("scripts/03-stan_model.R")
source("scripts/04-run_file.R")
```

> ⚠️ Note: All models were originally run on a High Performance Cluster using 8 CPUs. Running the full pipeline on a personal machine may be computationally intensive. Precomputed results are available on [Zenodo](https://doi.org/10.5281/zenodo.15090785).

---

## Figures & Visualizations

Each figure in the manuscript can be regenerated using the corresponding script:

- **Figure 2:** `scripts/fig2-posterior_analysis.R`
- **Figure 3:** `scripts/fig3-params_fig.R`
- **Figure 4 (part 1):** `scripts/fig4a-niche_space.R`
- **Figure 4 (part 2):** `scripts/fig4b-posterior_analysis.R`
- **Simulation walkthrough:** `Appendix_2.Rmd`

---

## Contributing

This project is not currently accepting outside contributions, but feel free to fork or open an issue if you have questions or suggestions.

---

## Contact

For questions, feedback, or collaboration inquiries:  
📧 [shubhi.sharma@yale.edu]  

---

