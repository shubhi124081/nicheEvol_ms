# nicheEvol_ms
Code from a manuscript titled "Measuring the evolution of n-dimensional environmental niches", currently in review at Ecography

# Introduction 

In our manuscript, we present a novel species’ distribution modeling framework that connects spatial and phylogenetic data to jointly estimate species’ niches and the underlying niche evolution model. 
The environmental niches of species, i.e. the set of abiotic conditions in which they can sustain a population, underpin the distribution of species in space and time. Understanding how such environmental niches evolve is therefore crucial for addressing some of the biggest challenges in ecology and evolutionary biology; which species are most at risk due to changing abiotic conditions? Will species be able to persist in the face of changing conditions? How will species respond to climate change? 

We introduce a joint framework that simultaneously estimates species’ *n*-dimensional environmental niches together with the underlying evolutionary process. The framework follows a Bayesian hierarchical structure where the relationship between species and environment is estimated while accounting for errors and biases in the raw occurrence datasets. Species’ responses to environmental variables are fitted jointly for the whole clade of species as a latent Gaussian process in which the relationship between niche similarity and phylogenetic distance is estimated. The advantage of this joint framework is two-fold: i. the evolutionary process (distance decay of niche similarity in phylogenetic space) is estimated jointly with the actual niches propagating uncertainty from occurrence data to model parameter estimates; ii. the approach integrates species distribution modeling with phylogenetics effectively connecting robust niche estimation from noisy occurrence data to measuring the rate and mode of niche evolution through traditional evolutionary models such as Ornstein-Uhlenbeck and Brownian Motion that. We derive mathematically how species distribution modeling can be connected to OU and BM models. Further, we show how the estimated parameters of the presented joint framework can enable specific hypothesis tests regarding niche evolution. Through a series of simulations, we emphasize how niche estimates can be improved across a species clade with our framework.

# Framework overview 

We provide all the scripts and functions required to replicate our analyses in the manuscript. Below is a description of the directory structure 

- nicheEvol_ms 
    - data
        - trees
    - raw_data
        - trees
    - res 
    - scripts 
        - 00-functions.R
        - 01-gen_tree.R 
        - 02-gen_data.R 
        - 03-stan_model.R 
        - 04-run_file.R
        - 05-analysis.R 
        - Appendix_2.Rmd
        - niche_space.R 
        - params_fig.R 
        - posterior_analysis.R 
- .lintR
- .gitignore 

# Installation 

# Reproduce figures and results  