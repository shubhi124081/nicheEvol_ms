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
- .lintr
- .gitignore 
- README.md

**Note**: there is no data stored in this repository. To access the data visit (https://zenodo.org/records/11497929)

The data for simulations is stored in `data` directory with all simulated trees stored in the `data/trees` sub-directory 

The data for real-world analysis will need to be stored in `raw_data` directory with trees stored in the `raw_data/trees` sub-directory

All results will be written to and read from `res` directory

The `scripts` directory stores all functions and scripts to fit the models and do the analysis. Below is a description of all the contents of the `scripts` directory. Overall, the numbered files run the model fitting and analysis. The non-numbered files are used for generating figures and appendicies. 

`00-functions.R` has all of the functions that need to be read in to support model fitting, analysis and plotting 

`01-gen_tree.R` is a script for generating a phylogeny for use in simulations. **Note**: If replicating analysis in the manuscript, skip this script. The simulated phylogeny is uploaded to the Zenodo repository

`02-gen_data.R` script for generating simulated data. For a brief overview of this script take a look at `Appendix_2.Rmd` which is a vigenette going through the simulation framework 

`03-stan_model.R` script with stan model written down 

`04-run_file.R` reads in data and stan models and runs files. Writes out results to `res` directory. **Note**: All models were run on the High Performance Cluster with 8 CPUs. These models may be too heavy to run on a personal machine. Therefore, all the result files have also been made available in the Zenodo repository. 

`05-analysis.R` empty 

`Appendix_2.Rmd` vigenette going through simulation framework 

`niche_space.R` Figure 4 part 1 

`params_fig.R` Figure 2

`posterior_analysis.R` Figure 3, Figure 4 part 2 

# Installation 



# Reproduce figures and results  