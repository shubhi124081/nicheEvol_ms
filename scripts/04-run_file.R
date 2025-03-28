# Set-up -----
root <- "your filepath here"
# You may choose your own file structure, the code below assumes that
# the root directory has the following sub-directories: "data", "scripts", "res"
# However, the zenodo directory file structure is different - it is organized
# by figures and sub-figures. Feel free to re-organize the files as you see fit.
data_directory <- file.path(root, "data")
scripts_directory <- file.path(root, "scripts")
res_directory <- file.path(root, "res")

# Source all the stan models
source(file.path(scripts_directory, "03-stan_model.R"))
dataset <- "your dataset name goes here" # for example, rafiki2_fig2_BM_as_0_1

# Load data
filename <- "your dataset name goes here" # for example, rafiki2_fig2_BM_as_0_1-var_alpha_sigma2_1_0_rep_1.Rdata
contents <- load(file.path(data_directory, dataset, filename))

# Organize file for STAN
standata <- everything$data
saved_args <- everything$config
data_sim <- saved_args$data_wrangling
stan_specs <- saved_args$stan_specs
model <- stan_specs$model
# Checl model
print(model)
# Select model
# If presence-absence data use multiprobit2_eff_off_nicheEvol
# If continuous data use GPtrun2_mod_newprior
code <- as.character(model)

# Run STAN
fit <- rstan::stan(
    model_code = code,
    data = standata,
    iter = stan_specs$iter,
    thin = stan_specs$thin,
    warmup = stan_specs$warmup,
    chains = stan_specs$chains,
    control = list(
        adapt_delta = stan_specs$delta,
        stepsize = stan_specs$stepsize,
        stepsize_jitter = stan_specs$stepsize_jitter
    ),
    cores = 1
)
modelname <- paste0(dataset, "_", stan_specs$model)

# Chuck "useless" columns - file sizes are seriously bloated
fit <- as.matrix(fit)
b_ind <- grep("beta", colnames(fit))
st_ind <- grep("st_devs", colnames(fit))
y_ind <- grep("y_devs", colnames(fit))
r_ind <- grep("rho", colnames(fit))
a_ind <- grep("alpha", colnames(fit))
indices <- c(b_ind, st_ind, y_ind, a_ind, r_ind)
fit <- fit[, indices]

# Create a data directory and filepath
res_dir2 <- file.path(res_directory, dataset)
if (!dir.exists(res_dir2)) {
    dir.create(res_dir2)
}
all <- list(
    fit = fit,
    args = saved_args
    # data = standata
)
save(all, file = file.path(res_dir, paste0(modelname, ".Rdata")))
