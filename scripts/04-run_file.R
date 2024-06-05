if (Sys.info()["login"] == "root") {
    root <- "~/nicheEvol"
} else {
    root <- "~/project/nicheEvol"
}
data.directory <- file.path(root, "data")
scripts.directory <- file.path(root, "scripts")
res.directory <- file.path(root, "res")
job.directory <- file.path(root, "jobs")
log.directory <- file.path(root, "log")
args <- commandArgs(trailingOnly = TRUE)

print(args)
model <- args[1]
dataset_main <- args[2]
dataset <- args[3]

# Soruce all the stan models
source(file.path(scripts.directory, "023-stan_models.R"))


# Load data
filename <- paste0(dataset, ".Rdata")
load(file.path(data.directory, dataset_main, filename))

standata <- everything$data
saved_args <- everything$config
data_sim <- saved_args$data_wrangling
stan_specs <- saved_args$stan_specs
code <- as.character(mostCommonlyUsed[model])
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
modelname <- paste0(dataset, "_", model)

# Chuck "useless" columns - file sizes are seriously bloated
# fit <- as.matrix(fit)
# b_ind <- grep("beta", colnames(fit))
# st_ind <- grep("st_devs", colnames(fit))
# y_ind <- grep("y_devs", colnames(fit))
# r_ind <- grep("rho", colnames(fit))
# a_ind <- grep("alpha", colnames(fit))
# indices <- c(b_ind, st_ind, y_ind, a_ind, r_ind)
# fit <- fit[, indices]

# Create a data directory and filepath
res.dir <- file.path(res.directory, data_sim$experiment_name)
if (!dir.exists(res.dir)) {
    dir.create(res.dir)
}
all <- list(
    fit = fit,
    args = saved_args
    # data = standata
)
save(all, file = file.path(res.dir, paste0(modelname, ".Rdata")))
