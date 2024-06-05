rm(list = ls())

# Some libs
library(ape)
library(phytools)

# Root directory
root <- "~/nicheEvol_ms"

# Source required functions
source(file.path(root, "scripts/00-functions.R"))

# Filepaths
data_directory <- file.path(root, "data")
raw_directory <- file.path(root, "raw_data")

# Config params
SIM <- TRUE
TREE_NAME <- "standard_sim"
TREE_TIME <- 10
MEAN_ROOT_VECTOR <- as.numeric(c(0, 0))
COVAR_TRAITS <- matrix(c(0.5, 1, 1, 2.5), nrow = 2, ncol = 2)
NTRAITS <- 2
SIGMA2 <- 1
ALPHA <- 1
USER_GRID <- FALSE
USER_PATH <- NULL
NSAMPLES <- 100
BINARY <- FALSE
EPATH <- "~/env"

if (SIM) {
    # Load in object tree
    contents <- load(file.path(root, "data/trees", paste0(TREE_NAME, ".Rdata")))


    # Create a spatial grid
    GRID <- createGrid(
        USER_GRID,
        nsamples = 100, USER_PATH,
        xmin = -10, xmax = 10,
        ymin = -10, ymax = 10
    )

    # Simulates an OU/BM process
    all_data <- simulateAbundanceData_new(
        grid = GRID,
        nsamples = NSAMPLES,
        binary = BINARY,
        print = TRUE,
        tree,
        vars = c("CHELSA_bio_1.tif", "CHELSA_bio_13.tif"),
        ntraits = NTRAITS,
        sigma2 = SIGMA2,
        alpha = ALPHA,
        epath = EPATH
    )

    Ytrue <- all_data$Ytrue
    Ytrue <- Ytrue[, -c(which(colnames(Ytrue) == "rwid"))]
    X <- all_data$X
    X <- X[, -c(which(colnames(X) == "rwid"))]

    # Subsetting matrices into test and train
    Ytrain <- all_data$Ytrain
    Ytrain <- Ytrain[, -c(which(colnames(Ytrain) %in% c("rwid", "tot")))]
    Y_occ <- Ytrain
    Y_occ[Y_occ < 0] <- 0
    Xtrain <- all_data$Xtrain
    Xtrain <- Xtrain[, -c(which(colnames(Xtrain) == "rwid"))]
    true_betas <- all_data$true_betas
    train_index <- all_data$train_index
    true_tree <- tree

    # Making matrices for stan code
    ns <- ncol(Ytrain)
    nsamples_train <- nrow(Ytrain)
    L <- matrix(ape::vcv(tree, corr = TRUE), ns, ns)
    distance_matrix <- ape::cophenetic.phylo(tree)
    I <- matrix(0, ns, ns)
    diag(I) <- 1

    I_y <- matrix(0, nsamples_train, nsamples_train)
    diag(I_y) <- 1

    standata <- list(
        N = nsamples_train,
        K = ncol(Xtrain),
        J = ns,
        E = 1,
        N_occ = nsamples_train,
        y = as.matrix(Ytrain),
        y_occ = as.matrix(Y_occ),
        x = Xtrain,
        x_occ = Xtrain,
        x_tilde = Xtrain,
        N_tilde = nrow(Xtrain),
        L_dist = distance_matrix,
        L_corr = L,
        beta_mu = rep(0, ncol(Ytrain)),
        true_betas = true_betas,
        I = I,
        Iy = I_y,
        alpha = 1,
        st_devs = rep(1, ncol(Ytrain)),
        U = 6,
        a = 3,
        b = 0.5,
        tr = tree,
        Y_true = Ytrue,
        X_true = X,
        train_index = train_index,
        rho_mean = 2.5,
        rho_sd = 0.2,
        alpha_mean = 0,
        alpha_sd = 0.2
    )

    # filename to write out
    filename <- paste0("example", ".Rdata")

    # Adding config params to the data file
    config <- list()
    config$sim <- SIM
    config$tree_name <- TREE_NAME
    config$tree_time <- TREE_TIME
    config$mean_root_vector <- MEAN_ROOT_VECTOR
    config$covar_traits <- COVAR_TRAITS
    config$ntraits <- NTRAITS
    config$sigma2 <- SIGMA2
    config$alpha <- ALPHA
    config$user_grid <- USER_GRID
    config$user_path <- USER_PATH
    config$nsamples <- NSAMPLES
    config$binary <- BINARY
    config$epath <- EPATH

    everything <- list(
        data = standata,
        config = config
    )

    # Save
    save(everything, file = paste0(
        data.directory, "/",
        filename
    ))
} else {
    # For non-simulated data

    raw_data_file <- "Opisthoprora_raw"

    # Load in raw annotated data
    contents <- load(file.path(
        raw_directory, paste0(raw_data_file, ".Rdata")
    )) # loads object store

    # Load in tree
    contents <- load(file.path(raw_directory, "trees/Opisthoprora_tree.Rdata"))
    # loads object tree

    if (class(tree) == "phylo") {
        distance_matrix <- ape::cophenetic.phylo(tree)
        # L <- matrix(ape::vcv(tree, corr = T), ncol(Y),ncol(Y))
    } else {
        distance_matrix <- tree
    }

    X <- store$x
    Y <- store$y
    I <- matrix(0, ncol(Y), ncol(Y))
    diag(I) <- 1
    # I_y <- matrix(0, nrow(Y), nrow(Y)); diag(I_y) <- 1

    standata <- list(
        N = nrow(X),
        K = ncol(X),
        J = ncol(Y),
        y = as.matrix(Y),
        x = X,
        # x_tilde = all_data$XData_new,
        # N_tilde = nrow(all_data$XData_new),
        L_dist = distance_matrix,
        # L_corr = L,
        beta_mu = rep(0, ncol(Y)),
        # true_betas = all_data$true_betas,
        I = I,
        # Iy = I_y,
        E = 3,
        alpha = 1,
        st_devs = rep(1, ncol(Y)),
        U = 4,
        a = 1,
        b = 1,
        tr = tree,
        rho_mean = 0,
        rho_sd = 0.4,
        alpha_mean = 0,
        alpha_sd = 0.4
    )
    # tiptraits = tip_traits)

    # Create a data directory and filepath
    filename <- "nicheEvol_Opisthoprora_newprior5"
    everything <- list(
        data = standata,
        config = config
    )

    save(everything, file = paste0(
        data.directory, "/",
        paste0(filename, ".Rdata")
    ))
}
