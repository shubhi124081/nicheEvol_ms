# Required functions

fixTree <- function(tree_time) {
    tr <- ape::rbdtree(0.1, 0, Tmax = tree_time)
    tr$edge.length <- tr$edge.length / max(phytools::nodeHeights(tr)[
        ,
        2
    ]) * 1
    return(tr)
}

checkTree <- function(TREE1 = tree, ATLEASTN, NOT_MORE_THAN, TREE_TIME) {
    if (!is.null(ATLEASTN) && is.null(NOT_MORE_THAN)) {
        startTime <- Sys.time()
        while (length(tree$tip.label) < ATLEASTN) {
            tree <- fixTree(tree_time = TREE_TIME)
            endTime <- Sys.time()
            if (as.numeric(endTime - startTime) > 120) {
                print("Warning: you might want to consider changing tree_time instead") # nolint
                break
            }
        }
    }
    if (is.null(ATLEASTN) && !is.null(NOT_MORE_THAN)) {
        startTime <- Sys.time()
        while (length(TREE_TIME) > NOT_MORE_THAN) {
            tree <- fixTree(tree_time = TREE_TIME)
            endTime <- Sys.time()
            if (as.numeric(endTime - startTime) > 120) {
                print("Warning: you might want to consider changing tree_time instead") # nolint
                break
            }
        }
    }
    if (!is.null(ATLEASTN) && !is.null(NOT_MORE_THAN)) {
    }
    return(tree)
}


# A couple of functions
simOU_new <- function(
    TREE,
    NTRAITS, SIGMA2, ALPHA) {
    if (ALPHA == 0) {
        bb <- phytools::fastBM(TREE, sig2 = SIGMA2, a0 = 1, nsim = NTRAITS)
    } else {
        bb <- phytools::fastBM(TREE,
            alpha = ALPHA, theta = 1,
            sigma2 = SIGMA2, nsim = NTRAITS
        )
    }
    colnames(bb) <- paste0("trait-", seq_len(ncol(bb)))
    return(bb)
}

GPcovar <- function(TREE, PARAMS) {
    # Collect data and set up empty covariance matrix
    sigma2 <- PARAMS$sigma2$value
    K <- PARAMS$alpha$value
    dist <- ape::cophenetic.phylo(TREE)
    ntips <- length(TREE$tip.label)
    covar <- matrix(0, ncol = ncol(dist), nrow = nrow(dist))
    ntraits <- PARAMS$ntraits

    # GP covariance
    for (i in 1:ntips) {
        for (j in 1:ntips) {
            covar[i, j] <- sigma2 * exp(-(dist[i, j]) / K)
        }
    }

    means <- c(sample(c(3, 4), ntraits, replace = TRUE))

    betas <- matrix(0, ncol = length(means), nrow = ntips)
    for (i in seq_len(length(means))) {
        est <- MASS::mvrnorm(n = 1000, mu = c(rep(means[i], ntips)), covar)
        betas[, i] <- colMeans(est)
    }

    return(betas)
}

traitTrend <- function(
    tree, ntraits, mean_root_vector, covar_traits, sigma2,
    alpha) {
    tr_table <- as.data.frame(tr$edge)
    colnames(tr_table) <- c("node_from", "node_to")
    tr_table$length <- tr$edge.length
    trait_matrix <- matrix(0, ntraits, nrow = nrow(tr_table))
    root_num <- tr_table[1, 1]
    root_traits <- MASS::mvrnorm(1, mean_root_vector, covar_traits)
    for (i in seq_len(nrow(tr_table))) {
        if (tr_table[i, 1] == root_num) {
            next_traits <- alpha * root_traits + MASS::mvrnorm(1,
                mu = rep(0, ntraits), Sigma = diag(sigma2 * tr_table[
                    i,
                    "length"
                ], nrow = ntraits)
            )
            trait_matrix[i, ] <- next_traits
        } else {
            old_state_num <- which(tr_table[i, 1] == tr_table[
                ,
                2
            ])
            old_state <- trait_matrix[old_state_num, ]
            next_traits <- alpha * old_state + MASS::mvrnorm(1,
                mu = rep(0, ntraits), Sigma = diag(sigma2 * tr_table[
                    i,
                    "length"
                ], nrow = ntraits)
            )
            trait_matrix[i, ] <- next_traits
        }
    }
    trait_matrix <- rbind(root_traits, trait_matrix)
    colnames(trait_matrix) <- paste0("trait-", seq_len(ncol(trait_matrix)))
    node_values <- c(root_num, tr_table$node_to)
    trait_df <- as.data.frame(trait_matrix)
    trait_df$node <- node_values

    return(list(trait_df = trait_df, tr = tr))
}

splitByNodesTips <- function(table, tree) {
    Ntips <- length(tree$tip.label)
    Nnodes <- tree$Nnode
    tip_index <- which(table$node %in% 1:Ntips)
    node_index <- which(!(seq_along(nrow(table)) %in% tip_index))
    tip_table <- table[tip_index, ]
    node_table <- table[node_index, ]
    return(list(tip_table = tip_table, node_table = node_table))
}

simulateDataAbundance1 <- function(
    tree, tree_time, ntraits,
    mean_root_vector, covar_traits,
    sigma2, alpha) {
    trait_tree <- traitTrend()

    tr <- trait_tree$tr
    tr$tip.label <- LETTERS[seq_len(length(tr$tip.label))]
    trait_tbl <- trait_tree$trait_df
    ns <- length(tr$tip.label)
    split_tbl <- splitByNodesTips(trait_tbl, tr)
    tip_traits <- split_tbl$tip_table
    rownames(tip_traits) <- tr$tip.label
    return(list(
        tip_traits = tip_traits, trait_tbl = trait_tbl,
        tr = tr
    ))
}

simulateAbundanceData_new <- function(grid,
                                      nsamples = 100,
                                      binary = TRUE,
                                      print = TRUE,
                                      tree,
                                      ntraits,
                                      sigma2,
                                      alpha,
                                      epath,
                                      vars = c("CHELSA_bio_1.tif", "CHELSA_bio_13.tif")) { # nolint

    # Get environmental data
    env <- wrangleEnv(vars1 = vars, epath1 = epath, grid1 = grid)

    # Blank names
    sp <- tree$tip.label

    #### Simulate responses ####

    ns <- length(sp)

    # +1 for intercept - intercept removed
    # mu <- matrix(0, nrow = (nc), ncol = ns)

    # Expected value
    # mu[1, ] <- 0.2 # Intercept
    # mu[1, ] <- tip_traits$`trait-1`
    # mu[2, ] <- tip_traits$`trait-2`

    # Note vcv assumes BM covariance structure
    # For OU a correlation structure needs to be specified
    # L <- ape::vcv(tree)
    # beta <- apply(mu, 1, function(x) {
    #   MASS::mvrnorm(n = 1, mu = x, Sigma = diag(1, ns, ns))
    # # })

    beta <- simOU_new(
        tree,
        ntraits, sigma2, alpha
    )
    # beta <- cbind(rep(1, ns), rep(1, ns))

    beta <- beta + 3

    # Systematic sampling grid first
    Xmat <- makeXmat(env)
    # Xmat <- Xmat[, -c(which(colnames(Xmat) == "Effort"))]
    M <- Xmat %*% t(beta)
    pool2 <- nrow(Xmat)

    # Independently simulating species observations
    covar <- diag(ns)
    diag(covar) <- 2
    Y <- M + MASS::mvrnorm(n = nrow(M), mu = rep(0, ns), Sigma = covar)

    # Fineaggling
    Y2 <- matrix(as.numeric(Y > 0), ncol = ncol(Y))
    Y2 <- data.frame(Y2)
    Y <- data.frame(Y)
    colnames(Y2) <- colnames(Y) <- sp
    Xmat <- data.frame(Xmat)
    Y2$rwid <- Xmat$rwid <- Y$rwid <- paste0("rw_", seq_len(pool2))

    # Create temporary frame so Y can remain the "groud truth"
    Y_temp <- Y2
    Y_temp$tot <- rowSums(Y_temp[, sp])


    k <- sample(seq_len(nrow(Y_temp)), size = nsamples)
    Y_fin <- Y_temp[k, ]


    # More checks
    X_fin <- Xmat[match(Y_fin$rwid, Xmat$rwid), ]
    Yc_fin <- Y[match(Y_fin$rwid, Y$rwid), ]

    # Another check
    if (any(X_fin$rwid != Y_fin$rwid)) {
        warning("Misalignment in X and Y dataframes")
    }

    if (any(nrow(X_fin) != nrow(Y_fin))) {
        warning("Number of rows in X and Y differ")
    }

    if (nrow(Y_fin) != nsamples) {
        warning("Number of observations != number of samples")
    }

    if (nrow(Y2) != nrow(Xmat)) {
        warning("The full datasets are mismatched")
    }

    # Create train and test datasets
    ntrain <- round(0.7 * nsamples, 0)
    ntest <- nsamples - ntrain
    ind_train <- sample(seq_len(nsamples), ntrain)
    Y_train <- Y_fin[ind_train, ]
    Yc_train <- Yc_fin[ind_train, ]
    Y_test <- Y_fin[-ind_train, ]
    Yc_test <- Yc_fin[-ind_train, ]
    X_train <- X_fin[ind_train, ]
    X_test <- X_fin[-ind_train, ]

    if (binary) {
        Yr <- Y2
        Yr_train <- Y_train
        Yr_test <- Y_test
    } else {
        Yr <- Y
        Yr_train <- Yc_train
        Yr_test <- Yc_test
    }

    if (print) {
        printDetails(Y_train[, sp], Y_test[, sp], digits = 4L, width = 72L)
    }

    e <- list(
        "Ytrue" = Yr, "X" = Xmat, "Ytrain" = Yr_train, "Xtrain" = X_train,
        "true_betas" = beta, "train_index" = ind_train,
        "Ytest" = Yr_test, "Xtest" = X_test
    )
    return(e)
}

createGrid <- function(
    USER_GRID, nsamples = 100, USER_PATH, xmin = -10, xmax = 10,
    ymin = -10, ymax = 10) {
    pool <- nsamples * 50
    if (!USER_GRID) {
        pool2 <- ceiling(sqrt(pool))
        x_seq <- seq(xmin, xmax, length.out = pool2)
        y_seq <- seq(ymin, ymax, length.out = pool2)
        sys_cood <- expand.grid(x = x_seq, y = y_seq)
        grid <- sys_cood
    } else {
        load(file = USER_PATH)
        user_cood <- cood[, c("x", "y")]
        rm(cood)
        if (pool > nrow(user_cood)) {
            error("Not enough user input coordinates")
        }
        grid <- user_cood[sample(nrow(user_cood), size = pool), ]
    }
    return(grid)
}

cleanEnv <- function(env1, vars2) {
    w <- which(complete.cases(env1) == FALSE)
    env1 <- env1[-w, ]
    if (any(grepl("\\.tif", vars2))) {
        ind <- grep("\\.tif", vars2)
        for (i in ind) vars2[i] <- gsub("\\.tif", "", vars2[i])
    }
    env1 <- as.matrix(env1[, c(vars2)])
    return(env1)
}

wrangleEnv <- function(
    vars1 = c("CHELSA_bio_1.tif", "CHELSA_bio_13.tif"),
    extent = c(-81, -60, -40, 0),
    epath1,
    grid1) {
    env_files <- file.path(epath1, vars1)
    ras <- terra::rast(env_files)
    env <- terra::extract(ras, grid1)
    env <- cleanEnv(env, vars2 = vars1)
    return(env)
}

makeXmat <- function(env1) {
    x <- scale(env1)
    xmat <- matrix(cbind(x), ncol = ncol(x))
    colnames(xmat) <- c(colnames(x))
    return(xmat)
}

printDetails <- function(..., digits = getOption("digits"), width = getOption("width")) {
    op <- options(digits = digits, width = width)
    on.exit(options(op))
    call <- match.call()
    call[c("digits", "width")] <- NULL
    f1 <- function(x, name) {
        paste0(
            "For ", name, " specs are as follows: Dimensions are ",
            paste0(format(dim(x)), collapse = ", "), " and number of presence points per species are ",
            paste0(format(colSums(x)), collapse = ", "), "."
        )
    }
    l <- Map(f1,
        x = list(...), name = lapply(call[-1], deparse),
        USE.NAMES = FALSE
    )
    s <- do.call(paste, c(l, list(sep = "\n\n")))
    writeLines(strwrap(s))
}

getExtentDf <- function(DF, NAMES) {
    DF2 <- DF[, NAMES]
    if (all(c("lat", "lon") %in% NAMES)) {
        xmin <- min(DF2[, "lon"])
        xmax <- max(DF2[, "lon"])
        ymin <- min(DF2[, "lat"])
        ymax <- max(DF2[, "lat"])
    }
    if (all(c("x", "y") %in% NAMES)) {
        xmin <- min(DF2[, "x"])
        xmax <- max(DF2[, "x"])
        ymin <- min(DF2[, "y"])
        ymax <- max(DF2[, "y"])
    }
    ext <- c(xmin, xmax, ymin, ymax)
    names(ext) <- c("xmin", "xmax", "ymin", "ymax")
    return(ext)
}
