# Set-up
rm(list = ls())

# Some libs
library(ape)
library(phytools)
# library(phyloGenie)

# Root directory
root <- "~/nicheEvol_ms"

# Source required functions
source(file.path(root, "scripts/00-functions.R"))

# Simulation parameters
MEAN_ROOT_VECTOR <- as.numeric(c(0, 0))
COVAR_TRAITS <- matrix(c(0.5, 1, 1, 2.5), nrow = 2, ncol = 2)
ATLEASTN <- 8
NOT_MORE_THAN <- NULL
TREE_TIME <- 10
NTIPS <- NULL
TREE_NAME <- "new_tree"
# Code to simulate tree

tree <- fixTree(tree_time = phylo$tree_time)
tree <- checkTree()

print(paste0("# Tips : ", length(tree$tip.label)))
plot(tree)
tree$tip.label <- LETTERS[seq_len(length(tree$tip.label))]

save(tree, file = file.path(
    getwd(), "trees",
    paste0(TREE_NAME, ".Rdata")
))


# Alternative approach - not using birth-death process

if (!is.null(NTIPS)) {
    tree <- ape::rtree(ntips)
}
