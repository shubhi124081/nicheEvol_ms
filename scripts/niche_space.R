# devtools::install_github("lvmuyang/MVNH")
DPATH <- "~/phylo-sdms/phyloproj/raw_data/hummingbirdSA"
NAME <- "hummingbird_sa"


X <- read.csv(file = paste0(DPATH, "/x_raw.csv"))
X <- X[, 1:7]
Y <- read.csv(file = paste0(DPATH, "/y.csv"))
COOD <- read.csv(file = paste0(DPATH, "/cood.csv"))
contents <- load(file = paste0(DPATH, "/tree.Rdata")) # loads obj tree
# TREE <- tree
dist <- ape::cophenetic.phylo(tree)

# For right now, just going to use temperature and precipQuart

X2 <- X[, c("meanTemp", "precipQuart")]

# annotates Y i.e. returns covariates for every presence point
annotateY <- function(x = X, y = Y, focal) {
  y_main <- y[, focal]
  y1_main <- which(y_main == 1)
  xcol <- ncol(x)
  xy_main <- cbind(x, y_main)
  xy_main <- xy_main[y1_main, ]
  if (length(y1_main) == 0) xonly <- NULL
  if (length(y1_main) == 1) xonly <- xy_main[c(1:xcol)]
  if (length(y1_main) > 1) xonly <- xy_main[, c(1:xcol)]

  return(list("xonly" = xonly, "pres_id" = y1_main))
}

focalsp <- "Opisthoprora_euryptera"
distf <- dist[, focalsp]
# tree2 <- ape::keep.tip(tree, grep("Coeligena", tree$tip.label))

# We're going to plot genus by genus
focalg <- "Opisthoprora"
# focalg_nm <- tree2$tip.label
focalg_nm <- c(
  "Opisthoprora_euryptera", "Oreotrochilus_adela",
  "Oreotrochilus_chimborazo", "Oreotrochilus_estella",
  "Oreotrochilus_leucopleurus", "Oreotrochilus_melanogaster", 
  "Polyonymus_caroli"
)
tree1 <- ape::keep.tip(tree, which(tree$tip.label %in% focalg_nm))

focal_a <- list()
for (i in seq_len(length(focalg_nm))) {
  w <- annotateY(x = X2, y = Y, focal = focalg_nm[i])
  focal_a[[i]] <- cbind(rep(focalg_nm[i], nrow(w$xonly)), w$xonly, w$pres_id)
}

a <- do.call(rbind, focal_a)
colnames(a) <- c("sp", "meanTemp", "precipQuart", "pres_id")

X3 <- X2
X3$high <- rep("none", nrow(X3))
X3$high[a$pres_id] <- a$sp

# ggplot() +
#   stat_density_2d(data = a[a$sp == "Coeligena_torquata", ], aes(x = meanTemp, y = precipQuart,
#   fill = after_stat(level)), geom = "polygon") +
#   theme_bw()

# ggplot() +
# stat_density_2d(data = X2, aes(x = meanTemp, y = precipQuart), bins = 30) +
# stat_density_2d(data = a2, aes(x = meanTemp, y = precipQuart,
#   fill = after_stat(level)), geom = "polygon", bins = 3) +
#   theme_bw()

order_names <- data.frame("sp" = tree1$tip.label, 
"order" = LETTERS[seq_len(length(tree1$tip.label))])
order_names$sync <- paste0(order_names$order, "-", order_names$sp)
order_names <- rbind(order_names, c(rep("none", ncol(order_names))))
X3$high <- order_names[match(X3$high, order_names$sp), "sync"]


ggplot(data = X3) +
  geom_point(aes(x = meanTemp, y = precipQuart), col = "grey90", size = 2) +
  geom_point(
    data = X3[X3$high != "none", ],
    aes(x = meanTemp, y = precipQuart, col = high), alpha = 0.5
  ) +
  stat_ellipse(
    data = X3[X3$high != "none", ],
    aes(x = meanTemp, y = precipQuart, col = high), lwd = 1
  ) +
  theme_bw() + 
  theme(panel.grid = element_blank())
