library(ggplot2)

stretchPal <- function(user_cols, num_cols, interpolation = "linear") {
  colfunc <- grDevices::colorRampPalette(user_cols, interpolate = interpolation)
  hex <- colfunc(num_cols)
  return(hex)
}

root <- "~/nicheEvol/analysis"

fp <- file.path(root, "zazu", "summaries")

param <- "rho_rho_zazu"

files <- dir(fp)
files <- files[grep(param, files)]

# parameter _ experiment _ model _ parameter1 _ parameter2 _ value1 _ value2

x <- lapply(strsplit(files, "_"), function(x) {
  return(c(x[3], x[7], x[8]))
})
x <- do.call(rbind, x)
colnames(x) <- c("model", "alpha", "sigma2")
x[, 3] <- gsub("\\.csv", "", x[, 3])

alpha_vals <- list()
for (i in seq_len(length(files))) {
  w <- read.csv(paste0(fp, "/", files[i]))
  alpha_vals[[i]] <- w$V1
}
alpha_val2 <- do.call(c, alpha_vals)

df <- data.frame(
  "est" = alpha_val2,
  "true" = rep(x[, 2], each = 5),
  "model" = rep(x[, 1], each = 5)
)

# Rho plot
df <- data.frame(
  "est_raw" = alpha_val2,
  "true" = rep(x[, 2], each = 5),
  "model" = rep(x[, 1], each = 5)
)
df$est <- 1 / df$est_raw
df <- df[-c(which(df$true == 0)), ]
df[, "est"] <- df[, "est"] + abs(rnorm(nrow(df), 0, 0.001))
df[, "est_raw"] <- df[, "est_raw"] + abs(rnorm(nrow(df), 0, 1))

dict2 <- data.frame("id" = unique(df$true))
dict2$num <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.5, 2, 2.5, 3)
df$true <- dict2[match(df$true, dict2$id), "num"]

line <- data.frame(y = seq(0, 3, 0.0001))
line$x <- 1 / line$y
line$grp <- rep("1", nrow(line))

rho_colors <- c("#b9d6bf", "#91b3b4", "#577994")
rho_colors <- stretchPal(rho_colors, length(unique(as.factor(df$true))))

ggplot(df) +
  geom_boxplot(aes(
    x = true, y = est_raw,
    group = true,
    fill = as.factor(true),
    color = as.factor(true)
  )) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_line(
    data = line, aes(x = y, y = x, group = grp),
    linetype = "dashed"
  ) +
  # geom_abline(slope = 1, linetype = "dashed") +
  ylim(c(0, 70)) +
  scale_fill_manual(values = c(rho_colors)) +
  scale_color_manual(values = c(rho_colors)) +
  theme(legend.position = "none")


# Alpha plot
param <- "alpha_rho_zazu"
df <- data.frame(
  "est" = alpha_val2,
  "true" = rep(x[, 3], each = 5),
  "model" = rep(x[, 1], each = 5)
)

df <- df[-c(which(df$true %in% c("0e21", "0e31", "0e41"))), ]
df <- df[1:30, ]
dict2 <- data.frame("id" = unique(df$true))
dict2$num <- c(0.1, 0.5, 1, 1.5, 2)
df$true <- dict2[match(df$true, dict2$id), "num"]



alpha_colors <- c("#b5c0da", "#91b3b4", "#5b9178")
alpha_colors <- stretchPal(alpha_colors, length(unique(as.factor(df$true))))

ggplot(df) +
  geom_boxplot(aes(
    x = true, y = est,
    group = true,
    fill = as.factor(true),
    color = as.factor(true)
  )) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_abline(slope = 1, linetype = "dashed") +
  ylim(c(0, 2)) +
  scale_fill_manual(values = c(alpha_colors)) +
  scale_color_manual(values = c(alpha_colors)) +
  theme(legend.position = "none")


# Beta estimation figure
# I think we should just do all the rafiki experiments
# This might be a huge mess
# Idk let's go

root <- "~/nicheEvol/"
fp <- file.path(root, "res")
exp <- "rafiki"

files <- dir(fp)
files <- files[grep(exp, files)]
rmv <- grep("\\.tar.gz", files)
files <- files[-rmv]

files_data <- dir(file.path(root, "data"))

beta_est <- list()
beta_true <- list()
for (i in seq_len(length(files))) {
  # est
  files_in <- dir(file.path(root, "res", files[i]))
  obj <- load(file.path(root, "res", files[i], files_in[1]))
  fit <- as.matrix(all$fit)
  x <- colMeans(fit[, grep("beta", colnames(fit))])
  estB <- matrix(c(
    x[seq(1, length(x) - 1, 2)],
    x[seq(2, length(x), 2)]
  ), ncol = 2)
  estB <- matrix(estB, ncol = 1)
  beta_est[[i]] <- data.frame(
    "est" = estB,
    "file" = rep(files_in[1], nrow(estB))
  )

  # true
  files_in_data <- dir(file.path(root, "data", files[i]))
  contents <- load(file.path(root, "data", files[i], files_in_data[1]))
  beta_true[[i]] <- matrix(everything$data$true_betas, ncol = 1)
}

beta_est2 <- do.call(rbind, beta_est)
beta_est2$file <- lapply(strsplit(beta_est2$file, "_"), function(x) {
  if (any(grepl("real", x))) {
    return(paste0(x[2], "_", x[10], "_", x[11]))
  } else {
    if (any(grepl("r2", x))) {
      return(paste0(x[2], "_", x[10], "_", x[11]))
    }
    return(paste0(x[2], "_", x[9], "_", x[10]))
  }
})
beta_true <- do.call(rbind, beta_true)
df <- cbind(beta_true, beta_est2)

dict <- data.frame("id" = do.call(c, unique(df$file)))
dict$str <- c(
  rep("Structured", 12),
  rep("Not structured", 1), rep("Clustered", 4),
  rep("Structured", 1), rep("Less structured", 10),
  rep("Structured", 10), rep("Clustered", 11)
)

df$str <- dict[match(df$file, dict$id), "str"]

df3 <- df[-c(which(grepl("BM_1.5_0", df$file))), ]
df3 <- df3[sample(seq_len(nrow((df3)))), ]
ggplot(df3) +
  geom_point(aes(x = est, y = beta_true, color = str), size = 2) +
  geom_abline(slope = 1, linetype = "dashed") +
  theme_bw() +
  scale_color_manual(values = c("#b9d6bf", "#577994", "#b5c0da", "#5b9178")) +
  theme(panel.grid = element_blank()) #+
