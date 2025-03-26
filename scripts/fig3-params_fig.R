# Set up ---
# Load required libraries
library(ggplot2)

# Required functions
stretchPal <- function(user_cols, num_cols, interpolation = "linear") {
  colfunc <- grDevices::colorRampPalette(user_cols, interpolate = interpolation)
  hex <- colfunc(num_cols)
  return(hex)
}

# NOTE: Code for generating the summary file for this experiment is at the end
# of this script

# FIGURE 3A --------
# File directory
root <- "~/Downloads/new_nicheEvolMs" # Change to your directory
filename <- "fig3a_betaEst_summary.csv"
# Read the summary file for Figure 3A
df <- read.csv(file.path(root, filename))

# Filter and shuffle data for plotting
df3 <- df[-c(which(grepl("BM_1.5_0", df$file2))), ]

# Plot beta estimates vs true values
ggplot(df3) +
  geom_point(aes(x = est, y = beta_true, color = str), size = 2) +
  geom_abline(slope = 1, linetype = "dashed") +
  theme_bw() +
  scale_color_manual(values = c("#b9d6bf", "#577994", "#b5c0da", "#5b9178")) +
  theme(panel.grid = element_blank()) +
  theme(legend.position = "none")

# FIGURE 3B --------
# Set root directory
filename <- "fig3b_rhoEst_summary.csv"
df <- read.csv(file.path(root, filename))

# Create line data for plot
line <- data.frame(y = seq(0, 3, 0.0001))
line$x <- 1 / line$y
line$grp <- rep("1", nrow(line))

# Generate colors for rho plot
rho_colors <- c("#b9d6bf", "#91b3b4", "#577994")
rho_colors <- stretchPal(rho_colors, length(unique(as.factor(df$true))))

# Plot rho values
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
  ylim(c(0, 70)) +
  scale_fill_manual(values = c(rho_colors)) +
  scale_color_manual(values = c(rho_colors)) +
  theme(legend.position = "none")

# FIGURE 3C --------
filename <- "fig3c_omegaEst_summary.csv"
df <- read.csv(file.path(root, filename))

# Generate colors for omega plot
omega_colors <- c("#b5c0da", "#91b3b4", "#5b9178")
omega_colors <- stretchPal(omega_colors, length(unique(as.factor(df$true))))

# Plot omega values
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
  scale_fill_manual(values = c(omega_colors)) +
  scale_color_manual(values = c(omega_colors)) +
  theme(legend.position = "none")

# BELOW ARE THE CODES FOR GENERATING THE SUMMARY FILES FOR FIGURES 3A, 3B, AND 3C
# THESE CODES ARE COMMENTED OUT BECAUSE THEY ARE NOT NEEDED FOR THE PLOTS

# For Figure 3a - IF USING RAW DATA

# Beta estimation figure
# Set root directory and file path
# Follow the filestructure listed in the README
# root <- "~/nicheEvol/"
# fp <- file.path(root, "res")
# exp <- "rafiki"

# # Get list of files matching the experiment
# files <- dir(fp)
# files <- files[grep(exp, files)]
# rmv <- grep("\\.tar.gz", files)
# files <- files[-rmv]

# # Initialize lists for beta estimates and true values
# beta_est <- list()
# beta_true <- list()
# for (i in seq_len(length(files))) {
#   # Load estimated beta values
#   files_in <- dir(file.path(root, "res", files[i]))
#   obj <- load(file.path(root, "res", files[i], files_in[1]))
#   fit <- as.matrix(all$fit)
#   x <- colMeans(fit[, grep("beta", colnames(fit))])
#   estB <- matrix(c(
#     x[seq(1, length(x) - 1, 2)],
#     x[seq(2, length(x), 2)]
#   ), ncol = 2)
#   estB <- matrix(estB, ncol = 1)
#   beta_est[[i]] <- data.frame(
#     "est" = estB,
#     "file" = rep(files_in[1], nrow(estB))
#   )

#   # Load true beta values
#   files_in_data <- dir(file.path(root, "data", files[i]))
#   contents <- load(file.path(root, "data", files[i], files_in_data[1]))
#   beta_true[[i]] <- matrix(everything$data$true_betas, ncol = 1)
# }

# # Combine beta estimates and true values into a data frame
# beta_est2 <- do.call(rbind, beta_est)
# beta_est2$file2 <- lapply(strsplit(beta_est2$file, "_"), function(x) {
#   if (any(grepl("real", x))) {
#     return(paste0(x[2], "_", x[10], "_", x[11]))
#   } else {
#     if (any(grepl("r2", x))) {
#       return(paste0(x[2], "_", x[10], "_", x[11]))
#     }
#     return(paste0(x[2], "_", x[9], "_", x[10]))
#   }
# })
# beta_true <- do.call(rbind, beta_true)
# df <- cbind(beta_true, beta_est2)

# # Map file names to structure types
# dict <- data.frame("id" = do.call(c, unique(df$file2)))
# dict$str <- c(
#   rep("Structured", 12),
#   rep("Not structured", 1), rep("Clustered", 4),
#   rep("Structured", 1), rep("Less structured", 10),
#   rep("Structured", 10), rep("Clustered", 11)
# )

# df$str <- dict[match(df$file2, dict$id), "str"]
# df$file <- as.character(df$file)
# df$file2 <- as.character(df$file2)
# write.csv(df, "~/Downloads/fig3a_summary.csv")


# For Figure 3b - IF USING RAW DATA
# Set file path
# fp <- file.path(root, "zazu", "summaries")

# # Set parameter to search for in filenames
# param <- "rho_rho_zazu"

# # Get list of files matching the parameter
# files <- dir(fp)
# files <- files[grep(param, files)]

# # Extract model, alpha, and sigma2 from filenames
# x <- lapply(strsplit(files, "_"), function(x) {
#   return(c(x[3], x[7], x[8]))
# })
# x <- do.call(rbind, x)
# colnames(x) <- c("model", "alpha", "sigma2")
# x[, 3] <- gsub("\\.csv", "", x[, 3])

# # Note, rho is summary files refers to alpha in manuscript
# # Read alpha values from files
# alpha_vals <- list()
# for (i in seq_len(length(files))) {
#   w <- read.csv(paste0(fp, "/", files[i]))
#   alpha_vals[[i]] <- w$V1
# }
# alpha_val2 <- do.call(c, alpha_vals)

# # Create data frame for rho plot
# df <- data.frame(
#   "est_raw" = alpha_val2,
#   "true" = rep(x[, 2], each = 5),
#   "model" = rep(x[, 1], each = 5)
# )
# df$est <- 1 / df$est_raw
# df <- df[-c(which(df$true == 0)), ]
# df[, "est"] <- df[, "est"] + abs(rnorm(nrow(df), 0, 0.001))
# df[, "est_raw"] <- df[, "est_raw"] + abs(rnorm(nrow(df), 0, 1))

# # Map true values to numeric values
# dict2 <- data.frame("id" = unique(df$true))
# dict2$num <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.5, 2, 2.5, 3)
# df$true <- dict2[match(df$true, dict2$id), "num"]


# Figure 3c - IF USING RAW DATA
# Set parameter for alpha plot
# param <- "alpha_rho_zazu"

# # Get list of files matching the parameter
# files <- dir(fp)
# files <- files[grep(param, files)]

# # Extract model, alpha, and sigma2 from filenames
# x <- lapply(strsplit(files, "_"), function(x) {
#   return(c(x[3], x[7], x[8]))
# })
# x <- do.call(rbind, x)
# colnames(x) <- c("model", "alpha", "sigma2")
# x[, 3] <- gsub("\\.csv", "", x[, 3])

# # Note, rho is summary files refers to alpha in manuscript
# # Read alpha values from files
# alpha_vals <- list()
# for (i in seq_len(length(files))) {
#   w <- read.csv(paste0(fp, "/", files[i]))
#   alpha_vals[[i]] <- w$V1
# }
# alpha_val2 <- do.call(c, alpha_vals)

# # Create data frame for alpha plot
# df <- data.frame(
#   "est" = alpha_val2,
#   "true" = rep(x[, 3], each = 5),
#   "model" = rep(x[, 1], each = 5)
# )

# # Filter and map true values to numeric values
# df <- df[-c(which(df$true %in% c("0e21", "0e31", "0e41"))), ]
# df <- df[1:30, ]
# dict2 <- data.frame("id" = unique(df$true))
# dict2$num <- c(0.1, 0.5, 1, 1.5, 2)
# df$true <- dict2[match(df$true, dict2$id), "num"]
