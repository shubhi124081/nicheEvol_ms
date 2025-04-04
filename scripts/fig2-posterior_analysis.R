# Set-up ---
library(ggplot2)

root <- "your_file_path_here"
scripts_directory <- "your_file_path_here"
results_directory <- "your_file_path_here"
data_directory <- "your_file_path_here"
source(file.path(scripts_directory, "00-functions.R"))

exp <- "rafiki2_fig2_OU_as_0e21_1"
BM <- FALSE # The model is Brownian for sub figure 1 ONLY
# NOTE: Commented below are the experiment names and the y axes limits used in
# the manuscript figure:
# Make sure to swap out the experiment name and y axes limits as needed
# Sub figure 1: "rafiki2_fig2_BM_as_0_1", ylim= c(0, 3)
# Sub figure 2: "rafiki2_fig2_OU_siga_25_1_r2", ylim = c(0, 10)
# Sub figure 3: "rafiki2_fig2_OU_siga_1_2", ylim = c(0, 2)
# Sub figure 4: "rafiki2_fig2_OU_as_0e21_1", ylim = c(0, 0.1)

res_files <- dir(file.path(results_directory, "fig2", exp))
contents <- load(file.path(results_directory, "res/fig2", exp, res_files[1]))
data_files <- dir(file.path(data_directory, "data", exp))
contents <- load(file.path(data_directory, "data", exp, data_files[1]))
fit <- as.matrix(all$fit)
tr <- everything$data$tr
dist_phylo <- ape::cophenetic.phylo(tr)

dist <- seq(0, 1, 0.01)
if (BM) {
  cov_df <- calculate_covariance_BM(fit, everything, dist)
} else {
  cov_df <- calculate_covariance(fit, everything, dist)
}

# Plot
ggplot(
  cov_df,
  aes(x = dist, y = cov_50, ymin = cov_05, ymax = cov_95)
) +
  geom_line(col = "dodgerblue4", lwd = 1) +
  geom_ribbon(alpha = 0.2, fill = "dodgerblue") +
  geom_line(aes(x = dist, y = cov_fixed),
    col = "red" # , linetype = "dashed"
  ) +
  ylim(c(0, 0.1)) +
  theme_bw() +
  ylab("Covariance") +
  xlab("Phylogenetic distance") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
