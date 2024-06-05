library(bayesplot)
library(ggplot2)

root <- "~/nicheEvol"
exp <- "rafiki2_fig2_OU_siga_1_2"

res_files <- dir(file.path(root, "res", exp))
contents <- load(file.path(root, "res", exp, res_files[1]))
# contents <- load("/Users/ss4224/phylo-sdms/phyloproj/res/grumpy_Opisthoprora_e2021/grumpy_Opisthoprora_e2021_rep1_multiprobit2_effort.Rdata")
# contents <- load("/Users/ss4224/phylo-sdms/phyloproj/res/grumpy_Coeligena_e2021/grumpy_Coeligena_e2021_rep1_multiprobit2_effort.Rdata")
data_files <- dir(file.path(root, "data", exp))
contents <- load(file.path(root, "data", exp, data_files[1]))
fit <- as.matrix(all$fit)
tr <- everything$data$tr
dist_phylo <- ape::cophenetic.phylo(tr)
# Histograms of posterior
bayesplot::color_scheme_set("blue")
mcmc_hist(fit, pars = c("rho", "alpha")) +
  theme_bw()

rho <- mean(fit[, "rho"])
alpha <- mean(fit[, "alpha"])

dist <- seq(0, 1, 0.01)

rho_05 <- quantile(fit[, "rho"], 0.05)
rho_95 <- quantile(fit[, "rho"], 0.95)
rho <- mean(fit[, "rho"])

alpha_05 <- quantile(fit[, "alpha"], 0.05)
alpha_95 <- quantile(fit[, "alpha"], 0.95)
alpha <- mean(fit[, "alpha"])

cov_95 <- matrix(0, nrow = length(dist))
for (i in seq_len(length(dist))) {
  cov_95[i] <- ((alpha_95^2 * exp(-0.5 * dist[i] / rho_95^2)))
}

cov_05 <- matrix(0, nrow = length(dist))
for (i in seq_len(length(dist))) {
  cov_05[i] <- (alpha_05^2 * exp(-0.5 * dist[i] / rho_05^2))
}

cov_50 <- matrix(0, nrow = length(dist))
for (i in seq_len(length(dist))) {
  cov_50[i] <- ((alpha^2 * exp(-0.5 * dist[i] / rho^2)))
}
# From simulation
cov_fixed <- matrix(0, nrow = length(dist))
alpha_fixed <- everything$config$phylo$sigma2$value
# alpha_fixed <- alpha
rho_fixed <- everything$config$phylo$alpha$value
# rho_fixed <- 1
cov_fixed <- matrix(0, nrow = length(dist))
for (i in seq_len(length(dist))) {
  cov_fixed[i] <- alpha_fixed^2 * exp(-rho_fixed * (dist[i]))
  # cov_fixed[i] <- (alpha_fixed^2 * exp(-(dist[i])))
}

cov_df <- data.frame(
  "cov_fixed" = cov_fixed,
  "cov_50" = cov_50,
  "cov_05" = cov_05,
  "cov_95" = cov_95,
  "dist" = dist
)

ggplot(
  cov_df,
  aes(x = dist, y = cov_50, ymin = cov_05, ymax = cov_95)
) +
  geom_line(col = "dodgerblue4", lwd = 1) +
  geom_ribbon(alpha = 0.2, fill = "dodgerblue") +
  geom_line(aes(x = dist, y = cov_fixed),
    col = "red" # , linetype = "dashed"
  ) +
  ylim(c(0, 2)) +
  theme_bw() +
  ylab("Covariance") +
  xlab("Phylogenetic distance") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
