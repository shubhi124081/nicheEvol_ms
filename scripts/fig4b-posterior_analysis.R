library(ggplot2)

root <- "~/nicheEvol_ms"
exp <- "nicheEvol_Coeligena_newprior5"
# Two empirical experiments
# 4 (C) "nicheEvol_Coeligena_newprior5"
# 4 (D) "nicheEvol_Opisthoprora_newprior5"

res_files <- dir(file.path(root, "res", exp))
contents <- load(file.path(root, "res", exp, res_files[1]))
data_files <- dir(file.path(root, "data", exp))
contents <- load(file.path(root, "data", exp, data_files[1]))
fit <- as.matrix(all$fit)
tr <- everything$data$tr
dist_phylo <- ape::cophenetic.phylo(tr)

dist <- seq(0, max(tr$edge.length), 0.1)
cov_df <- calculate_covariance_realDATA(fit, everything, dist)

# Plot
ggplot(
    cov_df,
    aes(x = dist, y = cov_50, ymin = cov_05, ymax = cov_95)
) +
    geom_line(col = "dodgerblue4", lwd = 1) +
    geom_ribbon(alpha = 0.2, fill = "dodgerblue") +
    geom_line(aes(x = dist, y = cov_fixed),
        col = "red", linetype = "dashed"
    ) +
    ylim(c(0, 0.6)) +
    theme_bw() +
    ylab("Covariance") +
    xlab("Phylogenetic distance") +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    )
