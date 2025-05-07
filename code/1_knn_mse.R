###############################################################################
### Author: Alec McClean
### Purpose: Estimate functional varying number of neighbors and dataset size
###############################################################################

set.seed(20230523)
source("0_functions.R")

NUM_ITERS <- 500
output <- expand.grid(n = c(50, 100, 200, 500, 1000, 2000),
                      k = seq(1, 30, 1),
                      ITER = seq(1, NUM_ITERS, 1))


output$A <- output$Y <- output$psi_dcdr <- output$psi_scdr <- NA

for (row in 1:nrow(output)) {
  cat("\nRow : ", row, " out of ", nrow(output))
  MSEs <- GetCCEstimate_kNN(FOLD_SIZE = output$n[row], 
                            K = output$k[row],
                            PSI = 0.1)
  output$A[row] <- MSEs[[1]]
  output$Y[row] <- MSEs[[2]]
  output$psi_dcdr[row] <- MSEs[[3]]
  output$psi_scdr[row] <- MSEs[[4]]
}

plot_data <- output %>%
  gather(estimand, MSE, psi_scdr:A) %>%
  group_by(n, k, estimand) %>%
  summarize(avg = mean(MSE), .groups = "drop") %>%
  group_by(n, estimand) %>%
  filter(avg == min(avg)) %>%
  ungroup() %>%
  mutate(estimand = factor(estimand,
                           levels = c("A", "Y", "psi_scdr", "psi_dcdr")))

y_lab <- expression(hat(mu) ~ " estimating " ~ mu)
a_lab <- expression(hat(pi) ~ " estimating " ~ pi)

p1a <- plot_data %>%
  ggplot(aes(x = n, y = k, shape = estimand, linetype = estimand)) +  # Single mapping for one legend
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  scale_shape_manual(values = c("Y" = 16,  # Filled circle
                                "A" = 17,  # Filled triangle
                                "psi_dcdr" = 15,  # Filled square
                                "psi_scdr" = 18),
                     labels = c("Y" = y_lab,
                                "A" = a_lab,
                                "psi_dcdr" = "DCDR estimator for ECC",
                                "psi_scdr" = "SCDR-MSE estimator for ECC")) +
  scale_linetype_manual(values = c("Y" = "solid",
                                   "A" = "dashed",
                                   "psi_dcdr" = "dotted",
                                   "psi_scdr" = "dotdash"),
                        labels = c("Y" = y_lab,
                                   "A" = a_lab,
                                   "psi_dcdr" = "DCDR estimator for ECC",
                                   "psi_scdr" = "SCDR-MSE estimator for ECC")) +
  scale_x_log10() +
  theme_bw() +
  labs(x = "Fold Size",
       y = "Optimal Number of Neighbors",
       shape = "",  # Single legend title
       linetype = "") +
  guides(shape = guide_legend(
    override.aes = list(linetype = c("solid", "dashed", "dotted", "dotdash")))
  )

ggsave(plot = p1a, filename = "../figures/optimal_k_doppler.png",
       height = 4, width = 6)

plot_data_filtered <- plot_data %>%
  filter(estimand != "A")  

p1b <- plot_data_filtered %>%
  ggplot(aes(x = n, y = k, shape = estimand, linetype = estimand)) +  # Single mapping for one legend
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  scale_shape_manual(values = c("Y" = 16,  # Filled circle
                                "psi_dcdr" = 15,  # Filled square
                                "psi_scdr" = 18),  # Filled diamond
                     labels = c("Y" = "Nuisance estimator",
                                "psi_dcdr" = "DCDR estimator",
                                "psi_scdr" = "SCDR-MSE estimator")) +
  scale_linetype_manual(values = c("Y" = "solid",
                                   "psi_dcdr" = "dotted",
                                   "psi_scdr" = "dotdash"),
                        labels = c("Y" = "Nuisance estimator",
                                   "psi_dcdr" = "DCDR estimator",
                                   "psi_scdr" = "SCDR-MSE estimator")) +
  scale_x_log10() +
  theme_bw() +
  labs(x = "Fold Size",
       y = "Optimal Number of Neighbors to minimize error",
       shape = "",  # Single legend title
       linetype = "") +
  guides(shape = guide_legend(
    override.aes = list(linetype = c("solid", "dotted", "dotdash")))
  )

ggsave(plot = p1b, filename = "../figures/optimal_k_doppler_intro.png",
       height = 4, width = 6)


example <- DGP1(N = 1000, PSI = 0.1)
p2 <- ggplot(data = example, aes(x =X, y = Y)) +
  geom_point(alpha = 0.5) +
  geom_line(aes(x = X, y = mu), color = "black", linewidth = 1) +
  theme_bw()

ggsave(plot = p2, filename = "../figures/doppler_function.png",
       height = 4, width = 5)

# Clean up
rm(list = ls(all = T))
gc()
