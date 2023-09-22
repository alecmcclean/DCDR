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
  summarize(avg = mean(MSE)) %>%
  ungroup() %>%
  group_by(n, estimand) %>%
  filter(avg == min(avg)) %>%
  ungroup() %>%
  mutate(estimand = factor(estimand,
                           levels = c("A", "Y", "psi_scdr", "psi_dcdr")))

y_lab <- expression(paste(hat(mu), " estimating ", mu))
a_lab <- expression(paste(hat(pi), " estimating ", pi))
p1 <- plot_data %>%
  ggplot(aes(x = n, y = k, color = estimand)) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  scale_color_colorblind(labels = c("Y" = y_lab, # Nuisance Function Estimation",
                                    "A" = a_lab,
                                    "psi_dcdr" = "DCDR Estimator for ECC",
                                    "psi_scdr" = "SCDR Estimator for ECC")) + 
  scale_x_log10() +
  theme_bw() +
  labs(x = "Fold Size",
       y = "Optimal Number of Neighbors",
       color = "")

ggsave(plot = p1, filename = "../figures/optimal_k_doppler.png",
       height = 4, width = 6)

example <- DGP1(N = 1000, PSI = 0.1)
p2 <- ggplot(data = example, aes(x =X, y = Y)) +
  geom_point(alpha = 0.5) +
  geom_line(aes(x = X, y = mu), color = "red", linewidth = 1) +
  theme_bw()

p2

ggsave(plot = p2, filename = "../figures/doppler_function.png",
       height = 4, width = 5)

# Clean up
rm(list = ls(all = T))
gc()
