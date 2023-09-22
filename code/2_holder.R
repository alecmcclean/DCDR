###############################################################################
### Author: Alec McClean
### Purpose: Analysis of SCDR and DCDR estimators with Holder smooth functions
###############################################################################

set.seed(20230526)
source("0_functions.R")

##################################
### Generate example function

PARAM_DAT <- expand.grid(N = c(100, 1000, 5000),
                         beta = c(0.1, 0.35, 0.6))

data <- data.frame()
for (i in 1:nrow(PARAM_DAT)) {
  data %<>% bind_rows(
    SimulateData(psi = 10, 
                 dimension = 1, 
                 beta = PARAM_DAT$beta[i], 
                 N = PARAM_DAT$N[i], 
                 signal_amp = 50) %>%
      mutate(N = PARAM_DAT$N[i],
             smoothness = PARAM_DAT$beta[i])
  )
}

p0 <- data %>%
  mutate(N_label = paste0("N = ", N),
         smooth_label = paste0("Smoothness = ", smoothness)) %>%
  ggplot(aes(x = X1)) +
  geom_point(aes(y = Y), color = "red", alpha = 0.3) +
  geom_line(aes(y = mu)) +
  theme_bw() + 
  facet_grid(reorder(N_label, N) ~ reorder(smooth_label, smoothness), scales = "free_y")

ggsave(plot = p0, filename = "../figures/example_holder.png",
       width = 8, height = 6)

#################################################
### Analyze several Holder smooth functions
#################################################

### Construct datasets and run estimators
FOLD_SIZES <- c(100 , 200, 500 , 1000, 2000, 5000)
DIMENSIONS <- c(1, 4)
BETAS <- c(0.1, 0.35, 0.6, 1.25, 2.25)
PSIS <- c(10)
NUM_ITERS <- 100
ITERS <- seq(1, NUM_ITERS, 1)

PARAM_DAT <- expand.grid("FOLD_SIZES" = FOLD_SIZES, 
                         "DIMENSIONS" = DIMENSIONS,
                         "BETAS" = BETAS, 
                         "PSIS" = PSIS,
                         "ITERS" = ITERS) %>%
  filter(!(DIMENSIONS == 1 & BETAS > 1),
         !(DIMENSIONS == 4 & BETAS < 0.5)) 

numCores <- detectCores()

cl <- makeCluster(numCores-2)
registerDoSNOW(cl)
pb <- txtProgressBar(max = nrow(PARAM_DAT), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

inf_data <- foreach(
  index = seq(1, nrow(PARAM_DAT), 1),
  .combine=rbind,
  .options.snow = opts) %dopar% {
    GenerateEstimateHolderSmooth(FOLD_SIZE = PARAM_DAT$FOLD_SIZES[index],
                                 DIMENSION = PARAM_DAT$DIMENSIONS[index],
                                 BETA = PARAM_DAT$BETAS[index], 
                                 PSI = PARAM_DAT$PSIS[index], 
                                 ITER_NUMBER = PARAM_DAT$ITERS[index])
  }

parallel::stopCluster(cl)


#################################################
### Create figures to illustrate results 
#################################################

##############################
### Coverage and power

p1 <- inf_data %>%
  filter(!(dimension == 4 & estimator == "unknown")) %>%
  mutate(cilb = psihat - qnorm(0.975) * sqrt(var / n),
         ciub = psihat + qnorm(0.975) * sqrt(var / n)) %>%
  group_by(dimension, smoothness, n, psi, estimator) %>%
  summarize(Coverage = mean(cilb <= psi & psi <= ciub),
            Power = mean(!(cilb <= 0 & 0 <= ciub))) %>%
  gather(var, value, Coverage, Power) %>%
  mutate(cilb = value - sqrt(value * (1 - value) / NUM_ITERS),
         ciub = value + sqrt(value * (1 - value) / NUM_ITERS),
         dim_label = paste0("d = ", dimension),
         smooth_label = paste0("s = ", smoothness),
         estimator = case_when(
           estimator == "dcdr" ~ "DCDR known density and smoothness",
           estimator == "scdr" ~ "SCDR",
           T ~ "DCDR unknown density or smoothness"
         )) %>%
  ggplot(aes(x = n, y = value, color = estimator)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = cilb, ymax = ciub)) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  scale_x_log10(breaks = c(100, 300, 1000, 5000)) +
  scale_color_colorblind() +
  facet_grid(var ~  reorder(dim_label, dimension) + reorder(smooth_label, smoothness), scales = "free_y") +
  theme_bw() +
  labs(x = "Sample Size", y = "", color = "Estimator") +
  theme(legend.position = "top")

ggsave(plot = p1, filename = "../figures/holder_coverage_power.png",
       width = 8, height = 5)


##############################
### QQ plot

plot_data <- inf_data %>%
  mutate(approx_normal = (psihat - psi) / sqrt(var / n)) %>%
  arrange(dimension, smoothness, n, psi, estimator, approx_normal)

plot_data$normal_quantiles <- rep(qnorm(ppoints(NUM_ITERS)),
                                  nrow(unique(inf_data[c("dimension", "smoothness", "n", "psi", "estimator")])))

p2 <- plot_data %>%
  filter(!(dimension == 4 & estimator == "unknown")) %>%
  filter(dimension == 1) %>%
  mutate(dim_label = paste0("d = ", dimension),
         smooth_label = paste0("s = ", smoothness),
         n_label = paste0("n = ", n),
         estimator = case_when(
           estimator == "dcdr" ~ "DCDR known density and smoothness",
           estimator == "scdr" ~ "SCDR",
           T ~ "DCDR unknown density or smoothness"
         )) %>% 
  ggplot(aes(x = normal_quantiles, y = approx_normal, color = estimator)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  scale_color_colorblind() +
  facet_grid(reorder(dim_label, dimension) + reorder(smooth_label, smoothness) ~ reorder(n_label, n),
             scales = "free_y") +
  theme_bw() +
  labs(x = "Theoretical N(0, 1) Quantiles",
       y = "Quantiles of Standardized Estimate",
       color = "Estimator") +
  theme(legend.position = "top")


ggsave(plot = p2, filename = "../figures/holder_qq_plots.png",
       width = 12, height = 6)

