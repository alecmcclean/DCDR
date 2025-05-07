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
  data <- data %>% bind_rows(
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
  geom_point(aes(y = Y), color = "grey", alpha = 0.6) +
  geom_line(aes(y = mu)) +
  theme_bw() + 
  facet_grid(reorder(N_label, N) ~ reorder(smooth_label, smoothness), scales = "free_y")

ggsave(plot = p0, filename = "../figures/example_holder.png",
       width = 8, height = 6)

#################################################
### Analyze several Holder smooth functions
#################################################

### Construct datasets and run estimators
FOLD_SIZES <- c(100, 200, 350, 700, 1500, 3000)
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
saveRDS(inf_data, "../output/smoothness-sims-data.RDS")

#################################################
### Create figures to illustrate results 
#################################################

##############################
### Coverage

plot_data <- inf_data %>%
  filter(!(dimension == 4 & estimator == "unknown")) %>%
  mutate(cilb = psihat - qnorm(0.975) * sqrt(var / n),
         ciub = psihat + qnorm(0.975) * sqrt(var / n)) %>%
  group_by(dimension, smoothness, n, psi, estimator) %>%
  summarize(Coverage = mean(cilb <= psi & psi <= ciub),
            emp_var = mean((psihat-mean(psihat))^2),
            emp_var_var = var((psihat - mean(psihat))^2),
            emp_bias2 = (mean(psihat - psi))^2,
            emp_bias2_var = var(psihat - psi) * (2 * mean(psihat-psi))^2,
            emp_mse = mean((psihat-psi)^2),
            emp_mse_var = var((psihat - psi)^2)) %>%
  ungroup() 

plot_data <- bind_rows(
  plot_data %>% 
    mutate(
      var = "Coverage",
      cilb = Coverage - qnorm(0.975) * sqrt(Coverage * (1 - Coverage) / NUM_ITERS),
      ciub = Coverage + qnorm(0.975) * sqrt(Coverage * (1 - Coverage) / NUM_ITERS)
    ) %>%
    select(dimension:estimator, var, value = Coverage, cilb, ciub),
  plot_data %>% 
    mutate(
      var = "Monte carlo empirical variance",
      cilb = emp_var - qnorm(0.975) * sqrt(emp_var_var / NUM_ITERS),
      ciub = emp_var + qnorm(0.975) * sqrt(emp_var_var / NUM_ITERS)
    ) %>%
    select(dimension:estimator, var, value = emp_var, cilb, ciub),
  plot_data %>% 
    mutate(
      var = "Monte carlo empirical bias squared",
      cilb = pmax(emp_bias2 - qnorm(0.975) * sqrt(emp_bias2_var / NUM_ITERS), 0),
      ciub = emp_bias2 + qnorm(0.975) * sqrt(emp_bias2_var / NUM_ITERS)
    ) %>%
    select(dimension:estimator, var, value = emp_bias2, cilb, ciub),
  plot_data %>%
    mutate(
      var = "Monte carlo empirical MSE",
      cilb = pmax(emp_mse - qnorm(0.975) * sqrt(emp_mse_var / NUM_ITERS), 0),
      ciub = emp_mse + qnorm(0.975) * sqrt(emp_mse_var / NUM_ITERS)
    ) %>% 
    select(dimension:estimator, var, value = emp_mse, cilb, ciub),
  )

plot_data %<>% 
  mutate(
    dim_label = paste0("d = ", dimension),
    smooth_label = paste0("s = ", smoothness),
    estimator = case_when(
      estimator == "dcdr" ~ "DCDR known density and smoothness",
      estimator == "scdr" ~ "SCDR-MSE",
      T ~ "DCDR undersmoothed LPR"
      )
  ) 

p1 <- plot_data %>%
  filter(var == "Coverage") %>%
  ggplot(aes(
    x = n, 
    y = value, 
    shape = estimator, 
    color = estimator
  )) +
  geom_point(size = 2) +  # Increase point size here
  geom_line() +
  geom_errorbar(aes(ymin = cilb, ymax = ciub), alpha = 0.4) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  scale_x_log10(breaks = c(300, 1000, 3000, 9000)) +
  scale_color_colorblind() +
  scale_shape_manual(
    values = c(
      "DCDR known density and smoothness" = 16,  # solid square
      "DCDR undersmoothed LPR" = 15, # solid circle
      "SCDR-MSE" = 17                           # solid triangle
    )
  ) +
  facet_wrap(
    ~ reorder(dim_label, dimension) + reorder(smooth_label, smoothness),
    scales = "free"
  ) +
  theme_bw() +
  labs(x = "Sample Size", y = "", color = "Estimator", shape = "Estimator") +
  theme(legend.position = "top")

ggsave(plot = p1, filename = "../figures/holder_ci_coverage.png",
       width = 8, height = 6)


######################################
### QQ plots

qq_data <- inf_data %>%
  mutate(approx_normal = (psihat - psi) / sqrt(var / n)) %>%
  arrange(dimension, smoothness, n, psi, estimator, approx_normal)

qq_data$normal_quantiles <- 
  rep(qnorm(ppoints(NUM_ITERS)),
      nrow(unique(inf_data[c("dimension", "smoothness", "n", "psi", "estimator")])))

p2a <- qq_data %>%
  filter(!(dimension == 4 & estimator == "unknown")) %>%
  filter(dimension == 1) %>%
  filter(n %in% c(300, 1050, 4500, 9000)) %>%
  mutate(dim_label = paste0("d = ", dimension),
         smooth_label = paste0("s = ", smoothness),
         n_label = paste0("n = ", n),
         estimator = case_when(
           estimator == "dcdr" ~ "DCDR known density and smoothness",
           estimator == "scdr" ~ "SCDR-MSE",
           T ~ "DCDR undersmoothed LPR"
         )) %>% 
  ggplot(aes(x = normal_quantiles, y = approx_normal, color = estimator,
             shape = estimator)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  scale_color_colorblind() +
  scale_shape_manual(
    values = c(
      "DCDR known density and smoothness" = 16, 
      "SCDR-MSE" = 17,                             
      "DCDR undersmoothed LPR" = 15
    )
  ) +
  facet_grid(reorder(dim_label, dimension) + reorder(smooth_label, smoothness) ~ reorder(n_label, n),
             scales = "free_y") +
  theme_bw() +
  labs(x = "Theoretical N(0, 1) Quantiles",
       y = "Quantiles of Standardized Estimate",
       color = "Estimator",
       shape = "Estimator") +
  theme(legend.position = "top")

ggsave(plot = p2a, filename = "../figures/holder_qq_plots.png",
       width = 8, height = 5)

p2b <- qq_data %>%
  filter(dimension != 4 & estimator != "dcdr", smoothness == 0.35) %>%
  filter(n %in% c(300, 1050, 4500)) %>%
  mutate(dim_label = paste0("d = ", dimension),
         smooth_label = paste0("s = ", smoothness),
         n_label = paste0("n = ", n),
         estimator = case_when(
           estimator == "dcdr" ~ "DCDR known density and smoothness",
           estimator == "scdr" ~ "SCDR-MSE",
           TRUE ~ "DCDR undersmoothed LPR"
         )) %>% 
  ggplot(aes(x = normal_quantiles, y = approx_normal, color = estimator, shape = estimator)) +  
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  scale_color_manual(
    values = c(
      "DCDR known density and smoothness" = "#E69F00",  # Default orange
      "SCDR-MSE" = "#56B4E9",  # Light blue
      "DCDR undersmoothed LPR" = "#E69F00"  
    )
  ) +
  scale_shape_manual(
    values = c(
      "DCDR known density and smoothness" = 16, 
      "SCDR-MSE" = 17,                             
      "DCDR undersmoothed LPR" = 15
    )
  ) +
  facet_grid(~ reorder(n_label, n),
             scales = "free_y") +
  theme_bw() +
  labs(x = "Theoretical N(0, 1) Quantiles",
       y = "Quantiles of Standardized Estimate",
       color = "Estimator",
       shape = "Estimator") +
  theme(legend.position = "top")

ggsave(plot = p2b, filename = "../figures/holder_qq_plots_intro.png",
       width = 6, height = 3)

####################################################
### Examine bias, variance, and MSE when d=1

p3 <- plot_data %>%
  filter(dimension == 1, grepl("Monte", var)) %>%
  mutate(var = gsub("Monte carlo e", "E", var),
         var = factor(var, levels = c("Empirical bias squared", "Empirical variance", "Empirical MSE"),
                      labels = c("Squared bias", "Variance", "MSE"))) %>%
  ggplot(aes(x = n, y = value, color = estimator, shape = estimator)) +
  geom_point(size = 2.5) +
  geom_line() +
  geom_errorbar(aes(ymin = cilb, ymax = ciub), alpha = 0.6) +
  scale_x_log10(breaks = c(300, 1000, 3000, 9000)) +
  scale_y_log10() +
  scale_color_colorblind() +
  scale_shape_manual(
    values = c(
      "DCDR known density and smoothness" = 16,  # circle
      "DCDR undersmoothed LPR" = 15, # square
      "SCDR-MSE" = 17 # triangle
    )
  ) +
  facet_grid(reorder(dim_label, dimension) + reorder(smooth_label, smoothness) ~ var, scales = "free") +
  theme_bw() +
  labs(x = "Sample Size", y = "Estimator error", color = "Estimator", shape = "Estimator") +
  theme(legend.position = "top")

p3 <- plot_data %>%
  filter(dimension == 1, grepl("Monte", var)) %>%
  mutate(var = gsub("Monte carlo e", "E", var),
         var = factor(var, 
                      levels = c("Empirical bias squared", "Empirical variance", "Empirical MSE"),
                      labels = c("Squared bias", "Variance", "MSE"))) %>%
  ggplot(aes(x = n, y = value, color = estimator, shape = estimator)) +
  geom_point(size = 2.5) +
  geom_line() +
  
  # Always draw vertical line from the point to the upper CI bound
  geom_segment(aes(x = n, xend = n, y = value, yend = ciub), alpha = 1) +
  # Always draw the upper horizontal cap at the upper bound
  geom_segment(aes(x = n - n/5, xend = n + n/5, y = ciub, yend = ciub)) +
  # Onl draw vertical line from the point to the lower CI bound if cilb nonzero
  geom_segment(aes(x = n, xend = n, y = value, yend = ifelse(cilb == 0, value, cilb))) +
  # Only draw the lower horizontal cap if cilb is nonzero
  geom_segment(data = . %>% filter(cilb != 0),
               aes(x = n - n/5, xend = n + n/5, y = cilb, yend = cilb),
               alpha = 1) +
  scale_x_log10(breaks = c(300, 1000, 3000, 9000)) +
  scale_y_log10() +
  scale_color_colorblind() +
  scale_shape_manual(values = c(
    "DCDR known density and smoothness" = 16,  # solid square
    "DCDR undersmoothed LPR" = 15,             # solid circle
    "SCDR-MSE" = 17                           # solid triangle
  )) +
  facet_grid(reorder(dim_label, dimension) + reorder(smooth_label, smoothness) ~ var,
             scales = "free") +
  theme_bw() +
  labs(x = "Sample Size", y = "Estimator error", color = "Estimator", shape = "Estimator") +
  theme(legend.position = "top")


ggsave(plot = p3, filename = "../figures/holder_bias_var_mse.png",
       width = 8, height = 6)

### Check width for DCDR in d = 1, s = 0.1
extrap_dat <- inf_data %>% 
  filter(estimator == "dcdr", dimension == 1, smoothness == 0.1) %>%
  mutate(ci_width = 2 * qnorm(0.975) * sqrt(var / n))

mod <- lm(log(ci_width) ~ log(n), extrap_dat)
cat("Confidence interval will be smaller than 10 at roughly n = ", 
    signif(exp((coef(mod)[[1]] - log(10)) / abs(coef(mod)[[2]])), 2))

