###############################################################################
### Author: Alec McClean
### Purpose: Functions for functional undersmoothing
###############################################################################

#################################
### Helper functions and kernels

expit <- function(x) {
  return (1 / (1 + exp(-x)))
}

TsybakovKernel <- function(u) {
  ifelse(abs(2 * u) > 1, 0,
         exp ( -1 / (1 - (2 * u)^2)))
}

EpanechnikovKernel <- function(u) {
  ifelse (abs(u) > 1, 0,
          0.75 * (1 - u^2))
}

BoxcarKernel <- function(u) {
  return(0.5 * (abs(u) < 1))
}


#####################################################
### Functions for 1_knn_mse.R
#####################################################

DGP1 <- function(N, PSI) {
  # Inputs:
  #   N   - size of dataset
  #   PSI - Value of the true ECC
  
  # Create covariates
  dat <- data.frame(X = runif(N))
  
  # Define E(A | X) and A
  dat$pi <- sqrt( dat$X * (1 - dat$X) ) * sin( 2.1 * base::pi / (dat$X + 0.05))
  
  # Define E(Y | X) = E(A | X)
  dat$mu <- dat$pi
  
  # Add noise
  if (PSI > 0) {
    
    # If psi > 0, then give nuisance functions the same error
    dat$A <- dat$pi + rnorm(n = N, mean = 0, sd = sqrt(PSI))
    dat$Y <- dat$A
    
  } else {
    
    # Make independent if PSI = 0
    dat$A <- dat$pi + rnorm(n = N, mean = 0, sd = 0.05)
    dat$Y <- dat$mu + rnorm(n = N, mean = 0, sd = 0.05)
    
  }
  
  # Return dataset
  return(dat)
}

CalculateModels_kNN <- function(data, K) {
  # Input:
  # data  - dataset with X, A, Y
  # K    - number of neighbors
  
  # Randomly Split the dataset into three equal parts
  shuffled_vec <- sample(1:nrow(data))
  split_vec <- split(shuffled_vec, rep(1:3, each = nrow(data) / 3))
  
  train_A <- data[split_vec[[1]],]
  train_Y <- data[split_vec[[2]],]
  estimation <- data[split_vec[[3]],]
  
  # Estimate E(A | X)
  pimod <- knn.reg(train = as.matrix(train_A$X),
                   y = train_A$A,
                   test = as.matrix(estimation$X),
                   k = K)
  
  estimation$pihat <- pimod$pred
  
  # E(Y | X)
  mumod <- knn.reg(train = as.matrix(train_Y$X),
                   y = train_Y$Y,
                   test = as.matrix(estimation$X),
                   k = K)
  
  estimation$muhat <- mumod$pred
  
  # E(Y|X) trained on same fold as E(A|X) for SCDR estimator
  # So prediction will be the same.
  estimation$muhat_scdr <- estimation$pihat
  
  return(estimation)  
}


GetCCEstimate_kNN <- function(FOLD_SIZE, K, PSI) {
  # Inputs:
  # N - size of each fold in the data
  # K - Number of neighbors to use for nuisance function estimators
  # PSI - value of functional
  
  data <- DGP1(3 * FOLD_SIZE, PSI) 
  data <- CalculateModels_kNN(data, K)
  AMSE <- mean((data$A - data$pihat)^2)
  YMSE <- mean((data$Y - data$muhat)^2)
  psihat_dcdr <- mean((data$Y - data$muhat) * (data$A - data$pihat))
  
  # Estimate if both nuisance functions trained on same fold
  psihat_scdr <- mean((data$Y - data$muhat_scdr) * (data$A - data$pihat))
  
  psiMSE_dcdr <- (psihat_dcdr - PSI)^2
  psiMSE_scdr <- (psihat_scdr - PSI)^2
  
  return(list(AMSE, YMSE, psiMSE_dcdr, psiMSE_scdr))
}


################################################################
### Functions for 2_holder.R
################################################################

##############################################
### bump()
### Adds a bump at a given grid point x0 
### See Tsybakov pg. 92

bump <- function(x, # Covariates 
                 n, # Number of samples
                 beta, # Smoothness of function
                 x0, # Grid point for a bump
                 K,  # Kernel for bump
                 L # Constant controlling signal to noise ratio
) {
  
  h <- n^(-1 / (2 * beta + 1)) * 10 # "bandwidth value", implicitly setting c_0
  u <- (x - x0) / h # kernel value
  Kh <- K(u) * (h^beta)
  C <- ifelse(sum(Kh) == 0, 0, L / sum(Kh)) # L controls signal
  C <- C * 2 * (rbinom(1, 1, 0.5) - 0.5) # Randomly multiply by -1 or 1
  if (any(is.infinite(C))) {
    stop("Vector contains missing values (NA)")
  }
  return(C * Kh)
}

#########################################################
### Holder()
### Constructs a Holder smooth function

Holder <- function(x, # Covariates
                   beta, # Smoothness
                   K, # Kernel for bump function
                   L # Parameter controlling signal to noise ratio
) {
  
  # Create lots of bumps, evenly spaced across X
  grid <- seq(0, 1, length.out = 0.5 * length(x)^(1 / (2 * beta + 1))) # n^(1 / (2s + d)) bumps
  y <- rep(0, length(x))
  for (X0 in grid) {
    y <- y + bump(x, n = length(x), beta = beta, x0 = X0, K = K, L = L)
  }
  
  return(y)
}

#############################
### SimulateData()
### Simulated data according to specific value of psi = ECC
### and smoothness level of nuisance functions

SimulateData <- function(psi, # ECC
                         dimension, # dimension of data
                         beta, # Nuisance function smoothness
                         N, # Sample size
                         signal_amp # User parameter controlling signal to noise ratio
) {
  
  # Create X data 
  dat <- data.frame(sapply(1:dimension, function(i) runif(N)))
  if(dimension == 1) {
    colnames(dat) <- c("X1")
  }
  
  # Generate Holder smooth function in each dimension; sum is also Holder
  dat$mu <- rep(0, N)
  for (d in 1:dimension) {
    dat$mu <- dat$mu + Holder(dat[, d],
                              beta,
                              K = TsybakovKernel,
                              L = signal_amp * psi)  
  }
  
  if (psi > 0) { # A = Y implies E(cov) = variance of noise
    dat$Y <- dat$mu + rnorm(N, mean = 0, sd = sqrt(psi))
    dat$A <- dat$Y
  } else { # A _||_ Y implies E(cov) = 0.
    dat$Y <- dat$mu + rnorm(N, mean = 0, sd = 1)
    dat$A <- dat$mu + rnorm(N, mean = 0, sd = 1)
  }
  
  return(dat) # Return dataset
}

### The covariate-density-adapted LPR using Epanechnikov kernel and assuming
### X ~ Unif[0, 1] so f(X) = 1 
CDA_LPR <- function(est_cov, bandwidth, train_cov, train_outcome) {
  diffs <- sweep(as.matrix(train_cov), 2, as.matrix(est_cov), "-")
  norms <- apply(diffs, MARGIN = 1, FUN = function(row) sqrt(sum(row^2)))
  u <- norms / bandwidth
  return(mean( bandwidth^(-1 * ncol(train_cov)) * EpanechnikovKernel(u) * train_outcome))
}

### A simple LPR estimator with Boxcar kernel
### Assumes d = 1
LPR <- function(est_cov, train_cov, train_outcome) {
  
  ### Use the 10 nearest neighbors - implicitly choosing very small bandwidth
  ### kernel and local linear regression
  knn_index <- get.knnx(train_cov, est_cov, k = 10)$nn.index
  
  est_outcome <- c()
  for (k in 1:length(est_cov)) {
    train_cov_close <- train_cov[knn_index[k,]]
    train_outcome_close <- train_outcome[knn_index[k,]]
    mod <- lm(train_outcome_close ~ train_cov_close)
    est_outcome <- c(est_outcome, coef(mod)[[1]] + coef(mod)[[2]] * est_cov[k])
  }
  
  return(est_outcome)
}


###############################################
### GenerateEstimateHolderSmooth()
### Create Holder smooth functions and estimate
### ECC using various local polynomial regressions

GenerateEstimateHolderSmooth <- function(FOLD_SIZE, DIMENSION, BETA, PSI, ITER_NUMBER) {
  
  library(tidyverse)
  library(magrittr)
  library(FNN)
  
  ### 1) Create dataset of total size 3 * FOLD_SIZE
  data <- SimulateData(
    psi = PSI, 
    dimension = DIMENSION,
    beta = BETA, 
    N = 3 * FOLD_SIZE, 
    signal_amp = 50
  )
  
  # Shuffle row indices
  indices <- sample(nrow(data))
  
  # Create a fold index: 1,2,3 repeated FOLD_SIZE times
  # Make sure length matches exactly 3 * FOLD_SIZE
  fold_index <- rep(0:2, each = FOLD_SIZE)
  # Shuffle fold_index to match the randomly permuted `indices`
  fold_index <- fold_index[order(indices)]
  
  # Reorder `data` rows so that row 1 has fold_index[1], row 2 has fold_index[2], etc.
  data <- data[order(indices), ]
  
  ### 2) Define bandwidths (they are the same across all folds)
  if (BETA < DIMENSION / 4) {
    # Undersmooth pihat and make muhat consistent
    PI_BANDWIDTH <- FOLD_SIZE^(-2.01 / (4 * BETA + DIMENSION))
    MU_BANDWIDTH <- FOLD_SIZE^(-1 / DIMENSION)
  } else {
    # Undersmooth both with bandwidth like n^{-1/d}
    PI_BANDWIDTH <- MU_BANDWIDTH <- FOLD_SIZE^(-1 / DIMENSION)
  }
  
  # Construct optimal bandwidth for SCDR estimator
  OPTIMAL_BANDWIDTH <- FOLD_SIZE^(-1 / (2 * BETA + DIMENSION))
  
  ### 3) Initialize empty columns for all predictions in the full data
  data <- data %>%
    mutate(
      pihat          = NA_real_,
      muhat          = NA_real_,
      pihat_scdr     = NA_real_,
      muhat_scdr     = NA_real_,
      pihat_unknown  = NA_real_,
      muhat_unknown  = NA_real_
    )
  
  ### 4) Loop over each of the 3 folds
  for (i in 0:2) {

    # Define test/estimation indices as the rows in fold i
    test_indices  <- which(fold_index == i)
    
    # Define training indices
    train_A_indices <- which(fold_index == mod(i+1,3))
    train_Y_indices <- which(fold_index == mod(i+2,3))
    
    # Separate training sets for double cross-fitting
    train_A <- data[train_A_indices, ]
    train_Y <- data[train_Y_indices, ]
    
    # Same training set for single cross-fitting
    train_both <- bind_rows(train_A, train_Y)
    
    # Test (estimation) set
    test_set <- data[test_indices, ]
    stopifnot("Missing data" = nrow(train_A) > 0 & nrow(train_Y) > 0 & nrow(test_set) > 0)
    
    ### 4a) Compute (pihat, muhat) for each test point in this fold
    for (k in seq_len(nrow(test_set))) {
      
      # Covariate of the test point
      x_test <- select(test_set, contains("X")) %>% slice(k)
      
      # ========== DCDR known density and smoothness ========== #
      # Undersmoothed pihat
      pihat_val <-
        CDA_LPR(
          x_test,
          bandwidth     = PI_BANDWIDTH,
          train_cov     = select(train_A, contains("X")),
          train_outcome = train_A$A
        )
      
      # Consistent muhat
      muhat_val <-
        CDA_LPR(
          x_test,
          bandwidth     = MU_BANDWIDTH,
          train_cov     = select(train_Y, contains("X")),
          train_outcome = train_Y$Y
        )
      
      # ========== SCDR-MSE (optimal bandwidth) ========== #
      pihat_scdr_val <-
        CDA_LPR(
          x_test,
          bandwidth     = OPTIMAL_BANDWIDTH,
          train_cov     = select(train_both, contains("X")),
          train_outcome = train_both$A
        )
      
      muhat_scdr_val <- pihat_scdr_val
      
      # Save in test_set
      test_set$pihat[k]       <- pihat_val
      test_set$muhat[k]       <- muhat_val
      test_set$pihat_scdr[k]  <- pihat_scdr_val
      test_set$muhat_scdr[k]  <- muhat_scdr_val
    }
    
    ### 4b) Re-insert test_set predictions back into the main data
    data[test_indices, c("pihat","muhat","pihat_scdr","muhat_scdr")] <-
      test_set[, c("pihat","muhat","pihat_scdr","muhat_scdr")]
    
    if (isTRUE(DIMENSION == 1) && isTRUE(BETA <= 1)) {

      if (length(test_indices) > 0 && length(train_A$X1) > 0) {
        
        data[test_indices, c("pihat_unknown")] <-
          LPR(
            est_cov       = data[test_indices,]$X1,
            train_cov     = train_A$X1,
            train_outcome = train_A$A
          )
        
        data[test_indices, c("muhat_unknown")] <-
          LPR(
            est_cov       = data[test_indices,]$X1,
            train_cov     = train_Y$X1,
            train_outcome = train_Y$Y
        )
    
      
      } else {
        data[test_indices, c("pihat_unknown")] <- 0
        data[test_indices, c("muhat_unknown")] <- 0
      }
    }
  }
  

  
  ### 6) Now compute the influence-function–based estimates
  #     over the *entire* data set (all folds combined).
  data <- data %>%
    mutate(
      varphihat_dcdr    = (A - pihat) * (Y - muhat),
      varphihat_scdr    = (A - pihat_scdr) * (Y - muhat_scdr),
      varphihat_unknown = (A - pihat_unknown) * (Y - muhat_unknown)
    )
  
  # Safety check: if everything predicted is 0
  if (all(data$pihat == 0, data$muhat == 0)) {
    stop("All predictions (pihat, muhat) were zero.")
  }
  
  ### 7) Summarize into a final data.frame to return
  dat <- data.frame(
    estimator = "dcdr",
    psihat    = mean(data$varphihat_dcdr),
    var       = var(data$varphihat_dcdr)
  ) %>%
    bind_rows(
      data.frame(
        estimator = "scdr",
        psihat    = mean(data$varphihat_scdr),
        var       = var(data$varphihat_scdr)
      )
    ) %>%
    bind_rows(
      data.frame(
        estimator = "unknown",
        psihat    = mean(data$varphihat_unknown),
        var       = var(data$varphihat_unknown)
      )
    ) %>%
    mutate(
      dimension  = DIMENSION,
      smoothness = BETA, 
      n          = FOLD_SIZE * 3,
      psi        = PSI,
      iter       = ITER_NUMBER
    )
  
  return(dat)
}
  