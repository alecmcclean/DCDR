###############################################################################
### Author: Alec McClean
### Purpose: Master script for functional undersmoothing.  Loads packages
###          and calls other scripts.
###############################################################################

if (!require(pacman)) {
  install.packages("pacman")
  library(pacman)
}

p_load(assert,
       doParallel,
       doSNOW,
       FNN,
       foreach,
       ggthemes,
       latex2exp,
       locfit,
       KernSmooth,
       magrittr,
       MASS,
       np,
       progress,
       tidyverse)

options(stringsAsFactors = F)

source("1_knn_mse.R")
source("2_holder.R")