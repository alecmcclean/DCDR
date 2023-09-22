# DCDR

Code for McClean et al. 2023 ``Double Cross-fit Doubly Robust Estimators: Beyond Series Regression''

The code conducts two analyses:
    1. It generates data from Doppler nuisance functions and estimate the Expected Conditional Covariance
    with k-Nearest Neighbors.  This analysis shows that with undersmoothing (i.e., small k) will be optimal
    with double cross-fitting (i.e., independent training datasets).

    2. It generates Holder smooth functions and esitmated the Expected Conditional Covariance with local 
    polynomial regression and covariate-density-adapted local polynomial regression.  The figures
    confirm the theoretical results in Theorems 1, 2, and 3 of the paper.
