################################################################################
#-------------------------------------------------------------------------------
####Purpose: Configuration file, control the script variables.
# Configuration ----------------------------------------------------------------
set.seed(1)
datapath <- "./Airline.RData"
config_lambdas <- c(1, 5, 10, 0.001, 0.01, 0.1, 0.2, 0.5, 100)
config_ds <- c(2, 4)
config_gammas <- c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 10)
config_k <- 20
config_rs <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
config_years_to_predict_in_testset <- 3
