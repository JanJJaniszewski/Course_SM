# Configuration ----------------------------------------------------------------
set.
datapath <- "./Airline.RData"
config_lambdas <- c(5, 10)
config_d <- 2
config_gamma <- 1 / 2
config_k <- 10

# Packages ---------------------------------------------------------------------
#install.packages(c("SVMMaj", "ISLR", "plotrix"))
library(pacman)
p_load(MASS,mlbench,
       tidyverse, tidymodels, parsnip, mlbench,
       caret, modelr,
       rbenchmark,
       parallel,
       doParallel,
       dsmle,
       rdetools,
       tidyverse,
       ggpubr,
       SVMMaj,
       ISLR,
       plotrix)

# Functions --------------------------------------------------------------------
## Kernels training ------------------------------------------------------------
kernel_linear <- function(X, X2 = NULL) {
  if (is.null(X2)){X2 <- X}
  out <- X2 %*% t(X)
  return(out)
}

kernel_inhomogeneous <- function(X, d, X2 = NULL) {
  if (is.null(X2)){ X2 <- X}
  out <- (1 + X2 %*% t(X)) ^ d
  return(out)
}

kernel_rbf <- function (X1, gamma, X2 = NULL) {
  n <- nrow(X1)
  if (is.null(X2)) {
    XtX1 <- tcrossprod(X1)
    XX1 <- matrix(1, n) %*% diag(XtX1)
    D <- XX1 - 2 * XtX1 + t(XX1)
  }
  else {
    m <- nrow(X2)
    XX1 <- matrix(apply(X1 ^ 2, 1, sum), n, m)
    XX2 <- matrix(apply(X2 ^ 2, 1, sum), n, m, byrow = TRUE)
    X1X2 <- tcrossprod(X1, X2)
    D <- XX1 - 2 * X1X2 + XX2
  }
  
  k <- exp(-D * gamma)
  return(k)
}

## Ridge Regression ------------------------------------------------------------
krr <- function(y, X, lambda, kernel_function, ...) {
  n <- nrow(X)
  ones <- matrix(1, n, 1)
  mean_DM <-
    colMeans(X)#Storages mean of the columns for the training set design matrix
  std_DM <-
    apply(X, 2, sd)#Storages the std of the columns for the training set design matrix
  I <- diag(n)
  
  J <- I - (1 / n) * (ones %*% t(ones))
  X_tilde <- J %*% X #centralized X
  k <- kernel_function(X, ...)
  k <- J %*% k %*% J #centralized K
  
  w_0 <- (1 / n) * t(ones) %*% y
  q_tilde <- solve(I + lambda * ginv(k), J %*% y)
  preds <- w_0[1] + q_tilde
  w <- round(ginv(X_tilde) %*% q_tilde, 2)
  model = list(
    'preds' = preds,
    'q_tilde' = q_tilde,
    'w_0' = w_0,
    'w' = w,
    'mean_DM' = mean_DM,
    'std_DM' = std_DM
  )
  class(model) <- 'krr_model'
  return(model)
}

## Prediction ------------------------------------------------------------------
predict_oos <-
  function(X_u,
           X_t,
           res,
           kernel_function,
           ...) {
    w_0 <- res$w_0
    q_tilde <- res$q_tilde
    mean_DM <- res$mean_DM
    std_DM <- res$std_DM
    
    K_t <- kernel_function(X_t, ...)
    
    X_u_tilde <- X_u -  outer(rep(1, dim(X_u)[1]), mean_DM)
    X_u_tilde <- X_u_tilde / outer(rep(1, dim(X_u)[1]), std_DM)

    ku <- kernel_function(X_t, X2 = X_u_tilde, ...)
    
    predictions <- w_0[1] + ku %*% ginv(K_t) %*% q_tilde
    
    return(predictions)
  }

## Cross validation ------------------------------------------------------------
# TODO

# Script -----------------------------------------------------------------------
load(datapath)

df <- Airline

df <- fastDummies::dummy_cols(df, select_columns = c('airline')) %>% dplyr::select(-airline)

df_train <- df %>% filter(year < 12)
nrow(df_train)
df_test <- df %>% filter(year >= 12)
nrow(df_test)

# Preprocessing ----------------------------------------------------------------
df_test <- df_test %>% mutate(
  year = (year - mean(df_train$year)) / sd(df_train$year),
  cost = (year - mean(df_train$cost)) / sd(df_train$cost),
  pf = (year - mean(df_train$pf)) / sd(df_train$pf),
  lf = (year - mean(df_train$lf)) / sd(df_train$lf)
)

df_train[,c('year', 'cost', 'pf', 'lf')] <- df_train[,c('year', 'cost', 'pf', 'lf')] %>% scale

# Cross Validation -------------------------------------------------------------
filter_by_idx <- function(idx){
  df_train %>% filter(as.numeric(rownames(df_train)) %in% idx)
}

run_cv <- function(df_train, df_test){
  X_train <- dplyr::select(df_train, -output) %>% as.matrix()
  y_train <- df_train$output
  X_test <- dplyr::select(df_test, -output) %>% as.matrix()
  y_test <- df_test$output
  
  res <- krr(y_train, X_train, lambda, kernel_rbf, gamma=config_gamma)
  predictions <- predict_oos(X_test, X_train, res, kernel_rbf, gamma=config_gamma)
  mse <- mean((predictions - y_test)^2)
  return(mse)
}

cv  <- crossv_kfold(df_train, k = 5)
cv$train <- map(cv$train, ~ filter_by_idx(.$idx))
cv$test <- map(cv$test, ~ filter_by_idx(.$idx))

# Gridsearch -------------------------------------------------------------------
crossv_cv <- tibble()
for (lambda in config_lambdas){
  cv_lambda <- cv
  cv_lambda$lambda <- lambda
  cv_lambda$mse <- map2_dbl(cv$train, cv$test, ~ run_cv(.x, .y))
  crossv_cv <- bind_rows(crossv_cv, cv_lambda)
}

crossv_output <- crossv_cv %>% select(c(lambda, mse)) %>% group_by(lambda) %>% summarise_all(mean)
best_lambda <- crossv_output %>% filter(mse == min(crossv_output$mse)) %>% pull(lambda)

# Final Testset ----------------------------------------------------------------
X_train <- dplyr::select(df_train, -output) %>% as.matrix()
y_train <- df_train$output
X_test <- dplyr::select(df_test, -output) %>% as.matrix()
y_test <- df_test$output

res <- krr(y_train, X_train, best_lambda, kernel_rbf, gamma = config_gamma)
predictions <- predict_oos(X_test, X_train, res, kernel_rbf, gamma = config_gamma)
mse <- mean((predictions - y_test))

# Results ----------------------------------------------------------------------
# z <- krr(y, X, lambda, kernel_linear)
# a <- krr(y, X, lambda, kernel_rbf, gamma = 1 / 2)
# b <-
#   krr(y, X, lambda, rbfkernel, sigma = 1) # Comparison to standard kernel function
# 
# res <-  krr(y, X, lambda, kernel_inhomogeneous, d = 2)


# Comparison to DSMLE ----------------------------------------------------------
res_pkg <-
  dsmle::krr(y, X, config_k, kernel.type = "nonhompolynom", kernel.degree = 1)
y_hat_dsmle <- res_pkg$yhat

# # Out of Sample Predictions ------------------------------------------------------------------
# w_0 <- res$w_0
# q_tilde <- res$q_tilde
# 
# predsINH <-
#   predict_oos(w_0[1],
#               X_u,
#               X_t,
#               q_tilde,
#               kernel_inhomogeneous,
#               d = 2)
# predsRBF <-
#   predict_oos(w_0[1], X_u, X_t, q_tilde, "rbf", kernel_rbf, gamma = config_gamma)
# 
# 
# 
# #rdetools
# k <- polykernel(X, 2, Y = NULL)
# k <- rbfkernel(X)
# r <- rde(k,
#          y,
#          est_y = TRUE,
#          dim_rest = 1,
#          regression = TRUE)
# y_hat_rdtools <- r$yh
# ? rde()
