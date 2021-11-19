# Configuration ----------------------------------------------------------------
set.seet(1)
datapath <- "./Airline.RData"
config_lambdas <- c(0, 0.001, 0.01, 0.1, 1, 5, 10, 100)
config_d <- 2
config_gammas <- c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 10)
config_k <- 20
config_r <- c(2, 3, 4, 5, 6, 7, 8, 9, 10)

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
kernel_linear <- function(X, X2 = NULL, gamma=NULL, d=NULL) {
  if (is.null(X2)){X2 <- X}
  out <- X2 %*% t(X)
  return(out)
}

kernel_inhomogeneous <- function(X, d, X2 = NULL, gamma=NULL) {
  if (is.null(X2)){ X2 <- X}
  out <- (1 + X2 %*% t(X)) ^ d
  return(out)
}

kernel_rbf <- function (X1, gamma, X2 = NULL, d=NULL) {
 
  if (is.null(X2)) {
    n <- nrow(X1)
    XtX1 <- tcrossprod(X1)
    XX1 <- matrix(1, n) %*% diag(XtX1)
    D <- XX1 - 2 * XtX1 + t(XX1)
  }
  else {
    n <- nrow(X2)
    m <- nrow(X1)
    XX2 <- matrix(apply(X2 ^ 2, 1, sum), n, m)
    XX1 <- matrix(apply(X1 ^ 2, 1, sum), n, m, byrow = TRUE)
    X1X2 <- tcrossprod(X2, X1)
    D <- XX2 - 2 * X1X2 + XX1
  }
  
  k <- exp(-gamma * D)
  return(k)
}

kernel_polyspline <- function (X1, X2=NULL, r=2){
  # 
  
  if (is.null(X2)) {
    n <- nrow(X1)
    XtX1 <- tcrossprod(X1)
    XX1 <- matrix(1, n) %*% diag(XtX1)
    D <- XX1 - 2 * XtX1 + t(XX1)
  }  
  else {
    n <- nrow(X2)
    m <- nrow(X1)
    XX2 <- matrix(apply(X2 ^ 2, 1, sum), n, m)
    XX1 <- matrix(apply(X1 ^ 2, 1, sum), n, m, byrow = TRUE)
    X1X2 <- tcrossprod(X2, X1)
    D <- XX2 - 2 * X1X2 + XX1
  }
  
  if (r%%2==0){
    k = D ^ r
  }
  else {
    k = (D ^ (r-1)) * log(D ^ D)
  }
  
  return(k)
}

config_kernel <- kernel_polyspline

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
    #mean_DM <- res$mean_DM
    #std_DM <- res$std_DM
    
    K_t <- kernel_function(X_t, ...)
    
    X_u_tilde <- X_u
    #X_u_tilde <- X_u -  outer(rep(1, dim(X_u)[1]), mean_DM)
    #X_u_tilde <- X_u_tilde / outer(rep(1, dim(X_u)[1]), std_DM)

    ku <- kernel_function(X_t, X2 = X_u_tilde, ...)
    
    predictions <- w_0[1] + ku %*% ginv(K_t) %*% q_tilde
    
    return(predictions)
  }

# Script -----------------------------------------------------------------------
load(datapath)

df <- Airline

df <- fastDummies::dummy_cols(df, select_columns = c('airline')) %>% dplyr::select(-airline)

# don't we need more than 6 obs OOS?
df_train <- df %>% filter(year < max(year))
df_test <- df %>% filter(year == max(year))

# Preprocessing ----------------------------------------------------------------
# why is this scaling happening like this?
df_test <- df_test %>% mutate(
  year = (year - mean(df_train$year)) / sd(df_train$year),
  cost = (cost - mean(df_train$cost)) / sd(df_train$cost),
  pf = (pf - mean(df_train$pf)) / sd(df_train$pf),
  lf = (lf - mean(df_train$lf)) / sd(df_train$lf)
)

df_train[,c('year', 'cost', 'pf', 'lf')] <- df_train[,c('year', 'cost', 'pf', 'lf')] %>% scale

# Data analysis ----------------------------------------------------------------
t1 <- Airline %>% group_by(year) %>% summarise_all(c(mean, sd))
plot(t1$year, t1$output_fn1)

# Cross Validation -------------------------------------------------------------


run_cv <- function(df_train, df_test, mygamma, mylambda, myr){
  X_train <- dplyr::select(df_train, -output) %>% as.matrix()
  y_train <- df_train$output
  X_test <- dplyr::select(df_test, -output) %>% as.matrix()
  y_test <- df_test$output
  
  # why are there so many variables here that are not defined in the function?
  res <- krr(y_train, X_train, lambda=mylambda, config_kernel, gamma=mygamma, d=config_d)
  predictions <- predict_oos(X_test, X_train, res, config_kernel, gamma=mygamma, d=config_d)
  mse <- mean((predictions - y_test)^2)
  return(mse)
}

cv <- df_train %>% split(.$year)
cv <- tibble(year = c(1:length(names(cv))), test = cv)

filter_by_year <- function(filteryear, df){df %>% filter(df$year != filteryear)}
cv$train <- map(cv$test, ~ filter_by_year(mean(.$year), df_train))

# Gridsearch -------------------------------------------------------------------
crossv_cv <- tibble()
for (config_lambda in config_lambdas){
  for (config_gamma in config_gammas){
    cv_lambda <- cv
    cv_lambda$lambda <- config_lambda
    cv_lambda$gamma <- config_gamma
    cv_lambda$mse <- map2_dbl(cv$train, cv$test, ~ run_cv(.x, .y, config_gamma, config_lambda))
    crossv_cv <- bind_rows(crossv_cv, cv_lambda)
  }
}

crossv_output <- crossv_cv %>% select(c(lambda, gamma, mse)) %>% group_by(lambda, gamma) %>% summarise_all(mean)
best_lambda <- crossv_output %>% filter(mse == min(crossv_output$mse)) %>% pull(lambda)
best_gamma <- crossv_output %>% filter(mse == min(crossv_output$mse)) %>% pull(gamma)

print(best_lambda)
print(best_gamma)
print(min(crossv_output$mse))

# Final Testset ----------------------------------------------------------------
X_train <- dplyr::select(df_train, -output) %>% as.matrix()
y_train <- df_train$output
X_test <- dplyr::select(df_test, -output) %>% as.matrix()
y_test <- df_test$output

res <- krr(y_train, X_train, lambda = best_lambda, config_kernel, gamma = best_gamma, d=config_d)
predictions <- predict_oos(X_test, X_train, res, config_kernel, gamma = best_gamma, d=config_d)
mse <- mean((predictions - y_test)^2)
print(mse)

# Comparison -------------------------------------------------------------------
## LM --------------------------------------------------------------------------
lm_model <- lm(y_train ~ ., data=as.data.frame(cbind(X_train, y_train)))
lm_preds <- predict(lm_model, data=as.data.frame(cbind(X_test, y_test)))
mse <- mean((lm_preds - y_test)^2)
print(mse) # Higher than our model

## DSMLE -----------------------------------------------------------------------
res_pkg <-
  dsmle::krr(y_train, X_train, config_k, kernel.type = "nonhompolynom", kernel.degree = 1)
y_hat_dsmle <- res_pkg$yhat
mse <- mean((y_hat_dsmle - y_train)^2)
print(mse) # Higher than our model
