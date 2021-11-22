# Packages ---------------------------------------------------------------------
#install.packages(c("SVMMaj", "ISLR", "plotrix"))
library(pacman)
p_load(MASS,mlbench, Hmisc,
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
## ---------------------------------KERNELS-------------------------------------
#-------------------------------------------------------------------------------
# kernel_linear: constructs the kernel matrix based in the linear kernel. 
#                Can be applied to training data or to construct ku and allow 
#                out of sample predictions.
#
# Parameters:
#     X:  matrix, contains independent variables
#     X2: matrix, contains out-of-sample data, used in the construction of ku. If
#                 NULL, X is used and k is obtained.
# 
# Returns:
#     out: matrix, the kernel matrix
kernel_linear <- function(X, hyperpar = NULL, X2 = NULL, ...) {
  if (is.null(X2)){X2 <- X}
  out <- X2 %*% t(X)
  return(out)
}

#-------------------------------------------------------------------------------
# kernel_inhomogeneous: constructs the kernel matrix based in the inhomogeneous
#                       kernel. 
#                       Can be applied to training data or to construct ku and 
#                       allow out of sample predictions.
#
# Parameters:
#     X:  matrix, contains independent variables
#     X2: matrix, contains out of sample data, used in the construction of ku. If
#                 NULL, X is used and k is obtained.
#     d:  float, the degree of the transformation
# 
# Returns:
#     out: matrix, the kernel matrix
kernel_inhomogeneous <- function(X, hyperpar, X2 = NULL, ...) {
  if (is.null(X2)){ X2 <- X}
  out <- (1 + X2 %*% t(X)) ^ hyperpar
  return(out)
}

#-------------------------------------------------------------------------------
# kernel_rbf: constructs the kernel matrix based in the rbf kernel. 
#                Can be applied to training data or to construct ku and allow 
#                out of sample predictions.
#
# Parameters:
#     X:  matrix, contains independent variables
#     X2: matrix, contains out of sample data, used in the construction of ku. If
#                 NULL, X is used and k is obtained.
#     gamma: float, hyperparameter that penalized the distance among the x's
# 
# Returns:
#     out: matrix, the kernel matrix
kernel_rbf <- function (X1, hyperpar, X2 = NULL, ...) {
  
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
  
  k <- exp(-hyperpar * D)
  return(k)
}

#-------------------------------------------------------------------------------
# kernel_polyspline : constructs the kernel matrix based in the polyspline kernel. 
#                Can be applied to training data or to construct ku and allow 
#                out of sample predictions.
#
# Parameters:
#     X:  matrix, contains independent variables
#     X2: matrix, contains out of sample data, used in the construction of ku. If
#                 NULL, X is used and k is obtained.
#     r: float, degree of the polyspline.
# 
# Returns:
#     out: matrix, the kernel matrix
kernel_polyspline <- function (X1, hyperpar, X2=NULL, ...){
  
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
  
  if (hyperpar%%2==0){
    k = D ^ hyperpar
  }
  else {
    k = (D ^ (hyperpar-1)) * log(D ^ D)
  }
  
  return(k)
}

#-------------------------------------------------------------------------------
## Ridge Regression ------------------------------------------------------------
#krr: function fits the Ridge Regression. Supports kernels.
#
#
# Parameters:
#     y:  vector, labels for the observations
#     X:  matrix, contains the training observations to fit the model
#     lambda: float, penalization
#     kernel_function: the function to produce the kernel from the X matrix
# 
# Returns:
#     preds:    vector, predictions
#     q_tilde:  vector, equivalent to predictions without the intercept
#     w_0:      float, intercept
#     mean_DM:  vector, design matrix column means
#     std_DM:   vector, design matrix column standard deviations
krr <- function(y, X, kernel_function, hyperpar, lambda) {
  n <- nrow(X)
  ones <- matrix(1, n, 1)
  mean_DM <-
    colMeans(X)#Storages mean of the columns for the training set design matrix
  std_DM <-
    apply(X, 2, sd)#Storages the std of the columns for the training set design matrix
  I <- diag(n)
  
  J <- I - (1 / n) * (ones %*% t(ones))
  X_tilde <- J %*% X #centralized X
  k <- kernel_function(X, hyperpar=hyperpar)
  k <- J %*% k %*% J #centralized K
  
  w_0 <- (1 / n) * t(ones) %*% y
  q_tilde <-solve(I + lambda * ginv(k), J %*% y)
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

#-------------------------------------------------------------------------------
# predict_oos: Used the estimated model 
#
# Parameters:
#     X_u:  matrix, independent variables, out-of-sample
#     X_t:  matrix, independent variables, training data
#     kernel_function: the function to produce the kernel from the X matrix
#     ...:    extra parameters
#
# Returns:
#     predictions: vector, out-of-sample predictions
predict_oos <-
  function(X_u,
           X_t,
           res,
           kernel_function,
           hyperpar) {
    w_0 <- res$w_0
    q_tilde <- res$q_tilde
    
    K_t <- kernel_function(X_t, hyperpar=hyperpar)
    
    X_u_tilde <- X_u
    
    ku <- kernel_function(X_t, hyperpar=hyperpar, X2 = X_u_tilde)
    
    predictions <- w_0[1] + ku %*% ginv(K_t) %*% q_tilde
    
    return(predictions)
  }

# Cross Validation Functions ---------------------------------------------------
# Wrapped helper functions to facilitate the gridsearch implementation
run_cv <- function(df_cv_train, df_cv_test, my_hyperpar, kernel, lambda){
  X_train <- dplyr::select(df_cv_train, -output) %>% as.matrix()
  y_train <- df_cv_train$output
  X_test <- dplyr::select(df_cv_test, -output) %>% as.matrix()
  y_test <- df_cv_test$output
  
  res <- krr(y_train, X_train, kernel, my_hyperpar, lambda)
  predictions <- predict_oos(X_test, X_train, res, kernel, my_hyperpar)
  mse <- mean((predictions - y_test)^2)
  return(mse)
}


compare_params <- function(cv, hyperpars, my_kernel, lambdas, verbose = FALSE){
  crossv_cv <- tibble()
  for (hyperpar in hyperpars){
    for (my_lambda in lambdas){
      if(verbose){
        print('---------------------------------------------------------------')
        print('Running CV for params:')
        print(paste('Hyperpar:', hyperpar))
        print(paste('Lambda:', my_lambda))    
      }
      
      cv_run <- cv
      cv_run$hyperparameter <- hyperpar
      cv_run$lambda <- my_lambda
      cv_run$mse <- map2_dbl(cv$train, cv$test, ~ run_cv(.x, .y, hyperpar, my_kernel, my_lambda))
      crossv_cv <- bind_rows(crossv_cv, cv_run)
    }
  }
  return(crossv_cv)
}

compare_kernels <- function(cv, lambdas, my_verbose=FALSE){
  run_cv <- function(kernel, hyperpars){compare_params(cv, hyperpars, kernel, config_lambdas, verbose=my_verbose) %>% bind_rows(model_comparisons)}
  model_comparisons <- tibble()
  
  print('Running RBF')
  model1 <- run_cv(kernel_rbf, config_gammas)
  model1['kernel'] <- 'RBF'
  
  print('Running Linear')
  model2 <- run_cv(kernel_linear, NULL)
  model2['kernel'] <- 'linear'
  
  print('Running inhomogeneous')
  model3 <- run_cv(kernel_inhomogeneous, config_ds)
  model3['kernel'] <- 'inhomogeneous'
  
  print('Running polyspline')
  model4 <- run_cv(kernel_polyspline, config_rs)
  model4['kernel'] <- 'polyspline'
  
  model_comparisons <- model_comparisons %>% bind_rows(model1) %>% bind_rows(model2) %>% bind_rows(model3) %>% bind_rows(model4)
  
  return(model_comparisons)
}
