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
kernel_linear <- function(X, X2 = NULL, ...) {
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

kernel_inhomogeneous <- function(X, d, X2 = NULL, ...) {
  if (is.null(X2)){ X2 <- X}
  out <- (1 + X2 %*% t(X)) ^ d
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

kernel_rbf <- function (X1, gamma, X2 = NULL, ...) {
  
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

kernel_polyspline <- function (X1, r=2, X2=NULL, ...){
  
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

krr <- function(y, X, lambda, kernel_function, ...) {
  n <- nrow(X)
  ones <- matrix(1, n, 1)
  mean_DM <-
    colMeans(X) #Storages mean of the columns for the training set design matrix
  std_DM <-
    apply(X, 2, sd) #Storages the std of the columns for the training set design matrix
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

# Cross Validation Functions ---------------------------------------------------
# Wrapped helper functions to facilitate the gridsearch implementation
run_cv <- function(df_cv_train, df_cv_test, mygamma, mylambda, myr, myd){
  X_train <- dplyr::select(df_cv_train, -output) %>% as.matrix()
  y_train <- df_cv_train$output
  X_test <- dplyr::select(df_cv_test, -output) %>% as.matrix()
  y_test <- df_cv_test$output
  
  res <- krr(y_train, X_train, lambda=mylambda, config_kernel, gamma=mygamma, d=myd, r=myr)
  predictions <- predict_oos(X_test, X_train, res, config_kernel, gamma=mygamma, d=myd, r=myr)
  mse <- mean((predictions - y_test)^2)
  return(mse)
}