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
## Kernels training ------------------------------------------------------------
kernel_linear <- function(X, X2 = NULL, ...) {
  if (is.null(X2)){X2 <- X}
  out <- X2 %*% t(X)
  return(out)
}

kernel_inhomogeneous <- function(X, d, X2 = NULL, ...) {
  if (is.null(X2)){ X2 <- X}
  out <- (1 + X2 %*% t(X)) ^ d
  return(out)
}

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