# Configuration ----------------------------------------------------------------
datapath <- "./Airline.RData"
lambda <- 10
d <- 1
gamma <- 1/2

# Packages ---------------------------------------------------------------------
#install.packages(c("SVMMaj", "ISLR", "plotrix"))
library(pacman)
p_load("MASS")
p_load("dsmle")
p_load('rdetools')
p_load("tidyverse")
p_load("ggpubr")
p_load("fastDummies")

# Functions --------------------------------------------------------------------
## Loss Function ---------------------------------------------------------------
# L_ridge <- function(w_0, qtilde, y, lambda, q_tilde){
#   y_w0 <- y-ones%*%w_0
#   JI_qtilde <- J%*%y - q_tilde 
#   xtx <- X%*%t(X)
#   loss <- y_w0%*%t(y_w0) + JI_qtilde%*%t(JI_qtilde) + lambda*q_tilde%*%solve(xtx)%*%q_tilde
# }
## Kernels ---------------------------------------------------------------------
kernel_inhomogeneous <- function(X, d){
  A <- (1 + X%*%t(X))^d
  return(A)
}

kernel_linear <- function(X){
  xxt <- X%*%t(X)
}

kernel_rbf <- function (X1, gamma, X2=NULL) {
  n <- nrow(X1)
  if (is.null(X2)) {
    XtX1 <- tcrossprod(X1)
    XX1 <- matrix(1, n) %*% diag(XtX1)
    D <- XX1 - 2 * XtX1 + t(XX1)
  }
  else {
    m <- nrow(X2)
    XX1 <- matrix(apply(X1^2, 1, sum), n, m)
    XX2 <- matrix(apply(X2^2, 1, sum), n, m, byrow = TRUE)
    X1X2 <- tcrossprod(X1, X2)
    D <- XX1 - 2 * X1X2 + XX2
  }
  
  k <- exp(-D * gamma)
  return(k)
}

## Ridge Regression -------------------------------------------------------------
krr <- function(y, X, lambda, kernel_function, ...){
  n <- nrow(X)
  ones <- matrix(1, n, 1)
  mean_DM <- colMeans(X)#Storages mean of the columns for the training set design matrix
  std_DM <- apply(X, 2, sd)#Storages the std of the columns for the training set design matrix
  I <- diag(n)
  
  J <- I - (1/n)*(ones%*%t(ones)) 
  Xtilde <- J %*% X #centralized X
  k <- kernel_function(X, ...)
  k <- J %*% k %*% J #centralized K
  
  w_0 <-(1/n)*t(ones)%*%y
  q_tilde <- solve(I + lambda*ginv(k),J%*%y)
  preds <- w_0[1] + q_tilde
  w <- round(ginv(Xtilde)%*%q_tilde,2)
  res_list = list('preds'= preds, 'q_tilde'=q_tilde, 'w_0'= w_0, 'w' = w, 
                  'mean_DM' = mean_DM, 'std_DM' = std_DM)
  return(res_list)
}

## Prediction ------------------------------------------------------------------ 
predict_oos <- function(X_u, X_t, res, kernel.type, kernel_function, ...){
  w_0 <- res$w_0
  q_tilde <- res$q_tilde
  mean_DM <- res$mean_DM
  std_DM <- res$std_DM
  
  K_t <- kernel_function(X_t, ...)
  
  X_u_tilde <- X_u -  outer(rep(1, dim(X_u)[1]), mean_DM)
  X_u_tilde <- X_u_tilde / outer(rep(1, dim(X_u)[1]), std_DM)
  #-------------------------------------------------------------------------------
  if(kernel.type == 'linear'){
    ku <-  X_u_tilde%*%t(X_t)
  } else if(kernel.type == 'inhomogeneous'){
    ku <- (1 + X_u_tilde%*%t(X_t))^d 
  } else if(kernel.type == 'rbf'){
    ku <- kernel_function( X_u_tilde, gamma, X_t) 
  }
  
  preds <- w_0[1] + X_u_tilde %*% ginv(K_t) %*% q_tilde
}

## Cross validation ------------------------------------------------------------
# TODO

# Script -----------------------------------------------------------------------
load(datapath)

df <- Airline

# Preprocessing ----------------------------------------------------------------
y <- df$output
X_names <- fastDummies::dummy_cols(df, select_columns = c('airline')) %>% names
df <- fastDummies::dummy_cols(df %>% scale, select_columns = c('airline'))
names(df) <- X_names
X <- df %>% select(-c(output, airline, airline_1)) %>% as.matrix

#Test Data
X_u <- X[(dim(X)[1]*0.7):dim(X)[1],]
y_u <- y[(dim(X)[1]*0.7):dim(X)[1]]

#Train Data
X_t <- X[1:(dim(X)[1]*0.7),] 
y_t <- y[1:(dim(X)[1]*0.7)]

# Results ----------------------------------------------------------------------
z <- krr(y, X, lambda, kernel_linear) 
a <- krr(y, X, lambda, kernel_rbf, gamma=1/2) 
b <- krr(y, X, lambda, rbfkernel, sigma=1) # Comparison to standard kernel function

res <-  krr(y, X, lambda, kernel_inhomogeneous, d=2)
# Out of Sample Predictions ------------------------------------------------------------------


#Out of sample predictions
w_0 <- res$w_0
q_tilde <- res$q_tilde

predsINH <- predict_oos(w_0[1], X_u, X_t, q_tilde, "inhomogeneous", kernel_inhomogeneous, d=2)
predsRBF <- predict_oos(w_0[1], X_u, X_t, q_tilde, "rbf", kernel_rbf, gamma=0.5)


################################################################################
#-------------------------------------------------------------------------------
#In sample comparison with rdtools and dsmle
#-------------------------------------------------------------------------------
#dsmle
res_pkg <- dsmle::krr(y, X, 10, kernel.type = "nonhompolynom", kernel.degree = 1)
y_hat_dsmle <- res_pkg$yhat

#-------------------------------------------------------------------------------
#rdetools
k <- polykernel(X, 2, Y = NULL)
k <- rbfkernel(X)
r <- rde(k, y, est_y = TRUE, dim_rest = 1, regression=TRUE)
y_hat_rdtools <- r$yh
?rde()





