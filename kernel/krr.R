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
  I <- diag(n)
  
  J <- I - (1/n)*(ones%*%t(ones)) 
  # Xtilde <- J %*% X #centralized X
  kkt <- kernel_function(X, ...)
  
  w_0 <-(1/n)*t(ones)%*%y
  q_tilde <- solve(I + lambda*ginv(kkt),J%*%y)
  preds <- w_0[1] + q_tilde
  w <- round(ginv(X)%*%q_tilde,2)
  res_list = list('preds'= preds, 'q_tilde'=q_tilde, 'w_0'= w_0, 'w' = w)
  return(res_list)
}

## Prediction ------------------------------------------------------------------ 
predict_oos <- function(w_0, X_u, X_t, q_tilde, kernel.type, kernel_function, ...){
  K_t <- kernel_function(X_t, ...)
  #-------------------------------------------------------------------------------
  if(kernel.type == 'linear'){
    ku <- X_u%*%t(X_t)
  } else if(kernel.type == 'inhomogeneous'){
    ku <- (1 + X_u%*%t(X_t))^d 
  } else if(kernel.type == 'rbf'){
    ku <- kernel_function(X_u, gamma, X_t) 
  }
  
  preds <- w_0[1] + ku %*% ginv(K_t) %*% q_tilde
}

## Cross validation ------------------------------------------------------------
# TODO

# Script -----------------------------------------------------------------------
load(datapath)

df <- Airline

# Preprocessing ----------------------------------------------------------------
df[, paste0('airline', 1:6)] <- NA
df[, (dim(df)[2]-6+1):dim(df)[2]] <- sapply(names(df)[(dim(df)[2]-6+1):dim(df)[2]], function(x) as.integer(substr(x,8,8) == df$airline))

# Initialisation
X <- df %>% dplyr::select(-c("airline", "airline1", "output"))

X[,1:4] <- scale(X[,1:4])
X <- as.matrix(X)
y <- df$output

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
predsINH <- predict_oos(w_0[1], X_u, X_t, q_tilde, "inhomogeneous", kernel_inhomogeneous, d=2)
predsRBF <- predict_oos(w_0[1], X_u, X_t, q_tilde, "rbf", kernel_rbf, gamma=0.5)

dsmle::krr


################################################################################
#-------------------------------------------------------------------------------
#In sample comparison with rdtools and dsmle
#-------------------------------------------------------------------------------
#dsmle
res_pkg <- dsmle::krr(y, X, 10, kernel.type = "inhomopolynom", kernel.degree = 2)
y_hat_dsmle <- res_pkg$yhat

#-------------------------------------------------------------------------------
#rdetools
k <- polykernel(X, 2, Y = NULL)
k <- rbfkernel(X)
r <- rde(k, y, est_y = TRUE, dim_rest = 1, regression=TRUE)
y_hat_rdtools <- r$yh

