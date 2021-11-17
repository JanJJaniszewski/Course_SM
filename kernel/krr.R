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

# Functions --------------------------------------------------------------------
## Loss Function ---------------------------------------------------------------
L_ridge <- function(w_0, qtilde, y, lambda, q_tilde){
  y_w0 <- y-ones%*%w_0
  JI_qtilde <- J%*%y - q_tilde 
  xtx <- X%*%t(X)
  loss <- y_w0%*%t(y_w0) + JI_qtilde%*%t(JI_qtilde) + lambda*q_tilde%*%solve(xtx)%*%q_tilde
}

## Efficient Matrix Multiplication ---------------------------------------------
matrix_power <- function(A, d){
  result <- diag(dim(A)[1])
  while(d > 0){
    if(d %% 2 == 1){
      result <- result %*% A
    }
    A = A * A
    d = d / 2
  }
  return(result)
}

## Kernels ---------------------------------------------------------------------
kernel_inhomogeneous <- function(X, d){
  A <- 1 + X%*%t(X)
  return(A^d)
}

kernel_linear <- function(X){
  xxt <- X%*%t(X)
}

kernel_rbf <- function (X, gamma) {
  XXt <- tcrossprod(X)
  n <- nrow(X)
  XX <- matrix(1, n) %*% diag(XXt)
  k <- exp(-(XX - 2 * XXt + t(XX)) * gamma)
  return(k)
}

## Ridge Regression -------------------------------------------------------------
krr <- function(y, X, lambda, kernel_function, ...){
  n <- nrow(X)
  ones <- matrix(1, n, 1)
  I <- diag(n)
  
  J <- I - (1/n)*(ones%*%t(ones)) 
  Xtilde <- J %*% X #centralized X
  kkt <- kernel_function(X, ...)
  
  w_0 <-(1/n)*t(ones)%*%y
  q_tilde <- solve(I + lambda*ginv(kkt))%*%J%*%y
  w <- round(ginv(Xtilde)%*%q_tilde,2)
  return(w)
}

## Prediction ------------------------------------------------------------------ 
# TODO

## Cross validation ------------------------------------------------------------
# TODO

# Script -----------------------------------------------------------------------
load(datapath)

df <- Airline

# Preprocessing ----------------------------------------------------------------
df[, paste0('airline', 1:6)] <- NA
df[, (dim(df)[2]-6+1):dim(df)[2]] <- sapply(names(df)[(dim(df)[2]-6+1):dim(df)[2]], function(x) as.integer(substr(x,8,8) == df$airline))

# df[, paste0('year', 1:15)] <- NA
# df[, (dim(df)[2]-15+1):dim(df)[2]] <- sapply(names(df)[(dim(df)[2]-15+1):dim(df)[2]], function(x) as.integer(substr(x,5,6) == df$year))
# df <- subset(df, select = -c(airline1, year1))

# Initialisation
X <- df %>% dplyr::select(-c("airline", "airline1", "output"))

X[,1:4] <- scale(X[,1:4])
X <- as.matrix(X)
y <- df$output

# Results ----------------------------------------------------------------------
krr(y, X, lambda, kernel_rbf, gamma=1/2)
krr(y, X, lambda, rbfkernel, sigma=1) # Comparison to standard kernel function
krr(y, X, lambda, kernel_inhomogeneous, d=1)

f1 <- function()


