source('./Week_2/Caio_Elastic_net/functions.R')

library(tidyverse)
library(optimbase)
library(matlib)
library(fastDummies)
library(rdetools)

# Kernel function
kernel_rbf <- function (XtX, gamma) {
  XX <- matrix(1, n) %*% diag(XtX)
  k <- exp(-(XX - 2 * XtX + t(XX)) * gamma)
  return(k)
}

# Load and define data
load('./Week_3/Airline.RData')
input <- Airline
y <- input$output
X_names <- fastDummies::dummy_cols(input, select_columns = c('airline')) %>% names
input <- fastDummies::dummy_cols(input %>% scale, select_columns = c('airline'))
names(input) <- X_names
X <- input %>% select(-c(output, airline, year, airline_1)) %>% as.matrix 
n <- nrow(X)
lambda <- 0.5
gamma <- 1
ones <- matrix(1,n,1)
w0 <- mean(y)

# Transformation
XXt <- tcrossprod(X)
t(X) %*% X
X %*% t(X)

all(XtX == X %*% t(X))

X_rbf <- kernel_rbf(XtX, gamma)
X_rbf_X_rbf_t <- X_rbf%*%t(X_rbf)
q_tilde <- (I + lambda * (X_rbf_X_rbf_t)^(-1))^(-1) %*% J %*% y
y

I <- diag(n)
J <- I - n^(-1) * matrix(1, n, n)
ones_vector <- matrix(1, 1, n)
X_tilde <- J %*% X



L_ridge <- function(w0, q_tilde){
  
}