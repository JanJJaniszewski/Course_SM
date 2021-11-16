# Configuration ----------------------------------------------------------------
datapath <- "./Airline.RData"

# Packages ---------------------------------------------------------------------
#install.packages(c("SVMMaj", "ISLR", "plotrix"))
library("MASS")
library("dsmle")
library("tidyverse")

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

# Kernels ----------------------------------------------------------------------
## Inhomogeneous kernel --------------------------------------------------------
inho_krnl <- function(X, d){
  A <- 1 + X%*%t(X)
  return(A^d)
}
## linear kernel ---------------------------------------------------------------
xxt <- X%*%t(X)

## rbf kernel ------------------------------------------------------------------
kernel_rbf <- function (XtX, gamma) {
  XX <- matrix(1, n) %*% diag(XtX)
  k <- exp(-(XX - 2 * XtX + t(XX)) / gamma)
  return(k)
}

# Script -----------------------------------------------------------------------
load(datapath)

df <- Airline
#-------------------------------------------------------------------------------
#Preprocessing
df[, paste0('airline', 1:6)] <- NA
df[, (dim(df)[2]-6+1):dim(df)[2]] <- sapply(names(df)[(dim(df)[2]-6+1):dim(df)[2]], function(x) as.integer(substr(x,8,8) == df$airline))

# df[, paste0('year', 1:15)] <- NA
# df[, (dim(df)[2]-15+1):dim(df)[2]] <- sapply(names(df)[(dim(df)[2]-15+1):dim(df)[2]], function(x) as.integer(substr(x,5,6) == df$year))
# df <- subset(df, select = -c(airline1, year1))

#Initialisation
X <- df %>% select(-c("airline", "airline1", "output"))

X[,1:4] <- scale(X[,1:4])
X <- as.matrix(X)

Y <- df$output
#-------------------------------------------------------------------------------

lambda <- 10
kkr <- function(Y, X, lambda, kernel, linear = TRUE){
  n <- dim(X)[1]
  ones <- matrix(1, n, 1)

  I <- diag(n)
  J <- I - (1/n)*(ones%*%t(ones)) 
  Xtilde <- J%*%X #centralized X
  if(linear == TRUE){
   kkt <- Xtilde%*%t(Xtilde) 
  } else {kkt <- kernel}

  w_0 <-(1/n)*t(ones)%*%Y
  q_tilde <- solve(I + lambda*ginv(kkt))%*%J%*%Y
  w <- round(ginv(Xtilde)%*%q_tilde,2)
return(w)
}
#-------------------------------------------------------------------------------
#inhomogenous kernel

K <- inho_krnl(X, 2)

#Results
kkr(Y, X, 10, K, linear = FALSE)




