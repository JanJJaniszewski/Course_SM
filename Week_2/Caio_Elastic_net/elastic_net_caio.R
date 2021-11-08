setwd("~/AA BDS/1-2 Supervised Machine Learning/Week 2")
load("supermarket1996.RData")
set.seed(0)
library("tidyverse")

df <- supermarket1996

# RSS function
RSS <- function(X,B,y) {
  rss <- t(y - X %*% B) %*% (y - X %*% B)
  return(rss)
}
#Data
y <- df[,c("GROCERY_sum")]
X <- df[,-which(names(df) %in% c("GROCERY_sum","STORE", "CITY", "ZIP", "GROCCOUP_sum", "SHPINDX"))]

#Functions
MJF <- function(B, A, D, xtx = xtx, xty = xty, yty = yty, lambda = lambda, alpha = alpha,n=n){
  c <- 1/(2*n)*yty + (1/2)*lambda*alpha*sum(abs(B))
  B_u <- 1/2*t(B) %*% A %*% B - 1/n *t(B)%*%xty + c 
  
  return(B_u)
}


elastic_net <-function(y, X, lambda, alpha){
  X <- as.matrix(scale(X))
  X <- cbind(matrix(1, dim(X)[1], 1), X)
  
  #pre-calculations
  p <- dim(X)[2]
  n <- dim(X)[1]
  xtx <- t(X)%*%X
  xty <- t(X)%*%y
  yty <- t(y)%*%y
  
  #set beta zero
  B0 <- matrix(1, p, 1)
  Bk <- matrix(1, p, 1)
  epsilon <- 1e-6
  #eigenvalue max
  k <- 1
  L_b0 <- 10
  L_bk <- 9
  
  while (k == 1 | (L_b0 - L_bk)/ L_b0 > 1e-6){
    k <- k + 1
    B0 <- Bk
    inter <- apply(Bk, 1, function(x) max(abs(x),epsilon))
    I <- diag(p)
    D <- I/inter
    A <- (1/n)*xtx + (lambda * (1-alpha)) * I + (lambda * alpha) * D
    
    L_b0 = MJF(B0, A, D, xtx = xtx, xty = xty, yty = yty, lambda = lambda, alpha = alpha, n=n)
    
    Bk <- 1/n * solve(A) %*% xty
    
    L_bk = MJF(Bk, A, D, xtx = xtx, xty = xty, yty = yty, lambda = lambda, alpha = alpha, n=n)
    
    print(L_bk)
    print(L_b0)
    
  }
  
  #y_mean <- mean(y)
  #SST <- t(y-y_mean)%*%(y-y_mean)
  #R2 <- 1 - (RSS(X, Bk, y) / SST)
  #res_list <- list(a, Bk, R2,k)
  return(Bk)
}


Bk <- elastic_net(y, X, 2, 0.4)
Bk
mean(y)
