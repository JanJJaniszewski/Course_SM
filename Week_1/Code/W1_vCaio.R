setwd('../Week 1')
getwd()
load('Airq_numeric.RData')
set.seed(0)

#data prep
rows <- nrow(Airq)
cols <- ncol(Airq)
X <- matrix(Airq[c(1:rows),c(2:cols)], nrow = rows)
X <- scale(X)
y <- matrix(Airq[c(1:rows),c(1)])

#Functions
##Residuals Sums of Squares
RSS <- function(y,X,B){
  return(t(y - X %*% B)%*%(y - X %*% B))
}

##MM Algo
#Functions
##Calculates residuals sums of squares from estimated coefficients
##MM Algo
MM <- function(y,X,epsilon = NULL,B_0 = NULL){
  X <- cbind(1,X)
  if(is.null(epsilon)) epsilon <- 0.01
  if(is.null(B_0)) B_0 <- matrix(1, nrow=ncol(X))
  lambda <- eigen(t(X)%*%X)[[1]][1]
  value <- 0
  k <- 1
  while(k == 1 | value > epsilon){
    k <- k + 1
    pRss <- RSS(y, X, B_0)
    u <- B_0 - ((1/lambda)*(t(X)%*%X)%*%B_0) + (1/lambda)*(t(X)%*%y)
    cRss <- RSS(y, X, B_t)
    value <- (pRss-cRss)/pRss
    B_0 <- B_t
    print(cRss)}
  y_mean <- mean(y)
  SST <- t(y-y_mean)%*%(y-y_mean) 
  r2 <- 1 - cRss/SST
  return(B_t)
  return(r2)
  return(k)
  sprintf("Coefficients: %f",B_t)
  sprintf("R2: %f",r2)
}

#Initialise
epsilon <- 0.0001
value <- 100

#Results
MM(y,X)

