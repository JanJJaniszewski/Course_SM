# Loading files
source('./Week_2/config.R')
source(path_functions)

#Data
y<- Airq[1:30,c(1)]
X <- cbind(1, scale(Airq[1:30, c(2:6)]))
Better_Sub_Selec <- function(y, X, mMax){
  #set beta zero
  B0 <- matrix(1, dim(X)[2], 1)
  Bk <- matrix(1, dim(X)[2], 1)
  #eigenvalue max
  evMAX <- eigen(t(X) %*% X)[[1]][1]
  k <- 1
  mMax <- mMax
  while (k==1|(RSS(X,B0,y)-RSS(X,Bk,y))/RSS(X,B0,y)>1e-6){
    k <- k + 1
    B0 <- Bk
    xTx <- t(X) %*% X
    xTy <- t(X) %*% y
    u <- B0 - 1/evMAX * (xTx) %*% B0 + 1/evMAX * (xTy)
    if (k>2) a <- ((t(B0)%*%xTy)/(t(B0)%*%xTx%*%B0))
    2
    else a <- 1
    Border <- order(abs(u), decreasing = TRUE)
    for (i in 1:dim(B0)[1]){
      if (i<=mMax){
        Bk[Border[i]] <- u[Border[i]]*a[1]
      }
      else{
        Bk[Border[i]] <- 0.0
      }
    }
  }
  y_mean <- mean(y)
  SST <- t(y-y_mean)%*%(y-y_mean)
  R2 <- 1 - (RSS(X, Bk, y) / SST)
  res_list <- list(a, Bk, R2,k)
  return(res_list)
}
Bk <- Better_Sub_Selec(y, X, 6)
Bk
#print(Bk[3])
