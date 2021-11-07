setwd("~/AA BDS/1-2 Supervised Machine Learning/Week 2")
getwd()
#install.packages("tidyverse")
library("tidyverse")
set.seed(0)

#data prep
load("supermarket1996.RData")
df <- supermarket1996
y <- df[c(1:nrow(df)),c(1)]
x <- select(supermarket1996, -c(1,2,3,4,5,30))

X <- x[c(1:nrow(x)),c(1:ncol(x))]
X <- scale(X)
XT <- t(X)
XTX <- XT%*%X
print(XTX)
p= dim(X)[2]
B = matrix(1,p,1)
n = 77

#Functions

#LSS
LSS <- function(y,X,B, lambda = 2){
  ldist = 0
  for (i in B){ldist = ldist + abs(i)}
  return(t(y - X %*% B)%*%(y - X %*% B)+lambda*ldist)
}

MMel <- function(y,X,n, epsilon = 1e-6, lambda = 2,alpha = 0.4){
  #set beta zero
  B0 <- matrix(1,p, 1)
  Bk <- matrix(1,p,1)
  I <- diag(p)
  k = 1
  
  while (k==1|((LSS(y,X,Bk)-LSS(y,X,B0))%/%LSS(y,X,B0))>epsilon) {
    k <- k + 1
    inter = apply(B0,1,function(x) max(abs(x),epsilon))
    D=I/inter
    print(diag(D))
    A =(1/n)*XTX+lambda*(1-alpha)*I+lambda*alpha*D
    Bk <- n^(-1)*solve(A)%*%(XT%*%y)
    B0 <- Bk
    dim(Bk)
  }
  return(Bk)
}

MMel(y,X,77)
