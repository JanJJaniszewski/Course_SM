# RSS function
RSS <- function(X,B,y) {
  y_hat <- X %*% B
  rss <- t(y - y_hat) %*% (y - y_hat)
  return(rss)
}


MJF <- function(B, A, D, xtx = xtx, xty = xty, yty = yty, lambda = lambda, alpha = alpha,n=n){
  c <- 1/(2*n)*yty + (1/2)*lambda*alpha*sum(abs(B))
  B_u <- 1/2*t(B) %*% A %*% B - 1/n *t(B)%*%xty + c 
  
  return(B_u)
}

MSE <- function(X, Bk, y){
  mse <- RSS(X, Bk, y)/dim(X)[1]
  return(mse)
}

rsquared <- function(X, Bk, y, adjusted = FALSE){
  y_mean <- mean(y)
  SST <- t(y-y_mean)%*%(y-y_mean)
  n <- dim(X)[1]
  p <- dim(X)[2]
  R2 <- 1 - (RSS(X, Bk, y) / SST)
  if(adjusted == FALSE){
    R2 <- R2
  } else {
    R2 <- 1 - ((1 - R2)*(n-1))/(n - p -1)
  }
  return(R2)
}

adj.designmatrix <- function(X, scale = TRUE, intercept = TRUE){
  if(scale == TRUE) {X <- as.matrix(scale(X))}
  if(intercept == TRUE) {X <- cbind(matrix(1, dim(X)[1], 1), X)}
}

elastic_net <-function(y, X, lambda, alpha){
  
  X <- adj.designmatrix(X, scale = TRUE, intercept = TRUE)
  
  #pre-calculations
  p <- dim(X)[2]
  n <- dim(X)[1]
  xtx <- t(X)%*%X
  xty <- t(X)%*%y
  yty <- t(y)%*%y
  
  #set initial coefficient vectors
  B0 <- matrix(1, p, 1)
  Bk <- matrix(1, p, 1)
  epsilon <- 1e-6
  
  k <- 1
  L_b0 <- 10
  L_bk <- 9
  
  while (k == 1 | (L_b0 - L_bk)/ L_b0 > 1e-6){
    k <- k + 1
    B0 <- Bk
    inter <- apply(Bk, 1, function(x) max(abs(x),epsilon))
    I <- diag(p)
    I[1,1] <- 0 
    D <- I/inter
    A <- (1/n)*xtx + (lambda * (1-alpha)) * I + (lambda * alpha) * D
    
    L_b0 = MJF(B0, A, D, xtx = xtx, xty = xty, yty = yty, lambda = lambda, alpha = alpha, n=n)
    Bk <- 1/n * solve(A) %*% xty
    
    L_bk = MJF(Bk, A, D, xtx = xtx, xty = xty, yty = yty, lambda = lambda, alpha = alpha, n=n)
    
    #print(L_bk)
    #print(L_b0)
    
  }
  R2 <- rsquared(X, Bk, y, adjusted = FALSE)
  res_list <- list('coefficients'=Bk,
                   'R2'=R2,
                   'X'=X, 
                   'iterations'= k,
                   'loss' = L_bk)
  return(res_list)
}