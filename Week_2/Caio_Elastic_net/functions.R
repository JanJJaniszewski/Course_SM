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

gridkcv <- function(start_lambda, end_lambda, delta_lambda, start_alpha, end_alpha, delta_alpha){
  grid <- expand.grid(lambda = seq(from = start_lambda, by = delta_lambda, to = end_lambda),
                      alpha = seq(from = start_alpha, by = delta_alpha, to = end_alpha))
  return(grid)
}

kfold <- function(k, y, X, lambdas, alpha){
  indices <- sample(dim(X)[1])
  folds <- split(indices, ceiling(seq_along(indices)/(length(indices)/k)))
  
  grid_size <- seq(dim(param_grid)[1])
  
  fit_stat <- param_grid
  
  for(j in grid_size){
    temp <- matrix(1, grid_size, 1)
    temp2 <- matrix(1, grid_size, 1)
    for(i in 1:k){
      
      y_train <- y[-folds[[i]]]
      X_train <- X[-folds[[i]],]
      
      y_test <- y[folds[[i]]]
      X_test <- X[folds[[i]],]
      
      fit <- elastic_net(y_train, X_train, param_grid[[1]][j], param_grid[[2]][i])
      
      coef <- fit[[1]]
      R2 <- rsquared(adj.designmatrix(X_test, scale = TRUE, intercept = TRUE), fit[[1]], y_test, adjusted = FALSE)
      mse <- MSE(adj.designmatrix(X_test, scale = TRUE, intercept = TRUE), coef, y_test)
      #print(param_grid[[2]][i]) 
      temp[[i]] <- mse
      temp2[[i]] <- R2
    }
    
    fit_stat[j,'MSE'] <- mean(temp)
    fit_stat[j,'R2'] <- mean(temp2)
  }
  return(fit_stat)
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

