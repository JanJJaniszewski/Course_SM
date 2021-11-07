source("functions.R")
load("../supermarket1996.RData")
set.seed(0)
library("tidyverse")

df <- supermarket1996

#Data
y <- df[,c("GROCERY_sum")]
X <- df[,-which(names(df) %in% c("GROCERY_sum","STORE", "CITY", "ZIP", "GROCCOUP_sum", "SHPINDX"))]

start_lambda <- 0
end_lambda <- 3
start_alpha <- 0
end_alpha <- 1

delta_lambda = 1
delta_alpha = 0.2

gridkcv <- function(start_lambda, end_lambda, delta_lambda, start_alpha, end_alpha, delta_alpha){
  grid <- expand.grid(lambda = seq(from = start_lambda, by = delta_lambda, to = end_lambda),
                      alpha = seq(from = start_alpha, by = delta_alpha, to = end_alpha))
  return(grid)
  }

gridkcv(start_lambda, end_lambda, delta_lambda, start_alpha, end_alpha, delta_alpha)

kfold <- function(k, y, X, lambdas, alpha, paramgrid){
  indices <- sample(dim(X)[1])
  folds <- split(indices, ceiling(seq_along(indices)/(length(indices)/k)))
  
  fit_stat <- matrix(1,k,1)
  for(i in 1:k){
    y_train <- y[-folds[[i]]]
    X_train <- X[-folds[[i]],]
    
    y_test <- y[folds[[i]]]
    X_test <- X[folds[[i]],]
    
    fit <- elastic_net(y_train, X_train, 2, 0.4)
    
    coef <- fit[[1]]
    #R2 <- rsquared(adj.designmatrix(X_test, scale = TRUE, intercept = TRUE), fit[[1]], y_test, adjusted = FALSE)
    fit_stat[i] <- MSE(adj.designmatrix(X_test, scale = TRUE, intercept = TRUE), coef, y_test)
  }
  
}
a <- matrix(1,5,1)
a[1] = 5
a
