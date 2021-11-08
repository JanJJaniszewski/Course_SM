source("functions.R")
addTaskCallback(function(...) {set.seed(123);TRUE})
load("../supermarket1996.RData")
set.seed(10)
library("tidyverse")
library(ggplot2)
#install.packages("hrbrthemes")
library(hrbrthemes)

df <- supermarket1996

#Data
y <- df[,c("GROCERY_sum")]
X <- df[,-which(names(df) %in% c("GROCERY_sum","STORE", "CITY", "ZIP", "GROCCOUP_sum", "SHPINDX"))]

start_lambda <- 0.1
end_lambda <- 5
start_alpha <- 0
end_alpha <- 0.99

delta_lambda = 1
delta_alpha = 0.1

param_grid <- gridkcv(start_lambda, end_lambda, delta_lambda, start_alpha, end_alpha, delta_alpha)


a <- kfold(5, y, X, param_grid, verbose = FALSE)
a

ggplot(a, aes(lambda, alpha, fill= RMSE)) + 
  geom_tile() + 
  scale_fill_gradient(low = "white", high = "red") +
  ggtitle("Plot 1 : Heatmap of the Alphas vs Lambdas (Filled with the RMSE)")
  
