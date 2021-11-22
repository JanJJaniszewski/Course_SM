p_load('KRMM')

#Linear kernel
Linear_KRR_model_train = Kernel_Ridge_MM(Y_train=y,
                                         Matrix_covariates_train=X, method="RR-BLUP")
f_hat_target_Linear_KRR = Predict_kernel_Ridge_MM( Linear_KRR_model_train,
                                                   Matrix_covariates_target=X )

#Gaussian kernel
Gaussian_KRR_model_train = Kernel_Ridge_MM( Y_train=y,kernel = 'Gaussian',
                                            Matrix_covariates_train=X, method="RKHS", rate_decay_kernel=5.0)
f_hat_target_Gaussian_KRR = Predict_kernel_Ridge_MM( Gaussian_KRR_model_train, 
                                                     Matrix_covariates_target=X )

# Laplacian kernel
Laplacian_KRR_model_train = Kernel_Ridge_MM( Y_train=y,kernel = 'Laplacian',
                                            Matrix_covariates_train=X, method="RKHS", rate_decay_kernel=5.0)
f_hat_target_Laplacian_KRR = Predict_kernel_Ridge_MM( Gaussian_KRR_model_train, 
                                                     Matrix_covariates_target=X )

# Compare prediction based on linear (i.e. RR-BLUP) and Gaussian kernel
par(mfrow=c(1,3))
plot(f_hat_target_Linear_KRR, y)
plot(f_hat_target_Gaussian_KRR, y)
plot(f_hat_target_Laplacian_KRR, y)

mse_m1 <- mean((f_hat_target_Linear_KRR - y)^2)
mse_m2 <- mean((f_hat_target_Gaussian_KRR - y)^2)
mse_m3 <- mean((f_hat_target_Laplacian_KRR - y)^2)

comparison <- data.frame(mse_linear_krr, mse_gaussian_krr, mse_laplacian_krr)

comparison

# We now know that the gaussian kernel fits best. We will now start experimenting with different sigma's
#Gaussian kernel
gridsearch_analysis <- list()
par(mfrow=c(2,2))
for (sigma in 5 * 10 ^ c(-2, -1, 1, 2)){
  Gaussian_KRR_model_train = Kernel_Ridge_MM( Y_train=y,kernel = 'Gaussian', init_sigma2K = sigma,
                                              Matrix_covariates_train=X, method="RKHS", rate_decay_kernel=5.0)
  prediction = Predict_kernel_Ridge_MM( Gaussian_KRR_model_train, 
                                                       Matrix_covariates_target=X )
  plot(prediction, y)
  gridsearch_analysis[toString(sigma)] <- mean((prediction - y) ^ 2)
}

gridsearch_analysis

# source('./Week_2/Caio_Elastic_net/functions.R')
# 
# library(tidyverse)
# library(optimbase)
# library(matlib)
# library(fastDummies)
# library(rdetools)
# 
# # Kernel function
# kernel_rbf <- function (XtX, gamma) {
#   XX <- matrix(1, n) %*% diag(XtX)
#   k <- exp(-(XX - 2 * XtX + t(XX)) * gamma)
#   return(k)
# }
# 
# # Load and define data
# load('./Week_3/Airline.RData')
# input <- Airline
# y <- input$output
# X_names <- fastDummies::dummy_cols(input, select_columns = c('airline')) %>% names
# input <- fastDummies::dummy_cols(input %>% scale, select_columns = c('airline'))
# names(input) <- X_names
# X <- input %>% select(-c(output, airline, year, airline_1)) %>% as.matrix 
# n <- nrow(X)
# lambda <- 0.5
# gamma <- 1
# ones <- matrix(1,n,1)
# w0 <- mean(y)
# 
# # Transformation
# XXt <- tcrossprod(X)
# t(X) %*% X
# X %*% t(X)
# 
# all(XtX == X %*% t(X))
# 
# X_rbf <- kernel_rbf(XtX, gamma)
# X_rbf_X_rbf_t <- X_rbf%*%t(X_rbf)
# q_tilde <- (I + lambda * (X_rbf_X_rbf_t)^(-1))^(-1) %*% J %*% y
# y
# 
# I <- diag(n)
# J <- I - n^(-1) * matrix(1, n, n)
# ones_vector <- matrix(1, 1, n)
# X_tilde <- J %*% X
# 
# 
# 
# L_ridge <- function(w0, q_tilde){
#   
# }
