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

mse_linear_krr <- mean((f_hat_target_Linear_KRR - y)^2)
mse_gaussian_krr <- mean((f_hat_target_Gaussian_KRR - y)^2)
mse_laplacian_krr <- mean((f_hat_target_Laplacian_KRR - y)^2)

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


# Request of Caio
krr <- function(y, X, lambda, kernel_function, output_one_value = TRUE, ...){
  outputs <- list()
  n <- nrow(X)
  ones <- matrix(1, n, 1)
  I <- diag(n)
  
  J <- I - (1/n)*(ones%*%t(ones)) 
  Xtilde <- J %*% X #centralized X
  kkt <- kernel_function(X, ...)
  
  w_0 <-(1/n)*t(ones)%*%y
  q_tilde <- solve(I + lambda*ginv(kkt))%*%J%*%y
  w <- round(ginv(Xtilde)%*%q_tilde,2)
  
  outputs['w0'] <- w_0
  outputs['w'] <- w
  outputs['q_tilde'] <- q_tilde
  
  if(output_one_value){return(w)}else{return(outputs)}
}

outputs <- krr(y, X, lambda, kernel_rbf, gamma=1/2, output_one_value=FALSE)
w0 <- outputs$w0
q_tilde <- outputs$q_tilde
w <- outputs$w
X_train <- X
X_test <- X

# Function (throws error)
compute_q_test <- function(X_train, X_test, q_tilde, w0){
  n <- nrow(X)
  ones <- matrix(1, n, 1)
  k_u <- kernel_rbf(X_test, gamma=1/2)
  k <- kernel_rbf(X_train, gamma=1/2)
  q_u <- w0 %*% ones + k_u %*% k^(-1) %*% q_tilde # something wrong here as non conformable arguments
  q_u <- w0 + k_u %*% k^(-1) %*% q_tilde # also does not work
  return(q_u)
}