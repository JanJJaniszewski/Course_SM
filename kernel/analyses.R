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