pacman::p_load(mlbench,
               tidyverse, tidymodels, parsnip, mlbench,
               caret, modelr,
               rbenchmark,
               parallel,
               doParallel)

run_cv <- function(X, y){
  X <- as.matrix(X)
  #Gaussian kernel
  Gaussian_KRR_model_train = Kernel_Ridge_MM( Y_train=y,kernel = 'Gaussian',
                                              Matrix_covariates_train=X, method="RKHS", rate_decay_kernel=5.0)
  prediction = Predict_kernel_Ridge_MM( Gaussian_KRR_model_train, Matrix_covariates_target=X )
  mse_gaussian_krr <- mean((prediction - y)^2)
  return(mse_gaussian_krr)
}

# define the independent variables and dependent variable
set.seed(2021)
train_test_split <- initial_split(df, prop = 0.8)

df_train <- training(train_test_split)

df_test <- testing(train_test_split)

cv  <- crossv_kfold(df_train, k = 5)

cv$y_train <- map(cv$train, ~ .$data$airline)
cv$X_train <- map(cv$train, ~ .$data %>% dplyr::select(-airline))
cv$y_test <- map(cv$test, ~ .$data$airline)
cv$X_test <- map(cv$test, ~ .$data %>% dplyr::select(-airline))

cv$mse <- map2(cv$X_train, cv$y_train, ~ run_cv(.x, .y))