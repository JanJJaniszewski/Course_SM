# Info on GLMNET package (to include in the final report?):
# The glmnet algorithms use cyclical coordinate descent, which successively optimizes the objective function over each parameter with others fixed, and cycles repeatedly until convergence. The package also makes use of the strong rules for efficient restriction of the active set.

# Loading files
source('./Week_2/config.R')
source(path_functions)
source(path_packages)

comparison_alpha <- 1
comparison_lambda_cv <- c(1, 10, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9)
comparison_nfolds <- 5
comparison_should_standardize <- TRUE

# Load data
df_input <- loadRData(path_data)

X <- df_input %>% select(-c(STORE, CITY, ZIP, GROCCOUP_sum, SHPINDX, GROCERY_sum)) %>% scale
X <- cbind(X, 1)
y <- df_input %>% pull(GROCERY_sum)

# Using glmnet to directly perform CV
comparison_cv <- cv.glmnet(x=X, y=y,alpha=comparison_alpha,type.measure = 'mse', 
                nfolds = comparison_nfolds, lambda = comparison_lambda_cv,
                standardize=comparison_should_standardize)

# Best lambda according to the comparison CV
comparison_lambda <- comparison_cv$lambda.min

comparison_model <- glmnet(x=X, y=y,alpha=comparison_alpha, nfolds = 1, lambda = comparison_lambda,type.measure = 'mse',
                              standardize=comparison_should_standardize)

# Output of the comparison model
comparison_output <- tibble(Lambda=comparison_cv$lambda,MSE=comparison_cv$cvm)

# Final Pseudo R2
comparison_model$dev.ratio

# Final comparison Log(lambda) vs. MSE
plot(comparison_cv,xvar="lambda",label=TRUE)

# Coefficients of final model
coef(comparison_model)