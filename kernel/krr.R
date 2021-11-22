source('./config.R')
source('./functions.R')

# Script -----------------------------------------------------------------------
load(datapath)

df <- Airline
df <- fastDummies::dummy_cols(df, select_columns = c('airline')) %>% dplyr::select(-airline)

years_to_predict <- max(df$year) - config_years_to_predict_in_testset + 1
df_train <- df %>% filter(year < years_to_predict)
df_test <- df %>% filter(year >= years_to_predict)

# Preprocessing ----------------------------------------------------------------
df_test <- df_test %>% mutate(
  # Scaling test data based on training data averages
  year = (year - mean(df_train$year)) / sd(df_train$year),
  cost = (cost - mean(df_train$cost)) / sd(df_train$cost),
  pf = (pf - mean(df_train$pf)) / sd(df_train$pf),
  lf = (lf - mean(df_train$lf)) / sd(df_train$lf)
)

df_train[,c('year', 'cost', 'pf', 'lf')] <- df_train[,c('year', 'cost', 'pf', 'lf')] %>% scale


# Data analysis ----------------------------------------------------------------
t1 <- Airline %>% group_by(year) %>% summarise_all(c(mean, sd))
plot(t1$year, t1$output_fn1)

Hmisc::describe(df_train)

# Cross Validation -------------------------------------------------------------
cv <- df_train %>% split(.$year)
cv <- tibble(year = c(1:length(names(cv))), test = cv)

filter_by_year <- function(filteryear, df){df %>% filter(df$year != filteryear)}
cv$train <- map(cv$test, ~ filter_by_year(mean(.$year), df_train))

# Gridsearch -------------------------------------------------------------------
crossv_output <- compare_kernels(cv, config_lambdas, my_verbose=T)

crossv_output <- crossv_output %>% select(c(lambda, hyperparameter, mse, kernel)) %>% group_by(kernel, hyperparameter, lambda) %>% summarise_all(mean) 
best_run <- crossv_output %>% filter(mse == min(crossv_output$mse))

print('Best run')
print(best_run)
best_kernel_function <- kernel_rbf # To adjust manually, R is not able to save functions in dataframes

# Final Testset ----------------------------------------------------------------
X_train <- dplyr::select(df_train, -output) %>% as.matrix()
y_train <- df_train$output
X_test <- dplyr::select(df_test, -output) %>% as.matrix()
y_test <- df_test$output

res <- krr(y_train, X_train, hyperpar = best_run$hyperparameter, kernel_function = best_kernel_function, lambda=best_run$lambda)
predictions <- predict_oos(X_test, X_train, res, best_kernel_function, hyperpar = best_run$hyperparameter)
mse <- mean((predictions - y_test)^2)
print(mse)

# Comparison -------------------------------------------------------------------
## LM --------------------------------------------------------------------------
lm_model <- lm(y_train ~ ., data=as.data.frame(cbind(X_train, y_train)))
lm_preds <- predict(lm_model, data=as.data.frame(cbind(X_test, y_test)))
mse <- mean((lm_preds - y_test)^2)
print(mse) # Higher than our model

## DSMLE -----------------------------------------------------------------------
res_pkg <-
  dsmle::krr(y_train, X_train, config_k, kernel.type = "nonhompolynom", kernel.degree = 1)
y_hat_dsmle <- res_pkg$yhat
mse <- mean((y_hat_dsmle - y_train)^2)
print(mse) # Higher than our model
