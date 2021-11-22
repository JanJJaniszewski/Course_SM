source('./functions.R')
source('./config.R')

# Script -----------------------------------------------------------------------
load(datapath)

df <- Airline

df <- fastDummies::dummy_cols(df, select_columns = c('airline')) %>% dplyr::select(-airline)

df_train <- df %>% filter(year < max(year) - config_years_to_predict_in_testset)
df_test <- df %>% filter(year == max(year) - config_years_to_predict_in_testset)

# Preprocessing ----------------------------------------------------------------
# why is this scaling happening like this?
df_test <- df_test %>% mutate(
  year = (year - mean(df_train$year)) / sd(df_train$year),
  cost = (cost - mean(df_train$cost)) / sd(df_train$cost),
  pf = (pf - mean(df_train$pf)) / sd(df_train$pf),
  lf = (lf - mean(df_train$lf)) / sd(df_train$lf)
)

df_train[,c('year', 'cost', 'pf', 'lf')] <- df_train[,c('year', 'cost', 'pf', 'lf')] %>% scale

# Data analysis ----------------------------------------------------------------
t1 <- Airline %>% group_by(year) %>% summarise_all(c(mean, sd))
plot(t1$year, t1$output_fn1)

# Cross Validation -------------------------------------------------------------

run_cv <- function(df_train, df_test, mygamma, mylambda, myr){
  X_train <- dplyr::select(df_train, -output) %>% as.matrix()
  y_train <- df_train$output
  X_test <- dplyr::select(df_test, -output) %>% as.matrix()
  y_test <- df_test$output
  
  # why are there so many variables here that are not defined in the function?
  res <- krr(y_train, X_train, lambda=mylambda, config_kernel, gamma=mygamma, d=config_d)
  predictions <- predict_oos(X_test, X_train, res, config_kernel, gamma=mygamma, d=config_d)
  mse <- mean((predictions - y_test)^2)
  return(mse)
}

cv <- df_train %>% split(.$year)
cv <- tibble(year = c(1:length(names(cv))), test = cv)

filter_by_year <- function(filteryear, df){df %>% filter(df$year != filteryear)}
cv$train <- map(cv$test, ~ filter_by_year(mean(.$year), df_train))

# Gridsearch -------------------------------------------------------------------
crossv_cv <- tibble()
for (config_lambda in config_lambdas){
  for (config_gamma in config_gammas){
    cv_lambda <- cv
    cv_lambda$lambda <- config_lambda
    cv_lambda$gamma <- config_gamma
    cv_lambda$mse <- map2_dbl(cv$train, cv$test, ~ run_cv(.x, .y, config_gamma, config_lambda))
    crossv_cv <- bind_rows(crossv_cv, cv_lambda)
  }
}

crossv_output <- crossv_cv %>% select(c(lambda, gamma, mse)) %>% group_by(lambda, gamma) %>% summarise_all(mean)
best_lambda <- crossv_output %>% filter(mse == min(crossv_output$mse)) %>% pull(lambda)
best_gamma <- crossv_output %>% filter(mse == min(crossv_output$mse)) %>% pull(gamma)

print(best_lambda)
print(best_gamma)
print(min(crossv_output$mse))

# Final Testset ----------------------------------------------------------------
X_train <- dplyr::select(df_train, -output) %>% as.matrix()
y_train <- df_train$output
X_test <- dplyr::select(df_test, -output) %>% as.matrix()
y_test <- df_test$output

res <- krr(y_train, X_train, lambda = best_lambda, config_kernel, gamma = best_gamma, d=config_d)
predictions <- predict_oos(X_test, X_train, res, config_kernel, gamma = best_gamma, d=config_d)
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
