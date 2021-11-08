# Loading files
source('./Week_2/config.R')
source(path_functions)
source(path_packages)

# Load data
df_input <- loadRData(path_data)

X <- df_input %>% select(-c(STORE, CITY, ZIP, GROCCOUP_sum, SHPINDX, GROCERY_sum)) %>% scale
X <- cbind(X, 1)
y <- df_input %>% pull(GROCERY_sum)

# Initialization
p <- ncol(X)
n <- nrow(X)
beta <- matrix(1, p, 1)
lambda <- 10
alpha <- 0.4
identity <- diag(p)

epsilon = 1e-6
improvement <- epsilon + 1


# computation of stable terms
yty <- t(y) %*% y # y squared
Xty <- t(X) %*% y
XtX <- t(X) %*% X
  
n_1_XtX <- n^(-1) * XtX
lambda_1_alpha_I <- lambda * (1-alpha) * identity
A_p1 <- n_1_XtX + lambda_1_alpha_I
D <- diag(p)

L_new = improvement
while (improvement > epsilon){
  diag(D) <- 1 / (beta %>% abs %>% sapply(max, ...=epsilon))
  A <- A_p1 + (lambda * alpha * D)
  beta <- n^(-1) * solve(A) %*% Xty
  L_old <- L_new
  L_new <- compute_L(beta, A, n, Xty, yty)
  improvement <- (L_new - L_old) / L_old
}

beta %>% mean
