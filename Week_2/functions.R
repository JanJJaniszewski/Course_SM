#### RSS ####
RSS <- function(X,B,y) {
  rss <- t(y - X %*% B) %*% (y - X %*% B)
  return(rss)
}

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

compute_L <- function(beta, A, n, xty, yty){
  before_brackets <- (1/2) * t(beta) # 1/2 bT
  after_brackets <- n^(-1) * t(beta) %*% Xty # 
  
  constant <- (2*n)^(-1) %*% yty + (1/2) * lambda * alpha * sum(abs(beta))
  
  out <- before_brackets %*% A %*% beta - after_brackets + constant
  return(out)
}