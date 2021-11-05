#### RSS ####
RSS <- function(X,B,y) {
  rss <- t(y - X %*% B) %*% (y - X %*% B)
  return(rss)
}

