library(pacman)
p_load(dplyr)
p_load(datasets)
data(iris)
set.seed(1)
df <- iris
df$Species <- df$Species == 'versicolor'
get_gini <- function(df, y_name){
  filter_positives <- df[,y_name] == 1
  filter_negatives <- df[,y_name] == 0
  n <- nrow(df)
  n_positives <- df[filter_positives,] %>% nrow
  n_negatives <- df[filter_negatives,] %>% nrow
  g <- 1 - (n_positives/n)^2 - (n_negatives/n)^2
  return(g)
}
define_best_split <- function(df){
  X_names <- df %>% select(-Species) %>% names
  y <- df %>% pull(Species)
  n <- nrow(df)
  col_bests <- tibble(gini = numeric(), split = numeric(), column = character())
  for(colname in X_names){
    #find mid-point vector
    uniques <- df[,colname] %>% sort %>% unique
    column_splits <- c()
    column_ginis <- c()
    for(i in 1:length(uniques)-1){
      column_splits[i] = mean(uniques[i:i+1])
    }
    for(i in 1:length(column_splits)){
      split1 <- df[df[,colname] < column_splits[i],]
      split2 <- df[df[,colname] >= column_splits[i],]
      # count Gini information
      g_root <- get_gini(df, 'Species')
      g_split1 <- get_gini(split1, 'Species')
      g_split2 <- get_gini(split2, 'Species')
      g_gain <- g_root - g_split1 * (nrow(split1)/n) - g_split2 * (nrow(split2)/n)
      column_ginis[i] <- g_gain
    }
    column_comp <- tibble(split = column_splits, gini = column_ginis, column=colname)
    column_best <- column_comp %>% filter(gini == min(gini))
    col_bests <- col_bests %>% add_row(column_best)
  }
  best_column <- col_bests %>% filter(gini == min(col_bests$gini)) %>% select(column) 
  best_split <- col_bests %>% filter(gini == min(col_bests$gini)) %>% select(split) 
  split_list <- list('b_column' = best_column, 'b_split' = best_split)
  return(split_list)
}

grow_tree <- function(max_depth, min_samples_split, df){
  df0 <- df
  k <- 0
  while(k < max_depth){
    k <- k + 1
    for(i in 2^(k-1):2^(k)-1){
      res <- define_best_split(paste0("df",i))
      b_column <- res$bcolumn
      b_split <- res$bsplit
      paste0("df",i*2)[which(paste0("df",i)$bcolumn > b_split)]
      paste0("df",i*2 + 1)[which(paste0("df",i)$bcolumn > b_split)]
    }
    best
  }
  
}
define_best_split(df)

