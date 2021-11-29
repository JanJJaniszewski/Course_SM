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
  col_bests <- tibble(gini_gain = numeric(), split = numeric(), column = character())
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
    column_comp <- tibble(split = column_splits, gini_gain = column_ginis, column=colname)
    column_best <- column_comp %>% filter(gini_gain == max(gini_gain))
    col_bests <- col_bests %>% add_row(column_best)
  }
  split_list <- col_bests %>% filter(gini_gain == max(gini_gain)) %>% .[1,]
  return(split_list)
}

provide_best_split_results_with_predictions <- function(df, splitno=1, max_splitno=4^2){
  if(((splitno * 2 + 1) <= max_splitno) & (nrow(df) > 1)){
    best_split <- define_best_split(df)
    colname <- best_split$column[1]
    best_split <- best_split$split[1]
    split1 <- df[df[,colname] < best_split,]
    split2 <- df[df[,colname] >= best_split,]

    # Printing progress
    print('<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Split node')
    print(paste('Number:', splitno))
    print(paste('Split:', best_split))
    print(paste('Split1 N:', nrow(split1)))
    print(paste('Split2 N:', nrow(split2)))

    split1 <- provide_best_split_results_with_predictions(split1, splitno * 2, max_splitno)
    split2 <- provide_best_split_results_with_predictions(split2, (splitno * 2) + 1, max_splitno)
    out <- rbind(split1, split2)
    return(out)
  }else{
    out <- df
    out$majority_decision <- round(mean(out$Species))
    out$leaf_number <- splitno
    
    print('-------------------------------- Final node')
    print(paste('Number:', splitno))
    print(paste('Mean:', mean(out$Species) %>% round(2)))
    print(paste('N:', nrow(out)))
    
    return (out)
  }
}

output <- provide_best_split_results_with_predictions(df)
