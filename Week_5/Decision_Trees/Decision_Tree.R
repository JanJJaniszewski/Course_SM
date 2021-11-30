library(pacman)
library(caret)
p_load(dplyr)
p_load(datasets)
data(iris)
set.seed(1)

load("../data/smoking.RData")


################################################################################
#Functions----------------------------------------------------------------------
################################################################################
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
  X_names <- df %>% select(-target) %>% names
  y <- df %>% pull(target)
  n <- nrow(df)
  col_bests <- tibble(gini_gain = numeric(), split = numeric(), column = character())
  for(colname in X_names){
    #find mid-point vector
    uniques <- df[, colname] %>% sort %>% unique
    column_splits <- c()
    column_ginis <- c()
    for(i in 1:length(uniques)-1){
      column_splits[i] = mean(uniques[i:i+1])
    }
    for(i in 1:length(column_splits)){
      split1 <- df[df[,colname] < column_splits[i],]
      split2 <- df[df[,colname] >= column_splits[i],]
      # count Gini information
      g_root <- get_gini(df, 'target')
      g_split1 <- get_gini(split1, 'target')
      g_split2 <- get_gini(split2, 'target')
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

provide_best_split_results_with_predictions <- function(df, splitno=1, max_splitno=4^2, min_samples=10){
  if(((splitno * 2 + 1) <= max_splitno) & (nrow(df) > min_samples)){
    best_split <- define_best_split(df)
    colname <- best_split$column[1]
    best_split <- best_split$split[1]
    split1 <- df[df[,colname] < best_split,]
    split2 <- df[df[,colname] >= best_split,]
    
    # Printing progress
    print('<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Split node')
    print(paste('Number:', splitno))
    print(paste('Splitrule: Split 1 <', best_split, 'and Split 2 >=', best_split))
    print(paste('Splitcolumn:', colname))
    print(paste('First Split (', splitno * 2, '): N =', nrow(split1)))
    print(paste('Spl.1 Class 1 (', splitno * 2, '): N =', (sum(split1$target))/nrow(split1)))
    print(paste('Second Split (', splitno * 2 + 1, '): N =', nrow(split2)))
    print(paste('Spl.2 Class 1 (', splitno * 2, '): N =', (sum(split2$target))/nrow(split2)))
    
    split1 <- provide_best_split_results_with_predictions(split1, splitno * 2, max_splitno, min_samples)
    split2 <- provide_best_split_results_with_predictions(split2, (splitno * 2) + 1, max_splitno, min_samples)
    out <- rbind(split1, split2)
    return(out)
  }else{
    out <- df
    out$majority_decision <- round(mean(out$target))
    out$leaf_number <- splitno
    
    print('-------------------------------- Final node')
    print(paste('Number:', splitno))
    print(paste('Mean:', mean(out$target) %>% round(2)))
    print(paste('N:', nrow(out)))
    
    return (out)
  }
}

################################################################################
#Script-------------------------------------------------------------------------
################################################################################


################################################################################
#Preprocessing------------------------------------------------------------------
################################################################################
df <- smoking
df <- as.data.frame(df)
class(df)

df$intention_to_smoke <- as.integer(df$intention_to_smoke == 'yes')

df$friends_smoke <- as.integer(df$friends_smoke == 'one or more') 
df$lied_to_parents <- as.integer(df$lied_to_parents == 'yes') 

names(df)[names(df) == "intention_to_smoke"] <- "target"
##

df_2 <- iris
df_2$Species <- df_2$Species == 'versicolor'
names(df_2)[names(df_2) == "Species"] <- "target"
################################################################################

output <- provide_best_split_results_with_predictions(df, max_splitno=4^2, min_samples = 5)
#output2 <- provide_best_split_results_with_predictions(df_2, max_splitno=4^3, min_samples = 5)
confusionMatrix(table(output$majority_decision, output$target))
output %>% group_by(leaf_number) %>% count()


library("rpart")

mytree <- rpart(
  target ~ ., 
  data = df, 
  method = "class",
  minsplit = 5,
  maxdepth = 4^2
)

predictions <- predict(mytree, type = "class")

confusionMatrix(table(predictions,df$target))

