##################################################################################################
library(roxygen2)


#' Create table of differentially expressed tags using nonparametric method
#' 
#' @description 
#' 
#' @param data_matrix Matrix, the data in form a matrix with integer as data type
#' @param case String, the case condition
#' @param groups A vector of string, the group information of the column in data_matrix
#'
#' @usage function(data_matrix, case, groups, method=c('fisher', 'Wilcoxon'))
#'
#' @return a data frame with row name the same as the data_matrix, and p-value,
#' adjusted p-value using BH method
nonpar_test <- function(data_matrix, case, groups, method='fisher') {
  # collapse all other group category make it a two group design
  groups <- ifelse(groups==case, case, 'control group')
  design <- model.matrix(~0+groups)
  
  
  res <- data.frame(PValue = rep(NA, nrow(data_matrix)))
  rownames(res) <- rownames(data_matrix)
  
  
  # calculate fold-change
  # caution: the statistical test is not about fold-change it is for the odds ratio
  # which is whether the gene is more likely to be expressed in case group than
  # other genes
  average_matrix <- 1/ (matrix(1, nrow = nrow(data_matrix), ncol = nrow(design)) %*% design)
  average_counts <- (data_matrix %*% design) * average_matrix
  # whether the  gene-specific group-aggregated counts is 0 replace it with 1 for 
  # more meaningful log fold change value
  average_counts[average_counts == 0] <- 1
  res$logFC <- log(average_counts[ ,1]/average_counts[ ,2])
  
  
  if (method == 'fisher') {
    # collapse the matrix into a 2 by 2 table with 
    #                    case group            control group
    #===================================================================
    #    gene X       |        a           |           b
    #--------------------------------------------------------------------
    #    Other genes  |        c           |           d
    
    # sum the count in the same group, and now data_matrix should be a n * 2 matrix
    data_matrix <- data_matrix %*% design
    
    for (rowname in rownames(data_matrix)) {
      other_rows <- colSums(data_matrix[-which(rownames(data_matrix) == rowname), ], na.rm = T)
      test_matrix <- rbind(data_matrix[which(rownames(data_matrix) == rowname), ], other_rows)
      res[rowname, 'PValue'] <- fisher.test(test_matrix)$p.value
    }
  } else if (method == 'Wilcoxon') {
    # each row of x is the 
    x <- data_matrix[, design[ ,1] == 1]
    y <- data_matrix[, design[ ,2] == 1]
    for (rowname in rownames(data_matrix)) {
      res[rowname, 'PValue'] <- wilcox.test(x[rowname, ], y[rowname, ])$p.value
    }
  } else {
    stop('method ')
  }
  
  # If the expected number of observations in any category is too small,
  # the G-test may give inaccurate results, and you should use an exact test instead (fisher.test()).
  # if (method == 'G-test') {
  #     
  # }
  # 
  # if (method == 't-test') {
  #   
  # }

  
  res$PAdjusted <- p.adjust(res$PValue, method = 'BH', n = length(res$PValue))
  
  res <- res[, c('logFC', 'PValue', 'PAdjusted')]
  
  return(res)
}

# 
# # 
# test_res <- nonpar_test(df_missing_normalized_matrix, 'Ohter', col_annotation$disease, method = 'fisher')
# sum(test_res$PValue < 0.001)
# # # 
# adjust_p_value <- p.adjust(test_res$p.val, method = 'BH', n = length(test_res$PValue))
# sum(adjust_p_value < 0.05)


#' @param res summary table protein expression analysis including logFC, PValue, and PAdjusted
#' @param fold_change absolute log2 fold change threshold
#'
#'
topTable <- function(res, fold_change, PAjusted) {
  filter <- (abs(res$logFC) > fold_change) & (res$PAdjusted < PAjusted)
  return(res[filter, ])
}





#' create data matrix based on conditions to compared and number of valid values per condition
#' 
#' @param X
#' @param col_annotation
#' @param condition1
#' @param condition2
#' @param min_num_valids
#' @return
contrast_dependent_missing_handling <- function(X, col_annotation, condition1, condition2, min_num_valids) {
  col_conditions <- col_annotation[colnames(X), ]
  newX <- X[ , col_conditions == 'Breast' | col_conditions == 'Leukemia']
  # minimum number of valid per condition
  condition1_matrix <- newX[, col_annotation[colnames(newX), ] == condition1]
  condition1_filter <- apply(condition1_matrix, 1, function(c) sum(c!=0, na.rm = T) >= min_num_valids)
  
  condition2_matrix <- newX[, col_annotation[colnames(newX), ] == condition2]
  condition2_filter <- apply(condition2_matrix, 1, function(c) sum(c!=0, na.rm = T) >= min_num_valids)
  newX <- newX[condition1_filter & condition2_filter, ]
  return(newX)
}


