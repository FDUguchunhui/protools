#' Create table of differentially expressed tags using nonparametric method
#' @description Perform non-parametric test for sequencing data
#' 
#' @param x a object of SpC_list class created by SpC_list function
#' @param case String, the case condition
#' @param method the method used for statistical testing
#' @param treat_zero_as_missing whether treat 0 value as NA, this option only
#' works when method='Wilcoxon'
#'
#' @details if method='fisher' is used. 
#'   # collapse the matrix into a 2 by 2 table with 

#' |     | case | other    |
#'  | :---        |    :----:   |          ---: |
#'  | gene X      | a     | b |
#' | Other genes   | c        | c      |
#'  
#'  The Fisher's exact test is applied to each row in the y$matrix
#'  
#'  if method='Wilcoxon' is used, a two-sided Wilxoxon-sum-rank test will be used,
#'  the p-value is calculated using approximation because it is very likely there are
#'  ties in the hiigh-throughput data especially there are a lot of 0's
#'
#' @return a data frame with row name the same as the data_matrix, and p-value,
#' adjusted p-value using BH method
#' \itemize{
#'   \item logFC - log2 fold change 
#'   \item PValue - p-value calculated based on selected statistical method
#'   \item PAdjusted - adjusted p-value using BH method
#' }
#' @author Chunhui Gu
#' @references Ming Li (2010), Comparative Shotgun Proteomics Using Spectral Count Data and Quasi-Likelihood Modeling
nonpar_test <- function(x, case, method=c('fisher', 'Wilcoxon'), treat_zero_as_missing=TRUE) {
  # collapse all other group category make it a two group design
  data_matrix <- x$matrix
  groups <- ifelse(x$annotation$disease==case, case, 'control group')
  design <- model.matrix(~0+groups)
  
  
  res <- data.frame(PValue = rep(NA, nrow(data_matrix)))
  rownames(res) <- rownames(data_matrix)
  
  
  
  method <- match.arg(method)
  if (method == 'fisher') {
    
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
    if (treat_zero_as_missing) {
      data_matrix[which(data_matrix==0)] <- NA
    }
    # each row of x is the 
    x <- data_matrix[, design[ ,1] == 1]
    y <- data_matrix[, design[ ,2] == 1]
    for (rowname in rownames(data_matrix)) {
      x_current_protein <- x[rowname, ]
      y_current_protein <- y[rowname, ]
      if (all(is.na(x_current_protein)) || all(is.na(y_current_protein))) {
        res[rowname, 'PValue'] <- NA
        res[rowname, 'logFC'] <- NA
      } else {
        res[rowname, 'PValue'] <- wilcox.test(x[rowname, ], y[rowname, ], exact=FALSE)$p.value
        res[rowname, 'logFC']  <- mean(x_current_protein, na.rm = TRUE)/mean(y_current_protein, na.rm = TRUE) 
      }
    }
  } else {
    stop('method is not "fisher" or "Wilcoxon"')
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



#' Table of the Top Differentially Expressed tags
#' 
#' @description filter the top differentially expressed tags based on adjustment p-value
#' and fold-change
#' 
#' @param res summary table protein expression analysis including logFC, PValue, and PAdjusted
#' @param fold_change absolute log2 fold change threshold
#' @param PAdjusted filter tags that less than the adjusted p-value used
#'
#'
topTable <- function(res, logFC, PAjusted) {
  row_index_by_logFC <- which((abs(res$logFC) > logFC))
  row_index_by_PAjuted <- which((res$PAdjusted < PAjusted))
  filter <- intersect(row_index_by_logFC , row_index_by_PAjuted )
  return(res[filter, ])
}



#' @title One-versus-all test
#' @description this function perform a series of one-versus-all tests with respect to every category in
#' the y$annotation$disease
#' @param y a object of SpC_list class which contains a data matrix and annotation table for replicates. It is created by SpC_list function
#' @param logFC the log-fold-change that the returned proteins must be large than for each of the one-versus-all tests
#' @param PAjusted the adjusted p-value that the returned proteins must be large than for each of the one-versus-all tests
#' @param ... additional parameter passed into nonpar_test function
#' 
#' @return a list, each element in the list is a data.frame from a one-versus-all test result. The name shows what is the case used in the test.
#' For example, list[['diseaseX']] is the test result of diseaseX-versus-the-other-diseases
#' @author Chunhui Gu
one_vs_all_test <- function(x,  method = c('fisher', 'Wilcoxon'), logFC = 1, PAjusted = 0.05, ...) {
  method <- match.arg(method)
  
  res_lst <- list()
  diseases <- unique(x$annotation$disease)
  for(case in diseases) {
    res <- nonpar_test(x, case = case, method = method)
    top_res <- topTable(res, logFC = logFC, PAjusted = PAjusted)
    res_lst[[case]] <- top_res
  }
  return(res_lst)
}










#' create data matrix based on conditions to compared and number of valid values per condition
#' 
#' @param X x
#' @param col_annotation x
#' @param condition1 x
#' @param condition2 x
#' @param min_num_valids x
#' @return x
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


