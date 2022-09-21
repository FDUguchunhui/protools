
#' @export
#' @title SpC list constructor
#' @description Create SpC_list object for conviently conducting downstream analysis
#' @param df a data.frame with proteins accession number as row names, replicate (sample) identifier as column names, with
#' numeric data in each cell
#' @param annotation an annotation table provide extra information about the each replicate in the df variable,
#' such as disease type, subtype, and biology processing approach, etc. The rowname of annotation should be replicate identifiers that
#' is also used in the column names of df. Caution: the rownames of annotation are supposed to match with colnames of df.
#' @param NA_substitution how NA value in the df should be replace, can be NULL or a scalar value. If it is NULL
#' then all NA value will be kept as NA
#' @param proteins_filter A String vector, proteins that you want to keep in the data and should be compatible with rownames in df.
#' @param replicates_remove A String vector, replicates that you want to remove 
#' from the data (corresponding replicates will be remove both from df and annotation), 

#' @return a list with two elements: matrix and annotation
#' \itemize{
#'   \item matrix - the data matrix that is converted from the dataframe df using provided parameters
#'   \item annotation - the original annotation passed into the function after removing rownames in replicates_remove 
#' }


#' @example 
#' 
#' 
#' 
#' 

SpC_List <- function(df, annotation, NA_substitution=NULL, proteins_filter=NULL, replicates_remove=NULL) {
  SpC_matrix <- as.matrix(df)
  
  # if need to substitute NA 
  if (!is.null(NA_substitution)) {
    SpC_matrix[is.na(SpC_matrix)] <-  NA_substitution
  }
  
  # if need to filter the proteins included in the data
  if (!is.null(proteins_filter)) {
    before <- nrow(SpC_matrix)
    SpC_matrix <- SpC_matrix[rownames(SpC_matrix) %in% proteins_filter, ]
    after <- nrow(SpC_matrix)
    message('Number of rows before filtering: ', before)
    message('Number of rows after filtering: ', after)
    message()
  }
   
  # if need to remove replicate (sample) from the corresponding column
  if (!is.null(replicates_remove)) {
    before <- ncol(SpC_matrix)
    SpC_matrix <- SpC_matrix[ ,! (colnames(SpC_matrix) %in% replicates_remove)]
    annotation <- annotation[!(rownames(annotation) %in% replicates_remove), ,drop=FALSE]
    after <- ncol(SpC_matrix)
    message('Number of columns before filtering: ', before)
    message('Number of columns after filtering: ', after)
    message()
  }
  
  # create a class "SpC_list"
  
  out <- new('SpC_list', list(matrix=SpC_matrix, annotation=annotation))
  # return_lst[['matrix']] <- SpC_matrix
  # return_lst[['annotation']] <- annotation
  # class(return_lst) <- 'SpC_list'
  return(out)
}