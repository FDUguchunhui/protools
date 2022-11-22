
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

SpC_List <- function(df, annotation, NA_substitution=NULL, proteins_filter=NULL, replicates_remove=NULL, replicates_keep=NULL) {
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
  
  if (!is.null(replicates_keep)) {
    before <- ncol(SpC_matrix)
    SpC_matrix <- SpC_matrix[ ,colnames(SpC_matrix) %in% replicates_keep]
    annotation <- annotation[(rownames(annotation) %in% replicates_keep), ,drop=FALSE]
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



#' @export
#' @title 
#' @description 
#' @param x a matrix of raw counts with a gene id in as row name
#' @param length a data.frame that contains information about gene length for 
#' @param per_count 
#' @param na_fill the values used to fill if NA is encountered. Can be a scalar
#' for all NA or a vector which values will be sequentially filled into NA positions
#' each gene represented by the gene id in each row 
#' @return the TPM normalized matrix in the same shape as x
#' @references 
#' # https://www.youtube.com/watch?v=TTUrtCY2k-w&t=548s
#  https://www.biostars.org/p/335187/
TPM <- function(x, gene_length, per_count=10e6, na_fill=NULL) {
  if(!is.null(na_fill)) {
    x[is.na(x)] = na_fill
  }
  # reorder the gene id and gene length information as the same in x
  len <-  gene_length[rownames(x), ]
  x <- sweep(x, 1, len, FUN = '/')
  # divide cells in each row by the corresponding total RPK of the replicate it belongs to
  return(t(t(x) * per_count/colSums(x)))
}



#' @export
#' @title Calculate Normalized Spectral Abundance Factor (NASF)
#' @description NSAF_j = (Sc_j/len)/sum(Sc_i/len for all proteins)
#' @param x a matrix of raw counts
#' @param protein_length a data.frame that contains information about protein length for 
#' each protein represented by protein id in each row 
#' @return the NASF normalized matrix in the same shape as x
#' @references 
#' https://github.com/moldach/proteomics-spectralCount-normalization/blob/master/nsaf.R
#' 
#' McIlwain S, Mathews M, Bereman MS, Rubel EW, MacCoss MJ, Noble WS. 
#' Estimating relative abundances of proteins from shotgun proteomics data. 
#' BMC Bioinformatics. 2012 Nov 19;13:308. doi: 10.1186/1471-2105-13-308. PMID: 23164367; PMCID: PMC3599300.
nsaf <- function(x, protein_length, per_count, na_fill=NULL) {
  # reorder the gene id and gene length information as the same in x
  TPM(x, protein_length, per_count = per_count, na_fill = na_fill)
}





