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
#' @export


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



#' @title Calculate transcription per million (TPM)
#' @description
#' `r lifecycle::badge("stable")`
#'Calculate transcription per million (TPM) based on given information
#' @param x a matrix of raw counts with a gene ID in as rowname
#' @param length a data.frame that contains information about protein length
#' where the first column as the protein ID and the second colmun as the length
#' @param per_count 
#' @param na_fill the values used to fill if NA is encountered. Should be a scalar
#' each gene represented by the gene ID in each row 
#' @return the TPM normalized matrix in the same shape as x
#' @export
#' @references 
#' https://www.youtube.com/watch?v=TTUrtCY2k-w&t=548s
#'  https://www.biostars.org/p/335187/
#'  
#' @example
#' count <-  matrix(
#'   c(10, 12, 30,
#'     20, 25, 60,
#'     5, 8, 15,
#'     0, 0, 1),
#'   nrow = 4,
#'   byrow = T,
#'   dimnames = list(c('A', 'B', 'C', 'D'), c('Rep1', 'Rep2', 'Rep3'))
#' )
#' 
#' gene_length <- data.frame(ID=c('D', 'C', 'B', 'A'), length=c(10, 1, 4, 2))
#' 
#' transcript_per_million_norm(count, gene_length=gene_length, per_count = 10)
#'
#' 
#' 

transcript_per_million_norm <- function(x, gene_length, per_count=10e6, na_fill=NULL) {
  if(!is.null(na_fill)) {
    x[is.na(x)] = na_fill
  }
  # reorder the gene ID and gene length information as the same in x
  protein_length <-  mapping(rownames(x), gene_length)
  # check if any protein in the matrix does have length information
  if (any(is.na(protein_length))) {
    print(rownames(x)[is.na(protein_length)])
    stop('Not all protein have a corresponding length information in "length" file')
  }
  x <- sweep(x, 1, protein_length, FUN = '/')
  # divide cells in each row by the corresponding total RPK of the replicate it belongs to
  return(t(t(x) * per_count/colSums(x)))
}



#' @title Calculate Normalized Spectral Abundance Factor (NASF)
#' @description 
#' `r lifecycle::badge("stable")`
#' \deqn{
#' NSAF_j = \frac{SpC_j/proteinLength_j}{denominator} \\
#' denomiator = \sum{over all protein i} SpC_i/proteinLength_i
#' }
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
#' @export
nsaf <- function(x, protein_length, per_count, na_fill=NULL) {
  # reorder the gene id and gene length information as the same in x
  transcript_per_million_norm(x, protein_length, per_count = per_count, na_fill = na_fill)
}





#' Title Total spectral count normalization
#' 
#' `r lifecycle::badge("stable")`
#' The total spectral count (TSC) normalization involve three steps:
#' 1. Calculate the total number of spectra in each BioSample
#' 2. Calculate the average number of spectra across all BioSamples
#' 3. Multiply each spectrum count in each sample by the average count over the BioSample's total spectrum count
#'
#' @param x a matrix of raw counts with a element ID in as rowname
#'
#' @return a matrix of total spectral count normalized version of original input
#' @export
#'
#' @examples 
#' count <- matrix(
#'   c(12, 8,
#'     6, 3,
#'     4, 3),
#'   nrow = 3,
#'   byrow = T,
#'   dimnames = list(c('A', 'B', 'C'), c('Rep1', 'Rep2'))
#' )
#' total_spectral_count_norm(count)
#' 
#' @references https://support.proteomesoftware.com/hc/en-us/articles/
#' 115002739586-Spectrum-Count-Normalization-in-Scaffold#:~:text=The%
#' 20normalization%20process%20involves%20three,the%20BioSample's%20total%20spectrum%20count
total_spectral_count_norm <- function(x) {
  col_sum <- colSums(x)
  return(sweep(x, 2, mean(col_sum)/col_sum, FUN = '*'))
}
