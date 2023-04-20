
#' @title mapping current entry ID to another entry ID based on dictionary in form of `tibble` or `data.frame`
#' @description
#' `r lifecycle::badge("stable")`
#' @param keys a character vector contains IDs needed to be mapped from
#' @param dictionary a tibble or data.frame in which the first column are the keys and the
#'  second column will be used as they mapped values. The value in first
#' column must be unique. NA is allowed as a key but can only appear once in the first column.
#' @return a string vectors that contains new IDs after mapping
#' @examples
#' dictionary <- tibble::tibble(key=c('a', 'b', 'c'), value=c('A', 'B', 'C'))
#' keys <- c('a', 'a', 'c', 'b')
#' mapping(keys, dictionary)
#' @seealso 
#' tibble: https://tibble.tidyverse.org/index.html
#' @export
mapping <- function(keys, dictionary, na.rm=TRUE) {
  if(!tibble::is_tibble(dictionary)) {
    stop('the dictionary used is not a tibble')
  }
  
  
  
  # check whether key-value pair are uniquely defined
  # [[ return a vector while [ return a tibble
  dictionary_keys <- dictionary[[1]]
  if (any(duplicated(dictionary_keys))) {
    stop('Non-unique keys in dictionary: the same key can be mapped to multiple values')
  }
  values <- dictionary[match(keys, dictionary_keys), ][[2]]
  return(values)
}


#' mapping current entry ID in a `tibble` or `data.frame` column to another 
#' entry ID based on dictionary in form of `tibble` or `data.frame`
#' 
#' Mapping values in a specific column (as keys) in a `tibble` or `data.frame` to another value based on a
#' a given dictionary (in form of `tibble` or `data.frame`) and return the mapped values
#' in the same order as in the keys
#' 
#' A wrapper function to mapping values in a column of `tibble` based on a key-value dictionary simpler. 
#'
#' @param tbl a `tibble` object 
#' @param col a string, the name of the column that whose values will be mapped based a dictionary
#' @param dictionary a dictionary in form of `tibble`, the first column should be `key` and the second column should be `value`
#'
#' @return a string vectors that contains new IDs after mapping
#' @export
#'
#' @examples
#' dictionary <- tibble::tibble(key=c('a', 'b', 'c'), value=c('A', 'B', 'C'))
#' A <- tibble::tibble(letter=c('a', 'a', 'c', 'b'), noun=c('apple', 'apple', 'cab', 'bean'))
#' column_mapping(A, 'letter', dictionary)
#' 
#' @seealso 
#' mapping
#' 
#' tibble: https://tibble.tidyverse.org/index.html
column_mapping <- function(tbl, col, dictionary) {
  values <- mapping(tbl[[{col}]], dictionary)
  tbl[[{col}]] <- values
  return(tbl)
}


#' @title convert rownames between different ID commonly used in Omics data
#' `r lifecycle::badge("experimental")`
#' @description  This is a conceived "one-to-all" solution to map between different ID based on
#' a single dictionary file. However, due to the complex cases of mapping between
#' different ID system, such as one-to-many mappin, the mapping dictionary could be
#' extremely complex when there are more than two ID systems.
#' Please use with caution for when there are multiple-to-multiple mapping exist in the dictionary.
#' This function is not recommended to use, and we recommend to create a clear binary 
#' mapping dictionary and then use protools::mapping instead
#' 
#' 
#' 
#' @param x a numeric matrix, cannot be a data.frame since duplicated row names are not allowed in data.frame
#' @param dictionary a `tibble` or `data.frame` with columns as ID from different naming systems
#' @param map_from a character for specifying which column to use as keys
#' @param map_to a character for specifying which column to use as values
rownames_mapping <- function(x, dictionary, map_from, map_to) {
  # dictionary <- dplyr::as.tibble(dictionary)
  # dictionary <- dictionary %>% dplyr::select(map_from, map_to) %>% 
  #   dplyr::group_by(map_from) %>% 
  #   dplyr::arrange(map_to) %>% 
  #   dplyr::filter(row_number() == 1)
  dictionary <- dictionary %>% dplyr::select(map_to, map_from)
  # copy to a new variable, y and x will point to different object after modification on y later
  y <- x
  # clean any NA map_from 
  dictionary <-  dictionary[!is.na(dictionary[[map_from]]) & !is.na(dictionary[[map_to]]), ]
  # only keep unique pair
  if (any(duplicated(dictionary))) {
    warning('duplicate mapping key-value pairs exists, the mapping continues after removing duplicates key-value pairs\n')
    dictionary <- unique(dictionary)
  }
  
  num_values_per_key <- dictionary %>% dplyr::group_by(.data[[map_from]]) %>% 
    dplyr::arrange(.data[[map_to]]) %>% 
    dplyr::summarise(n=n()) %>% select(n) %>% pull()
  if (any(num_values_per_key > 1)) {
    warning('one-to-many mapping exists, only take the first pair of that map_from key, please be cautious about the result!\n',
            'Please check whether the multiple value being matched to one key is caused by NA or missing value for some keys
            since they are treated as the same\n', 
            'the number of "map_from" key with multiple "map_to" values is ', sum(num_values_per_key > 1), '\n')
    dictionary <- dictionary %>% dplyr::group_by(.data[[map_from]]) %>% 
      dplyr::arrange(.data[[map_to]]) %>% 
      dplyr::summarise(!!map_to := first(.data[[map_to]]))
  }
  # create the mapping data.frame
  mapping_df <- data.frame(map_to = dictionary[[map_to]] , row.names = dictionary[[map_from]])
  warning('If there is error about no duplicate rownames, please make sure the input x is a matrices rather than data.frame\n')
  rownames(y) <- mapping_df[rownames(y), ]
  return(y)
}



#' @export
#' @title Remove protein isoforms that have a name suffix.
#' 
#' @description  example, protein `P08181-1` is an isoform of protein `P08181`.
#' Be aware that not all isoforms are expressed with a suffix that consists of 
#' hyphen and number: https://www.uniprot.org/help/canonical_and_isoforms.
#' 
#' @param x a numeric matrix, with each row represents a protein and each column represents a sample
#' @param method how data is processed in the processing of removing isoforms.
#' Use `remove` if you want just discard isoform data. Use `combine` to sum data of isoform into
#' non-isoform entries.
#' @examples
#' # when use `method='remove'`
#' dat <- data.frame(a=c(1, 2, 3, 4, 5), b=c(1, 2, 3, 4, 5), row.names = c('P1001', 'P1001-1', 'P1002', 'P1002-2', 'P1003'))
#' mat <- as.matrix(dat)
#' remove_isoform_with_suffix(mat, method='remove')
#' 
#' # When use `method='combine'`
#' remove_isoform_with_suffix(mat, method='combine')
#' @seealso https://www.uniprot.org/help/canonical_and_isoforms
remove_isoform_with_suffix <- function(x, method=c('remove', 'combine')) {
  tbl <- x %>% dplyr::as_tibble(rownames = 'accession')
  method <- match.arg(method)
  if (method == 'remove') {
    tbl <- tbl %>% 
      dplyr::mutate('is_isoform'=stringr::str_detect(accession, pattern = '.+-.+')) %>% 
      dplyr::filter(is_isoform == FALSE) %>% 
      dplyr::select(-c('is_isoform')) 
  } else {
    tbl <- tbl %>% dplyr::mutate(accession=stringr::str_extract(accession, pattern='.+(?=-[0-9]+)')) %>% 
      dplyr::group_by(accession) %>% dplyr::summarise(across(.cols=everything(), .fns=sum))
    out <- as.matrix(tbl[-1])
  }
  
  # extract accession as the new rownames for matrix
  out <- as.matrix(tbl[-1])
  rownames(out) <- tbl[[1]]
  return(out)
}