#' @export
#' @title mapping current entry ID to another entry ID based on dictionary in form of `tibble`
#' @description
#' `r lifecycle::badge("stable")`
#' @param keys a string vector contains IDs needed to be mapped from
#' @param dictionary a dataframe in which the rownames are the keys and the first 
#' column contains values that keys will be mapped to in a one-to-one manner. The first column
#' is used as keys and second column will be used as they mapped values. The value in first
#' column must be unique. NA is allow as a key but can only appear once in the first column.
#' @return a string vectors that containes new IDs after mapping
#' @examples
#' dictionary <- tibble::tibble(key=c('a', 'b', 'c'), value=c('A', 'B', 'C'))
#' keys <- c('a', 'a', 'c', 'b')
#' mapping(keys, dictionary)
#' @seealso tibble: https://tibble.tidyverse.org/index.html
mapping <- function(keys, dictionary, na.rm=TRUE) {
  # check whether key-value pair are uniquely defined
  if (any(duplicated(dictionary[[1]]))) {
    rlang::abort('Non-unique keys in dictionary: the same key can be mapped to multiple values')
  }
  keys <- tibble::tibble(key=keys)
  dictionary <- tibble::tibble(key=dictionary[[1]], value=dictionary[[2]])
  values <- dplyr::left_join(keys, dictionary, by='key') %>% dplyr::select('value') %>% dplyr::pull()
  return(values)
}



#' Mapping values in a specific column of a `tibble` based on a dictionary (in form of tibble) and return the same
#' 
#' A helper function to mapping values in a column of `tibble` based on a key-value dictionary simpler. 
#'
#' @param tbl a `tibble` object 
#' @param col a string, the name of the column that whose values will be mapped based a dictionary
#' @param dictionary a dictionary in form of `tibble`, the first column should be `key` and the second column should be `value`
#'
#' @return
#' @export
#'
#' @examples
#' dictionary <- tibble::tibble(key=c('a', 'b', 'c'), value=c('A', 'B', 'C'))
#' A <- tibble::tibble(letter=c('a', 'a', 'c', 'b'), noun=c('apple', 'apple', 'cab', 'bean'))
#' column_mapping(A, 'letter', dictionary)
#' 
#' @seealso tibble: https://tibble.tidyverse.org/index.html
column_mapping <- function(tbl, col, dictionary) {
  values <- mapping(tbl[[{col}]], dictionary)
  tbl[[{col}]] <- values
  return(tbl)
}


# enhanced_mapping <- function(keys, dictionary, map_from, map_to) {
#   dictionary <- dictionary %>% dplyr::select({{map_to}}, {{map_from}})
#   
#   # only keep unique pair
#   if (any(duplicated(dictionary))) {
#     warning('duplicate mapping key-value pairs exists after extracting binary 
#     key-value mapping from the dictionary, 
#             the mapping continues after removing duplicates key-value pairs\n')
#     dictionary <- unique(dictionary)
#   }
#   
#   
#   num_values_per_key <- dictionary %>% dplyr::group_by(.data[[map_from]]) %>% 
#     dplyr::arrange(.data[[map_to]]) %>% 
#     dplyr::summarise(n=n()) %>% select(n) %>% pull()
#   if (any(num_values_per_key > 1)) {
#     warning('one-to-many mapping exists, only take the first pair of that map_from key, please be cautious about the result!\n',
#             'Please check whether the multiple value being matched to one key is caused by NA or missing value for some keys
#             since they are treated as the same\n', 
#             'the number of "map_from" key with multiple "map_to" values is ', sum(num_values_per_key > 1, '\n'))
#     dictionary <- dictionary %>% dplyr::group_by(.data[[map_from]]) %>% 
#       dplyr::arrange(.data[[map_to]]) %>% 
#       dplyr::summarise(!!map_to := first(.data[[map_to]]))
#   }
#   # create the mapping data.frame
#   mapping_df <- data.frame(map_to = dictionary[[map_to]] , row.names = dictionary[[map_from]])
#   warning('If there is error about no duplicate rownames, please make sure the input x is a matrices rather than data.frame\n')
#   rownames(y) <- mapping_df[rownames(y), ]
#   return(y)
# }



#' @title convert row names between different identification commonly used in Omics data
#' @description 
#' `r lifecycle::badge("experimental")`
#' use with caution for when there are multiple-to-multiple mapping exist in the dictionary
#' @param x a numeric matrix, cannot be a data.frame since duplicated row names are not allowed in data.frame
#' @param dictionary
#' @param map_from
#' @param map_to
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
#' Remove proteins isoforms that have a name suffix.
#' 
#' The isoform 
#' For example, protein `P08181-1` is a isoform of protein `P08181`.
#' Be aware that not all isoforms are expressed with a suffix that consists of 
#' hyphen and number: https://www.uniprot.org/help/canonical_and_isoforms.
#' 
#' @param x a numeric matrix, with each row represents a protein and each column represents a sample
#' 
#' 
#' @seealso https://www.uniprot.org/help/canonical_and_isoforms
remove_isoform_suffix <- function(x, method=c('remove', 'combine')) {
  tbl <- x %>% dplyr::as_tibble(rownames = 'accession')
  method <- match.arg(method)
  if (method == 'remove') {
    tbl <- tbl %>% 
      dplyr::mutate('is_isoform'=stringr::str_detect(accession, pattern = '.+-.+')) %>% 
      dplyr::filter(is_isoform == FALSE) %>% 
      dplyr::select(-c('is_isoform')) 
  } else {
    tbl <- tbl %>% dplyr::mutate(accession=stringr::str_extract(accession, pattern='.+(?=-[0-9]+)')) %>% 
      dplyr::group_by(accession) %>% dplyr::summarise_each(sum)
    out <- as.matrix(tbl[-1])
  }
  
  # extract accesson as the new rownames for matrix
  out <- as.matrix(tbl[-1])
  rownames(out) <- tbl[[1]]
  return(out)
}