#' @export
#' @title convert row names between different identification commonly used in Omics data
#' @description use with caution for when there are multiple-to-multiple mapping exist in the dictionary
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
            'the number of "map_from" key with multiple "map_to" values is ', sum(num_values_per_key > 1))
    dictionary <- dictionary %>% dplyr::group_by(.data[[map_from]]) %>% 
      dplyr::arrange(.data[[map_to]]) %>% 
      dplyr::summarise(!!map_to := first(.data[[map_to]]))
  }
  # create the mapping data.frame
  mapping_df <- data.frame(map_to = dictionary[[map_to]] , row.names = dictionary[[map_from]])
  rownames(y) <- mapping_df[rownames(y), ]
  return(y)
}


#' @export
#' @title Remove proteins isoforms from the proteomics data matrix
#' @description 
#' @param x a numeric matrix, with each row represents a protein and each column represents a sample
remove_isoforms <- function(x) {
  tbl <- x %>% as_tibble(rownames = 'accession') %>% 
    mutate('is_isoform'=str_detect(accession, pattern = '.+-.+')) %>% 
    filter(is_isoform == FALSE) %>% 
    dplyr::select(-c('is_isoform'))
  
  out <- as.matrix(tbl[-1])
  # need to use pull to make it as vector
  rownames(out) <- pull(tbl[1])
  return(out)
}