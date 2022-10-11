#' @export
#' @title convert row names between different identification commonly used in Omics data
#' @description use with caution for when there are multiple-to-multiple mapping exist in the dictionary
#' @param x
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
  if (any(duplicated(dictionary[[map_to]]))) {
    warning('one-to-many mapping exists, only take the first pair of that map_from key, please be cautious about the result!\n
            Please check whether the multiple value being matched to one key is caused by NA or missing value for some keys
            since they are treated as the same')
    dictionary <- dictionary %>% dplyr::group_by(.data[[map_from]]) %>% 
      dplyr::arrange(.data[[map_to]]) %>% 
      dplyr::summarise(to=first(.data[[map_to]]))
  }
  # create the mapping data.frame
  mapping_df <- data.frame(map_to = dictionary$to , row.names = dictionary[[map_from]])
  rownames(y) <- mapping_df[rownames(y), ]
  return(y)
}