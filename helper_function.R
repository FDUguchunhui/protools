combine_isoforms <- function(mat) {
  tbl <- mat %>% as_tibble(rownames = 'accession') %>% 
    mutate(accession=str_extract(accession, pattern = '[^-]+')) %>%
    group_by(accession) %>% 
    summarise_all(sum)
  
  out <- as.matrix(tbl[-1])
  # need to use pull to make it as vector
  rownames(out) <- pull(tbl[1])
}