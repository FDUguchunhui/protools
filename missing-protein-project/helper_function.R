combine_isoforms <- function(mat) {
  tbl <- mat %>% as_tibble(rownames = 'accession') %>% 
    mutate(accession=str_extract(accession, pattern = '[^-]+')) %>%
    group_by(accession) %>% 
    summarise_all(sum)
  
  out <- as.matrix(tbl[-1])
  # need to use pull to make it as vector
  rownames(out) <- pull(tbl[1])
}


summarize_by_case <- function(tbl, missing_protein=FALSE, remove_IPAS=NULL, group_by_IPAS=FALSE, TPM_cutoff = 0, NSAF_cutoff = 0, cases = c(1, 2, 4)) {
  
  if (missing_protein) {
    tbl <- tbl %>% filter(accession %in% missing_protein_df$`gene name(s)`)
  }
  
  if (!is.null(remove_IPAS)) {
    tbl <- tbl %>% filter(!(IPAS==remove_IPAS))
  }
  
  tbl <- tbl %>% rowwise() %>%  
    mutate(case=if_else(TPM > TPM_cutoff & NSAF > NSAF_cutoff, 1,
                        if_else(TPM <= TPM_cutoff & NSAF > NSAF_cutoff, 2, 
                                if_else(TPM > TPM_cutoff & NSAF <= NSAF_cutoff, 4, 3)
                        )
    )
    ) %>% 
    filter(case %in% cases) %>% 
    group_by(IPAS, case) %>% 
    summarise(count=n()) %>% 
    ungroup(IPAS, case)
  
  if (group_by_IPAS) {
    tbl <- tbl %>% group_by(IPAS)
  }
  
  tbl <- tbl %>%
    group_by(case, .add=TRUE) %>% 
    summarise(count= sum(count)) %>%
    mutate(countT=sum(count)) %>% 
    mutate(per=paste0(round(100*count/countT,2),'%'))
  
  return(tbl)
}