summarize_by_case <- function(tbl, proteins_subset, negate=FALSE, remove_IPAS=NULL, group_by_IPAS=FALSE, TPM_cutoff = 0, NSAF_cutoff = 0, cases = c(1, 2, 4)) {
  
  if(negate) {
    tbl <- tbl %>% filter(!(accession %in% proteins_subset))
  } else{
    tbl <- tbl %>% filter(accession %in% proteins_subset)
  }
  
  
  if (!is.null(remove_IPAS)) {
    tbl <- tbl %>% filter(!(IPAS %in% remove_IPAS))
  }
  
  tbl <- tbl %>% rowwise() %>%  
    mutate(case=if_else(TPM > TPM_cutoff & NSAF > NSAF_cutoff, 1,
                        if_else(TPM <= TPM_cutoff & NSAF > NSAF_cutoff, 2, 
                                if_else(TPM > TPM_cutoff & NSAF <= NSAF_cutoff, 3, 4)
                        )
    )
    ) %>% 
    filter(case %in% cases) %>% 
    group_by(IPAS, case) %>% 
    summarise(count=n(), ) %>% 
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


b <- summarize_by_case(combined_long, proteins_subset=missing_proteins, group_by_IPAS = FALSE, cases = c(1, 2))