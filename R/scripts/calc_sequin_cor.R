calc_sequin_cor <- function(sqanti = sqanti, files = x, 
                            isoform_count = NULL, requested_counts = requested_counts) {
  
  if (is.null(isoform_count)) {
    tmp <-   tximport(files,
                      type = "salmon",
                      txOut = T)
    
    isoform <- data.frame(counts = tmp$counts %>% rowSums(),
                          isoform = rownames(tmp$counts))
  } else {
    isoform <- isoform_count
  }
  
  test <- quantile_overlap(sqanti, isoform, 'test', requested_counts)
  
  # anno and anno2 for sequin gene and tx
  cor.tx <- test$lcpm_tx %>% 
    rownames_to_column %>% filter(str_detect(rowname, '^R')) %>% 
    left_join(anno2, by = c('rowname'='NAME')) %>% 
    mutate(logA = log2(MIX_A), logB = log2(MIX_B)) %>%
    select(# true_cpm, this is bambu quant, ignore
      assemble_cpm.max, assemble_cpm.sum, logA, logB) %>%
    cor
  
  cor.gene <- test$lcpm_gene %>% 
    filter(str_detect(gene_id, '^R')) %>% 
    left_join(anno, by = c('gene_id'='NAME')) %>% 
    mutate(logA = log2(MIX_A), logB = log2(MIX_B)) %>%
    select(# true_cpm, this is bambu quant, ignore
      assemble_cpm, logA, logB) %>%
    cor
  
  return(list(cor.tx = cor.tx[1,], cor.gene = cor.gene[1,]))
  
}
