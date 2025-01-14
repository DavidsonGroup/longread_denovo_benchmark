
get_cluster_exp <- function(x, y, requested_counts = requested_counts) {
  # x salmon quant: 
  # x <- salmon_count$merged_2m_quant_rnabloom2_ontp
  # y summary data from eval_cluster.R containing asm_map and sqanti_anno: 
  # y <- filelist$rnabloom2_corset_2m
  
  ## get cluster counts
  # sum counts when multiple cluster match to same gene
  asmcounts <- x %>% 
    left_join(y$sqanti_anno, by = 'isoform', suffix = c("", ".y")) %>% # in rattle annotation there is rattle count column conflicting
    group_by(geneid) %>%
    summarise(clustercounts = sum(counts)) %>% 
    ungroup() %>%
    left_join(y$asm_map, by = 'geneid') %>% ## obtain cluster counts, but some cluster could match to several genes
    group_by(assigned_gene) %>%
    summarise(assigned_gene_counts = sum(clustercounts)) %>%
    ungroup() 
  
  # or just pick highest count when multiple cluster match to same gene
  asmcounts2 <- x %>% 
    left_join(y$sqanti_anno, by = 'isoform', suffix = c("", ".y")) %>% 
    group_by(geneid) %>%
    summarise(clustercounts = sum(counts)) %>%
    ungroup() %>%
    left_join(y$asm_map, by = 'geneid') %>% ## obtain cluster counts, but some cluster could match to several genes
    group_by(assigned_gene) %>%
    dplyr::filter(clustercounts == max(clustercounts)) %>%
    slice_head(n = 1) %>% # get the first if there is still ties
    ungroup() %>%
    dplyr::rename(assigned_gene_counts = clustercounts) %>%
    dplyr::select(c(3,2)) 
  
  joincounts <- requested_counts %>% 
    filter(!str_detect(gene_id, '^BambuGene')) %>%
    group_by(gene_id) %>%
    summarise(bambucounts = sum(sim_counts)) %>%
    left_join(asmcounts, by = c('gene_id' = 'assigned_gene')) %>%
    left_join(asmcounts2, by = c('gene_id' = 'assigned_gene')) %>%
    replace_na(list(assigned_gene_counts.x = 0, 
                    assigned_gene_counts.y = 0)) 
  
  libsize <- colSums(joincounts[, -1])/1000000
  
  join.cpm <- joincounts %>% 
    mutate(logbambu = log2(bambucounts/libsize[1] + 1), 
           logassigned_sum = log2(assigned_gene_counts.x/libsize[2] + 1),
           logassigned_max = log2(assigned_gene_counts.y/libsize[3] + 1))
  
  # cor(join.cpm[, 5:7])  
  
  return(cluster_exp = join.cpm)
}




