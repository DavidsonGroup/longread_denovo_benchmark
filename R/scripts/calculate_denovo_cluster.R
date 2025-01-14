calculate_denovo_cluster <- function(cluster_summary = filelist$isonform_corset, sequin = F) {
  
  subset <- cluster_summary$sqanti_anno %>%
    mutate(gene_category = ifelse(structural_category %in% c('antisense','fusion','intergenic','genic_intron', NA), 
                                  structural_category,
                                  'known_gene')) %>%
    mutate(gene_category = ifelse(is.na(gene_category), 'not_mapped', gene_category)) %>%
    dplyr::filter(gene_category=='known_gene' ) ## any de novo transcript mapped to known gene (FSM/ISM, NNC/NIC, genic)
  
  nclst_denovo <- length(unique(cluster_summary$sqanti_anno$geneid)) # same number of cluster in the *mod.txt corset output, and same as asm_map
  nclst_denovo_known <- length(unique(subset$geneid)) 
  nclst_denovo_na <- nclst_denovo - nclst_denovo_known
  
  clst_count <- subset %>% 
    group_by(geneid) %>% 
    summarise(Unique_Count = n_distinct(associated_gene), # how many unique associated genes in cluster
              num_denovo_tx = dplyr::n()) %>%
    ungroup()
  
  if (sequin) {
    asm <- cluster_summary$asm_map %>%
      filter(str_detect(assigned_gene, '^R'))
  } else {
    asm <- cluster_summary$asm_map
    stopifnot(all.equal(dim(asm)[1], nclst_denovo))
    stopifnot(all.equal(dim(asm %>% filter(is.na(assigned_gene)))[1], nclst_denovo_na))
  }

  # asm_combined <- asm %>% 
  #   left_join(clst_count, by = 'geneid') %>% 
  #   group_by(assigned_gene) %>%
  #   arrange(desc(num_denovo_tx)) %>%
  #   mutate(dup = duplicated(assigned_gene)) %>%
  #   mutate(category = ifelse(is.na(Unique_Count), 'no_match', 
  #                            ifelse(Unique_Count == 1 & !dup, 'perfect_match', 
  #                                   ifelse(Unique_Count == 1 & dup, 'redundant_match', 'mixed')))
  #   ) %>% 
  #   mutate(category = factor(category, levels = c('mixed', 'no_match', 'perfect_match', 'redundant_match'))) %>%
  #   ungroup()
  
  # asm_combined <- asm %>% 
  #   left_join(clst_count, by = 'geneid') %>% 
  #   group_by(assigned_gene) %>% 
  #   arrange(desc(num_denovo_tx)) %>%
  #   add_count() %>%
  #   mutate(category = case_when(
  #     is.na(assigned_gene) ~ 'no_match',
  #     Unique_Count == 1 ~ 'match',
  #     Unique_Count > 1 ~ 'mixed',
  #     TRUE ~ NA_character_
  #   )) %>%
  #   # First filter out mixed entries for ranking
  #   group_modify(function(group_data, key) {
  #     # Create rank only considering non-mixed rows
  #     group_data %>%
  #       filter(category != 'mixed') %>%
  #       mutate(filtered_rank = row_number()) %>%
  #       right_join(group_data, by = colnames(group_data))
  #   }) %>%
  #   mutate(
  #     match_type = case_when(
  #       # Perfect match if category is match and it's the first non-mixed row
  #       category == 'match' & filtered_rank == 1 ~ 'perfect_match',
  #       # Redundant match for other matches
  #       category == 'match' & filtered_rank > 1 ~ 'redundant_match',
  #       TRUE ~ NA_character_
  #     )
  #   ) %>%
  #   ungroup()
  
  setDT(asm)
  setDT(clst_count)
  
  asm_combined <- merge(asm, clst_count, 
                        by = "geneid", 
                        all.x = TRUE)
  
  # Convert to data.table if not already
  setDT(asm_combined)
  
  # Sort by num_denovo_tx within each group
  setorder(asm_combined, assigned_gene, -num_denovo_tx)
  
  # Add count by group
  asm_combined[, n := .N, by = assigned_gene]
  
  # Add category
  asm_combined[, category := fcase(
    is.na(assigned_gene), "no_match",
    Unique_Count == 1, "match",
    Unique_Count > 1, "mixed",
    default = NA_character_
  )]
  
  # Create filtered rank
  asm_combined[, filtered_rank := {
    temp_rank <- rep(NA_integer_, .N)
    non_mixed_idx <- category != "mixed"
    if(any(non_mixed_idx)) {
      temp_rank[non_mixed_idx] <- seq_len(sum(non_mixed_idx))
    }
    temp_rank
  }, by = assigned_gene]
  
  # Add match_type
  asm_combined[, match_type := fcase(
    category == "match" & filtered_rank == 1, "perfect_match",
    category == "match" & filtered_rank > 1, "redundant_match",
    default = NA_character_
  )]
  
  cnts <- asm_combined %>% 
    mutate(category = ifelse(category == 'match', match_type, category)) %>%
    mutate(category = factor(category, levels = c('mixed', 'no_match', 'perfect_match', 'redundant_match'))) %>%
    dplyr::count(category, .drop = F)
  
  # stopifnot(all.equal(cnts %>% filter(category %in% c('perfect_match', 'redundant_match')) %>% pull(n) %>% sum(), 
  #           cluster_summary$pure_denovo %>% filter(Unique_Count==1) %>% pull(n)))
  
  return(list(asm = asm_combined, counts = cnts))
  
}



