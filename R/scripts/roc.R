## transcript level ROC for DTU-tx and DTE analysis
## transcript assignment based on SQANTI3 assocciated_transcript

dt_roc_all <- function(df = test[[x]], name = name, truth = true_dtutx) {
  
  if ('FDR' %in% colnames(df)) {
    df <- arrange(df, FDR)
  } else if ('adj.P.Val' %in% colnames(df)) {
    df <- arrange(df, adj.P.Val)
  }

  # this keeps all input clusters without removing, but when counting TPs, it only count unique
  df %>% 
    rowid_to_column(var = "rowid") %>%
    mutate(associated_transcript = ifelse(associated_transcript == 'novel', 
                                          paste(name, rowid, sep = '_'),
                                          associated_transcript),
           is.de = as.numeric(associated_transcript %in% truth)) %>% 
    select(associated_transcript, is.de, rowid) %>% 
    # Create a flag for unique DE transcripts
    group_by(associated_transcript) %>%
    mutate(is_first = row_number() == 1 & is.de == 1) %>%
    ungroup() %>%
    # Create cumulative sum of unique DE transcripts
    arrange(row_number()) %>%  # Preserve original order
    mutate(csum = cumsum(is_first),
           method = name)
}

dt_roc_dedup <- function(df = test[[x]], name = name, truth = true_dtutx) {
  
  if ('FDR' %in% colnames(df)) {
    df <- arrange(df, FDR)
  } else if ('adj.P.Val' %in% colnames(df)) {
    df <- arrange(df, adj.P.Val)
  }
  
  # this first remove novel and duplicate rows based on associated_transcript, creating a reduced version
  # then count cumulative TP and FP
  df %>% 
    filter(associated_transcript != 'novel') %>%
    distinct(associated_transcript) %>% 
    rowid_to_column(var = "rowid") %>%
    mutate(tp = as.numeric(associated_transcript %in% truth),
           fp = as.numeric(!(associated_transcript %in% truth))) %>%
    mutate(csum.tp = cumsum(tp),
           csum.fp = cumsum(fp),
           method = name)
}

## gene level ROC for DTU-gene and DGE analysis
## gene assignment based on asm_map from previous analysis of majority vote or Bambu gene_id

dg_roc_all <- function(df = test2[[x]], map = filelist[[x]], name = name, truth = true_dge) {
  
  if ('FDR' %in% colnames(df)) {
    df <- arrange(df, FDR)
  } else if ('adj.P.Val' %in% colnames(df)) {
    df <- arrange(df, adj.P.Val)
  }
  
  # this keeps all input clusters without removing, but when counting TPs, it only count unique
  
  if (!('assigned_gene' %in% colnames(df))) {
    if ('gene_id' %in% colnames(df) ) { # for Bambu and limma
      df <- df %>% mutate(assigned_gene = gene_id) 
    } else if ('geneid' %in% colnames(df)) { # for ams_map
      df <- df %>% left_join(map$asm_map, by= 'geneid')
    }
  }
  
  df %>% 
    rowid_to_column(var = "rowid") %>%
    mutate(is.de = as.numeric(assigned_gene %in% truth)) %>% 
    select(assigned_gene, is.de, rowid) %>% 
    # Create a flag for unique DE transcripts
    group_by(assigned_gene) %>%
    mutate(is_first = row_number() == 1 & is.de == 1) %>%
    ungroup() %>%
    # Create cumulative sum of unique DE transcripts
    arrange(row_number()) %>%  # Preserve original order
    mutate(csum = cumsum(is_first),
           method = name)
}

dg_roc_dedup <- function(df = test2[[x]], map = filelist[[x]], name = name, truth = true_dge) {
  
  if ('FDR' %in% colnames(df)) {
    df <- arrange(df, FDR)
  } else if ('adj.P.Val' %in% colnames(df)) {
    df <- arrange(df, adj.P.Val)
  }
  
  # this first remove novel and duplicate rows based on associated_transcript, creating a reduced version
  # then count cumulative TP and FP
  
  if (!('assigned_gene' %in% colnames(df))) {
    if ('gene_id' %in% colnames(df) ) { # for Bambu and limma
      df <- df %>% mutate(assigned_gene = gene_id) 
    } else if ('geneid' %in% colnames(df)) { # for ams_map
      df <- df %>% left_join(map$asm_map, by= 'geneid')
    }
  }
  
  df %>% 
    filter(!is.na(assigned_gene)) %>%
    distinct(assigned_gene) %>% 
    rowid_to_column(var = "rowid") %>%
    mutate(tp = as.numeric(assigned_gene %in% truth),
           fp = as.numeric(!(assigned_gene %in% truth))) %>%
    mutate(csum.tp = cumsum(tp),
           csum.fp = cumsum(fp),
           method = name)
}


