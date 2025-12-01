# this is based on salmon quant for subset data used for assembly

library(tidyverse)
library(ggpubr)
library(tximport)
library(cowplot)
library(data.table)

source('~/lab_davidson/yan.a/software/scripts_denovo/R/get_gene_tx.R')

# bambu reference
# tx_counts_bambu <- readRDS('')
# tx_counts_bambu <- tx_counts_bambu[rowSums(tx_counts_bambu$counts) != 0, ]
# dim(tx_counts_bambu)
# 
# gene_counts_bambu <- data.frame(tx_counts_bambu$counts, 
#                                 gene_id = tx_counts_bambu$genes$gene_id) %>%
#   filter(!str_detect(gene_id, '^Bambu')) %>%
#   group_by(gene_id) %>%
#   summarise_all(sum)
# dim(gene_counts_bambu)

gene_counts_bambu <- read.table('../bambu/PB_merged_dge/default/counts_gene.txt', 
                                sep = '\t', header = T)

gene_df_bambu <- data.frame(geneid = gene_counts_bambu$GENEID, 
                            log2(gene_counts_bambu[,-1] + 1), 
                            avgcpm = rowMeans(log2(gene_counts_bambu[,-1] + 1))) 
colnames(gene_df_bambu) <- c('geneid','rep1_10k_counts','rep1_20k_counts','rep2_20k_counts','avgcpm')

truecpm_gene <- gene_df_bambu %>% 
  select(-avgcpm) %>%
  pivot_longer(-geneid, values_to = 'log2CPM') %>%
  unite('rowname', geneid, name) 

# rattle

combos <- data.frame(dir = c(rep(c('../rattle', '../isonform'), each = 4), rep('../rnabloom2', each = 2)),
                 tx2gene = rep(c('../rattle/fx2tab.txt', '../rattle/PB_merged_corset/PB_merged-clusters_mod.txt',
                             '../isonform/fx2tab.txt', '../isonform/PB_merged_corset/PB_merged-clusters_mod.txt',
                             '../rnabloom2/PB_merged_corset/PB_merged-clusters_mod.txt'), each = 2),
                 cluster = rep(c('rattle_rattle', 'rattle_corset', 'isonform_isonclust', 'isonform_corset', 'rnabloom2_corset'), each = 2), 
                 quant = rep(c('ontp','onts'), 5))

cnts.list <- lapply(1:nrow(combos), function(x){
  
  dirs <- list.dirs(combos[x, 1], recursive = F, full.names = T)
  
  files <- list.files(dirs[grepl('dge', dirs)], 
                      pattern = 'quant.sf', full.names = T, recursive = T)
  
  files <- files[grepl(combos[x,4], files)]
  
  method <- str_split(combos[x,3], pattern = '_', simplify = T)[2]
  
  rt_gene <- get_gene_tx(combos[x,2], 
                         list.files(dirs[grepl('sqanti3$', dirs)], pattern = 'classification.txt', full.names = T, recursive = T), 
                         method = method)
  rownames(rt_gene) <- rt_gene$isoform
  
  salmon_gene <- tximport(files, 
                          type = "salmon", 
                          tx2gene = rt_gene[,c('isoform', 'geneid')])
  counts <- salmon_gene$counts
  colnames(counts) <- c('rep1_10k_counts','rep1_20k_counts','rep2_20k_counts')
  counts
})

names(cnts.list) <- combos %>% unite(name, cluster, quant) %>% pull(name)

# use cnts.list to replace x.rds

# for each combo of assember and clustering, quant
# dirs <- list.files('.', 'x.rds', recursive = T) 
# dirs <- dirs[str_detect(dirs, '/x')]

df <- data.frame(path = names(cnts.list)) %>%
  separate(path, c('assembler', 'clustering', 'quant'), remove = F, sep = '_|/') %>%
  unite("cluster_summary", assembler, clustering, sep = '_', remove = F) 

cor <- lapply(1:nrow(df), function(x){
  gene_counts <- cnts.list[[x]]
  summary <- readRDS(paste0(df$cluster_summary[x], '/summary.rds'))
  
  genecount_df <- data.frame(geneid = rownames(gene_counts), 
                             log2(gene_counts + 1), 
                             rowsum = rowSums(log2(gene_counts + 1))) 
  
  cluster_max <- genecount_df %>% 
    left_join(summary$asm_map, by = 'geneid') %>%
    filter(!is.na(assigned_gene)) %>%
    select(-geneid) %>%
    group_by(assigned_gene) %>%
    slice_max(order_by = rowsum, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(-rowsum) %>%
    pivot_longer(1:3, values_to = 'log2CPM') %>%
    unite('rowname', assigned_gene, name) 
  
  exp_df <- left_join(truecpm_gene, 
                      cluster_max, 
                      by = 'rowname', suffix = c('','.asm')) %>%
    replace_na(list(log2CPM.asm = 0)) %>%
    filter(!str_detect(rowname, '^R'))
  
  gene.cor.cluster <- cor(exp_df %>%
                            select(-1))[1,2]
  
  res <- data.frame(gene.cor.cluster = gene.cor.cluster)
  
})

df <- data.frame(df,
                 rbindlist(cor) )
write.csv(df, 'plot/cor_cluster_tiled.csv')



