# this is based on oarfish quant for complete data

library(tidyverse)
library(ggpubr)
library(tximport)
library(cowplot)
library(data.table)

# bambu reference
# tx_counts_bambu <- readRDS('../bambu/A549_MCF7_directRNA_dge/R/y.rds')
# tx_counts_bambu <- tx_counts_bambu[rowSums(tx_counts_bambu$counts) != 0, ]
# dim(tx_counts_bambu)
# tx_counts_bambu <- read.csv('bambu_bambu_oarfish/tx_counts_join.csv')
# rownames(tx_counts_bambu) <- tx_counts_bambu$tx
# tx_counts_bambu <- tx_counts_bambu[,3:5]
# tx_counts_bambu <- tx_counts_bambu[rowSums(tx_counts_bambu) != 0, ]
# dim(tx_counts_bambu)

# gene_counts_bambu <- data.frame(tx_counts_bambu$counts, 
#                                 gene_id = tx_counts_bambu$genes$gene_id) %>%
#   filter(!str_detect(gene_id, '^Bambu')) %>%
#   group_by(gene_id) %>%
#   summarise_all(sum)
# dim(gene_counts_bambu)
gene_counts_bambu <- read.csv('bambu_bambu_oarfish/gene_counts_join.csv')
rownames(gene_counts_bambu) <- gene_counts_bambu$gene
gene_counts_bambu <- gene_counts_bambu[,3:5]

gene_df_bambu <- data.frame(geneid = rownames(gene_counts_bambu), 
                            log2(gene_counts_bambu + 1), 
                            avgcpm = rowMeans(log2(gene_counts_bambu + 1))) 
colnames(gene_df_bambu) <- c('geneid','rep1_10k_counts','rep1_20k_counts','rep2_20k_counts','avgcpm')

truecpm_gene <- gene_df_bambu %>% 
  select(-avgcpm) %>%
  pivot_longer(-geneid, values_to = 'log2CPM') %>%
  unite('rowname', geneid, name) 

# for each combo of assember and clustering, quant
# dirs <- list.files('.', 'x.rds', recursive = T) 
# dirs <- dirs[str_detect(dirs, '/x')]

dirs <- list.files('.', 'gene_counts_join.csv', recursive = T)

df <- data.frame(path = dirs) %>%
  separate(path, c('assembler', 'clustering', 'quant', NA), remove = F, sep = '_|/') %>%
  unite("cluster_summary", assembler, clustering, sep = '_', remove = F) 

cor <- lapply(1:nrow(df), function(x){
  # gene_counts <- readRDS(df$path[x])
  gene_counts <- read.csv(df$path[x])
  counts <- gene_counts[, 3:5]

  summary <- readRDS(paste0(df$cluster_summary[x], '/summary.rds'))
  
  genecount_df <- data.frame(geneid = gene_counts$gene, 
                             log2(counts + 1), 
                             rowsum = rowSums(log2(counts + 1))) 
  
  if (df$clustering[x] != 'rattle') {
    genecount_df <- genecount_df %>% 
      mutate(geneid = str_replace_all(geneid, '_', '-'))
  }
  
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
    filter(!str_detect(rowname, '^R'))

})

cor1 <- sapply(cor, function(x) {
  df <- x %>%
    replace_na(list(log2CPM.asm = 0)) %>%
    select(-1)
  cor(df)[1,2]
}, simplify = T) 

df1 <- data.frame(df, 
                  gene.cor.cluster = cor1)
write.csv(df1, 'plot/cor_cluster_tiled_full.csv')

cor2 <- sapply(cor, function(x) {
  df <- x %>%
    filter(!is.na(log2CPM.asm)) %>%
    select(-1)
  cor(df)[1,2]
}, simplify = T) 

df2 <- data.frame(df,
                  gene.cor.cluster = cor2)
write.csv(df2, 'plot/cor_cluster_tiled_full2.csv')


  