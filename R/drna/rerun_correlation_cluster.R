
library(tidyverse)
library(ggpubr)
library(tximport)
library(cowplot)
library(data.table)

# bambu reference
tx_counts_bambu <- readRDS('../bambu/A549_MCF7_directRNA_dge/R/y.rds')
tx_counts_bambu <- tx_counts_bambu[rowSums(tx_counts_bambu$counts) != 0, ]
dim(tx_counts_bambu)

gene_counts_bambu <- data.frame(tx_counts_bambu$counts, 
                                gene_id = tx_counts_bambu$genes$gene_id) %>%
  filter(!str_detect(gene_id, '^Bambu')) %>%
  group_by(gene_id) %>%
  summarise_all(sum)
dim(gene_counts_bambu)

gene_df_bambu <- data.frame(geneid = gene_counts_bambu$gene_id, 
                            log2(gene_counts_bambu[,-1] + 1), 
                            avgcpm = rowMeans(log2(gene_counts_bambu[,-1] + 1))) 
colnames(gene_df_bambu) <- c('geneid','ctrl1','ctrl2','ctrl3','de1','de2','de3','avgcpm')

truecpm_gene <- gene_df_bambu %>% 
  select(-avgcpm) %>%
  pivot_longer(-geneid, values_to = 'log2CPM') %>%
  unite('rowname', geneid, name) 

# for each combo of assember and clustering, quant
dirs <- list.files('.', 'x.rds', recursive = T) 
dirs <- dirs[str_detect(dirs, '/x')]

df <- data.frame(path = dirs) %>%
  separate(path, c('assembler', 'clustering', 'quant', NA), remove = F, sep = '_|/') %>%
  unite("cluster_summary", assembler, clustering, sep = '_', remove = F) 

cor <- lapply(1:nrow(df), function(x){
  gene_counts <- readRDS(df$path[x])
  summary <- readRDS(paste0(df$cluster_summary[x], '/summary.rds'))
  
  genecount_df <- data.frame(geneid = rownames(gene_counts$counts), 
                             log2(gene_counts$counts + 1), 
                             rowsum = rowSums(log2(gene_counts$counts + 1))) 
  
  cluster_max <- genecount_df %>% 
    left_join(summary$asm_map, by = 'geneid') %>%
    filter(!is.na(assigned_gene)) %>%
    select(-geneid) %>%
    group_by(assigned_gene) %>%
    slice_max(order_by = rowsum, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(-rowsum) %>%
    pivot_longer(1:6, values_to = 'log2CPM') %>%
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


  