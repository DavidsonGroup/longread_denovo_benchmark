
library(tidyverse)
library(ggpubr)
library(tximport)
library(cowplot)
library(data.table)

# load sequin
anno.gene <- read.table("~/lab_davidson/yan.a/software/LongReadRNA/sequins/annotations/rnasequin_genes_2.4.tsv",
                        header = TRUE, stringsAsFactors = FALSE) %>%
  mutate(a = log2(MIX_A), b = log2(MIX_B))
# which(!(anno.gene$NAME %in% requested_counts$gene_id))

anno.tx <- read.table("~/lab_davidson/yan.a/software/LongReadRNA/sequins/annotations/rnasequin_isoforms_2.4.tsv", 
                      header = TRUE, stringsAsFactors = FALSE) %>%
  mutate(a = log2(MIX_A), b = log2(MIX_B))

sequin.gene.exp <- anno.gene[, c('NAME','a','a','a','b','b','b')] 
rownames(sequin.gene.exp) <- anno.gene$NAME
colnames(sequin.gene.exp) <- c('gene', 'ctrl1','ctrl2','ctrl3','de1','de2','de3')
sequin.gene.exp <- sequin.gene.exp %>%
  pivot_longer(-gene, values_to = 'log2CPM') %>%
  unite('rowname', gene, name) 

sequin.tx.exp <- anno.tx[, c('NAME','a','a','a','b','b','b')] 
rownames(sequin.tx.exp) <- anno.tx$NAME
colnames(sequin.tx.exp) <- c('tx', 'ctrl1','ctrl2','ctrl3','de1','de2','de3')
sequin.tx.exp <- sequin.tx.exp %>%
  pivot_longer(-tx, values_to = 'log2CPM') %>%
  unite('rowname', tx, name) 


# bambu reference
tx_counts_bambu <- readRDS('../bambu/merged_10m_dge/R/y.rds')
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
dirs <- c(dirs, 
          list.files('../bambu/', 'x.rds', recursive = T, full.names = T))
dirs <- dirs[str_detect(dirs, '/x')]

df <- data.frame(path = dirs) 

df$path2 <- df$path  
df$path2[(nrow(df)-2):nrow(df)] <- c('bambu_bambu_10m_count/x.rds', 'bambu_bambu_2m_count/x.rds', 'bambu_bambu_5m_count/x.rds')

df <- df %>%
  separate(path2, c('assembler', 'clustering', 'depth', 'quant', NA), remove = F, sep = '_|/') %>%
  unite("cluster_summary", assembler, clustering, depth, sep = '_', remove = F) 

cor <- lapply(1:nrow(df), function(x){
  gene_counts <- readRDS(df$path[x])
  colnames(gene_counts) <- c('ctrl1','ctrl2','ctrl3','de1','de2','de3')
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
  
  exp_df2 <- left_join(sequin.gene.exp, 
                      cluster_max, 
                      by = 'rowname', suffix = c('','.asm')) %>%
    replace_na(list(log2CPM.asm = 0)) %>%
    filter(str_detect(rowname, '^R'))
  
  gene.cor.cluster.sequin <- cor(exp_df2 %>%
                            select(-1))[1,2]
  
  res <- data.frame(gene.cor.cluster = gene.cor.cluster, gene.cor.cluster.sequin = gene.cor.cluster.sequin)
  
})

df1 <- data.frame(df,
                 rbindlist(cor) )
write.csv(df1, 'plot/cor_cluster_tiled.csv')

## remove unassembled?
cor2 <- lapply(1:nrow(df), function(x){
  gene_counts <- readRDS(df$path[x])
  colnames(gene_counts) <- c('ctrl1','ctrl2','ctrl3','de1','de2','de3')
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
    # replace_na(list(log2CPM.asm = 0)) %>%
    filter(!str_detect(rowname, '^R'), !is.na(log2CPM.asm))
  
  gene.cor.cluster <- cor(exp_df %>%
                            select(-1))[1,2]
  
  exp_df2 <- left_join(sequin.gene.exp, 
                       cluster_max, 
                       by = 'rowname', suffix = c('','.asm')) %>%
    #replace_na(list(log2CPM.asm = 0)) %>%
    filter(str_detect(rowname, '^R'), !is.na(log2CPM.asm))
  
  gene.cor.cluster.sequin <- cor(exp_df2 %>%
                                   select(-1))[1,2]
  
  res <- data.frame(gene.cor.cluster = gene.cor.cluster, gene.cor.cluster.sequin = gene.cor.cluster.sequin)
  
})

df2 <- data.frame(df,
                 rbindlist(cor2) )
write.csv(df2, 'plot/cor_cluster_tiled2.csv')

  