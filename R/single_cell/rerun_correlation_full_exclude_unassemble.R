# this is based on salmon quant for subset data used for assembly

library(tidyverse)
library(ggpubr)
library(tximport)
library(data.table)
library(cowplot)

source("~/lab_davidson/yan.a/software/scripts_denovo/R/quantile_overlap.R")

sqanti_summary <- readRDS('plot/sqanti_summary.rds')
names(sqanti_summary) <- c('bambu', 'isonform','rattle','rnabloom2')
requested_counts <- readRDS('plot/requested_counts.rds')

dirs <- list.files('.', 'tx_counts_join.csv', recursive = T)

df <- data.frame(path = dirs) %>%
  separate(path, c('assembler', 'clustering', 'quant', NA), remove = F, sep = '_|/') %>%
  unite("cluster_summary", assembler, clustering, sep = '_', remove = F) %>%
  unite("sqanti_summary", assembler, #depth, 
        sep = '_', remove = F) 
df <- df[c(1,3,5,6),]

df_list <- split(df, list( df$assembler, #df$depth, 
                           df$quant), drop = T)
# df_list <- Filter(function(x) nrow(x) > 0, df_list)

all_quant2 <- lapply(df_list, function(df){
  print(df)
  quant_summary <- lapply(1:3, function(x){
    
    counts <- read.csv(df[1,'path'])
    
    if (str_detect(df[1,'path'],'isonform|rnabloom2|rattle')) {
      counts <- data.frame(counts = counts[,2+x], 
                           isoform = str_replace_all(counts$tx, '-', '_'))
    } else {
      counts <- data.frame(counts = counts[,2+x], 
                           isoform = counts$tx)
    }
    
    quant <- quantile_overlap(sqanti_summary[[df[1,'sqanti_summary']]], 
                              counts, 'isoform', requested_counts)
    
  })
  
  stopifnot(length(unique(lapply(quant_summary, function(x) rownames(x$lcpm_tx)))) == 1)
  
  tx_exp_tiled <- lapply(quant_summary, function(y) {
    x <- y$lcpm_tx %>% filter(isassemble != 'not assembled')
    data.frame(tx = rownames(x),
               assemble_cpm = x$assemble_cpm.max, # take max if multiple de novo tx match to 1 ref tx
               counts.max = x$counts.max)
  }) %>% rbindlist()
  
  stopifnot(length(unique(lapply(quant_summary, function(x) rownames(x$lcpm_gene)))) == 1)
  
  gene_exp_tiled <- lapply(quant_summary, function(y) {
    x <- y$lcpm_gene %>% filter(isassemble)
    data.frame(gene = rownames(x),
               assemble_cpm = x$assemble_cpm, # sum up all tx with same sqanti3 associated gene 
               assemble_count = x$counts)
  }) %>% rbindlist()
  
  res <- list(tx_exp_tiled = tx_exp_tiled,
              gene_exp_tiled = gene_exp_tiled)
})

saveRDS(all_quant, 'plot/all_quant_full2.rds')

# these will be false as each data is filtered by there existing assembly
length(unique(lapply(all_quant2, function(x){
  x$tx_exp_tiled$tx
}))) == 1
length(unique(lapply(all_quant2, function(x){
  x$gene_exp_tiled$gene
}))) == 1

# no replicate
## this will be true
# all(sapply(1:3, function(i) {
#   identical(all_quant2[[i]]$tx_exp_tiled$tx, 
#             all_quant2[[i + 3]]$tx_exp_tiled$tx)
# }))
# all(sapply(1:3, function(i) {
#   identical(all_quant2[[i]]$gene_exp_tiled$gene, 
#             all_quant2[[i + 3]]$gene_exp_tiled$gene)
# }))

# load bambu reference

tx_counts_bambu <- read.csv('bambu_bambu_oarfish/tx_counts_join.csv')
rownames(tx_counts_bambu) <- tx_counts_bambu$tx
tx_counts_bambu <- tx_counts_bambu[,3:5]

gene_counts_bambu <- read.csv('bambu_bambu_oarfish/gene_counts_join.csv')
rownames(gene_counts_bambu) <- gene_counts_bambu$gene
gene_counts_bambu <- gene_counts_bambu[,3:5]

# for non sequin
cor_noseq2 <- lapply(all_quant2, function(x){
  
  txid <- unique(x$tx_exp_tiled %>%
                   filter(!str_detect(tx, '^R')) %>%
                   pull(tx))
  
  df_ass <- x$tx_exp_tiled %>%
    filter(tx %in% txid) 
  
  stopifnot(df_ass$txid != rep(txid, 3))
  
  cor.tx <- cor(df_ass%>%
                  pull(assemble_cpm), 
                log2(tx_counts_bambu + 1)[txid, ] %>%
                  as.matrix() %>%
                  c())
  
  geneid <- unique(x$gene_exp_tiled %>%
                     filter(!str_detect(gene, '^R')) %>%
                     pull(gene))
  
  df_ass2 <- x$gene_exp_tiled %>%
    filter(gene %in% geneid) 
  
  stopifnot(df_ass$gene != rep(geneid, 3))
  
  cor.gene <- cor(x$gene_exp_tiled %>%
                    filter(gene %in% geneid) %>%
                    pull(assemble_cpm), 
                  log2(gene_counts_bambu + 1)[geneid, ] %>% 
                    as.matrix() %>%
                    c())
  
  res <- data.frame(cor.tx = cor.tx, cor.gene = cor.gene)
  
}) %>% 
  rbindlist() %>%
  mutate(data = names(all_quant2))

saveRDS(cor_noseq2, 'plot/cor_noseq_full2.rds')


