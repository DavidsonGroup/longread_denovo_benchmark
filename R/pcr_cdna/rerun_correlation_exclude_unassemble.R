
library(tidyverse)
library(ggpubr)
library(tximport)
library(data.table)
library(cowplot)

source("~/lab_davidson/yan.a/software/scripts_denovo/R/quantile_overlap.R")

sqanti_summary <- readRDS('plot/sqanti_summary.rds')
requested_counts <- readRDS('plot/requested_counts.rds')

anno.gene <- read.table("~/lab_davidson/yan.a/software/LongReadRNA/sequins/annotations/rnasequin_genes_2.4.tsv",
                        header = TRUE, stringsAsFactors = FALSE) %>%
  mutate(a = log2(MIX_A), b = log2(MIX_B))
# which(!(anno.gene$NAME %in% requested_counts$gene_id))

anno.tx <- read.table("~/lab_davidson/yan.a/software/LongReadRNA/sequins/annotations/rnasequin_isoforms_2.4.tsv", 
                      header = TRUE, stringsAsFactors = FALSE) %>%
  mutate(a = log2(MIX_A), b = log2(MIX_B))
# which(!(anno.tx$NAME %in% requested_counts$transcript_id))

count_add <- data.frame(requested_counts = 0, 
                        sim_counts = 0 , 
                        transcript_id = anno.tx$NAME[which(!(anno.tx$NAME %in% requested_counts$transcript_id))], 
                        length = 0, 
                        gene_id = anno.tx$NAME[which(!(anno.tx$NAME %in% requested_counts$transcript_id))] %>% str_remove('_[0-9]$'), 
                        dtu = NA, dte = NA, truecount_quantile = NA, truecount_quantile_range = NA)

requested_counts <- rbind(requested_counts, count_add)

args <- c('../bambu/',
          '../bambudenovo/',
          '../rattle/', 
          '../rnabloom2/',
          '../isonform/',
          '../trinitystranded/',
          '../rnaspades/',
          '../rnabloom2hybrid/'
)

salmon_files <- lapply(args, function(x) {
  list.files(x, 'quant.sf', recursive = T, full.names = T)
}) %>% unlist()

salmon_files <- salmon_files[str_detect(salmon_files, 'dge')]

df <- data.frame(path = salmon_files) %>%
  separate(path, c(NA, 'assembler',NA, "dep1", "sample",NA), sep = '/', remove = F) %>%
  mutate(assembler = str_remove(assembler, 'stranded')) %>%
  separate(dep1, c(NA, "depth", NA), sep = '_') %>%
  mutate(quant = ifelse(str_detect(sample, 'map|onts'), 'sec', 'pri'),
         depth = ifelse(depth %in% c('ilu','merged'), '10m', depth)) %>%
  unite("sqanti_summary", assembler, depth, sep = '_', remove = F) 

df_list <- split(df, list( df$assembler, df$depth, df$quant), drop = T)

all_quant2 <- lapply(df_list, function(df){
  
  if (nrow(df) == 6) {
    print('Long/short')
    quant_summary <- lapply(1:6, function(x){
      
      counts <- tximport(df[x,'path'], 
                         type = "salmon", 
                         txOut = T)
      
      counts <- data.frame(counts = counts$counts[,1], 
                           isoform = rownames(counts$counts))
      
      quant <- quantile_overlap(sqanti_summary[[df[x,'sqanti_summary']]], 
                                counts, 'isoform', requested_counts)
      
    })
  } else if (nrow(df) == 12) {
    
    print('Hybrid')
    
    quant_summary <- lapply(1:6, function(x){
      # long read
      counts1 <- tximport(df[x,'path'], 
                          type = "salmon", 
                          txOut = T)
      # short read
      counts2 <- tximport(df[x+6,'path'], 
                          type = "salmon", 
                          txOut = T)
      
      stopifnot(all.equal(rownames(counts1$counts), rownames(counts2$counts)))
      
      counts <- data.frame(counts = counts1$counts[,1] + counts2$counts[,1], 
                           isoform = rownames(counts1$counts))
      
      quant <- quantile_overlap(sqanti_summary[[df[x,'sqanti_summary']]], 
                                counts, 'isoform', requested_counts)
      
    })
  }
  
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

saveRDS(all_quant2, 'plot/all_quant2.rds')

# these will be false as each data is filtered by there existing assembly
length(unique(lapply(all_quant2, function(x){
  x$tx_exp_tiled$tx
}))) == 1
length(unique(lapply(all_quant2, function(x){
  x$gene_exp_tiled$gene
}))) == 1

## this will be true
all(sapply(1:14, function(i) {
  identical(all_quant2[[i]]$tx_exp_tiled$tx, 
            all_quant2[[i + 14]]$tx_exp_tiled$tx)
}))
all(sapply(1:14, function(i) {
  identical(all_quant2[[i]]$gene_exp_tiled$gene, 
            all_quant2[[i + 14]]$gene_exp_tiled$gene)
}))

# load bambu reference
tx_counts_bambu <- readRDS('../bambu/merged_10m_dge/R/y.rds')
tx_counts_bambu <- tx_counts_bambu[rowSums(tx_counts_bambu$counts) != 0, ]
dim(tx_counts_bambu)

gene_counts_bambu <- data.frame(tx_counts_bambu$counts, 
                                gene_id = tx_counts_bambu$genes$gene_id) %>%
  filter(!str_detect(gene_id, '^Bambu')) %>%
  group_by(gene_id) %>%
  summarise_all(sum) %>% 
  data.frame()
rownames(gene_counts_bambu) <- gene_counts_bambu$gene_id

# for non sequin
cor_noseq2 <- lapply(all_quant2, function(x){
  
  txid <- unique(x$tx_exp_tiled %>%
                   filter(!str_detect(tx, '^R')) %>%
                   pull(tx))
  
  df_ass <- x$tx_exp_tiled %>%
    filter(tx %in% txid) 
  
  stopifnot(df_ass$txid != rep(txid, 6))
  
  cor.tx <- cor(log2(df_ass$counts.max + 1), 
                log2(tx_counts_bambu$counts + 1)[txid, ] %>%
                  as.matrix() %>%
                  c())
  
  geneid <- unique(x$gene_exp_tiled %>%
                     filter(!str_detect(gene, '^R')) %>%
                     pull(gene))
  
  df_ass2 <- x$gene_exp_tiled %>%
    filter(gene %in% geneid) 
  
  stopifnot(df_ass$gene != rep(geneid, 6))
  
  cor.gene <- cor(log2(df_ass2$assemble_count + 1), 
                  log2(gene_counts_bambu[,-1] + 1)[geneid, ] %>% 
                    as.matrix() %>%
                    c())
  
  res <- data.frame(cor.tx = cor.tx, cor.gene = cor.gene)
  
}) %>% 
  rbindlist() %>%
  mutate(data = names(all_quant2))

# left_join(cor_noseq, cor_noseq2, by= c('data')) %>% 
#   mutate(diff_tx=cor.tx.x-cor.tx.y, 
#          diff_gene= cor.gene.x-cor.gene.y) %>% View

# for sequin
sequin.gene.exp <- anno.gene[, c('a','a','a','b','b','b')] 
rownames(sequin.gene.exp) <- anno.gene$NAME

sequin.tx.exp <- anno.tx[, c('a','a','a','b','b','b')] 
rownames(sequin.tx.exp) <- anno.tx$NAME

cor_seq2 <- lapply(all_quant2, function(x){
  
  txid <- unique(x$tx_exp_tiled %>%
                   filter(str_detect(tx, '^R')) %>%
                   pull(tx))
  print(length(txid))
  
  df_ass <- x$tx_exp_tiled %>%
    filter(tx %in% txid) 
  
  stopifnot(df_ass$txid != rep(txid, 6))
  
  cor.tx <- cor(log2(df_ass$counts.max + 1), 
                sequin.tx.exp[txid, ] %>%
                  as.matrix() %>%
                  c())
  
  geneid <- unique(x$gene_exp_tiled %>%
                     filter(str_detect(gene, '^R')) %>%
                     pull(gene))
  print(length(geneid))
  
  df_ass2 <- x$gene_exp_tiled %>%
    filter(gene %in% geneid) 
  
  stopifnot(df_ass$gene != rep(geneid, 6))
  
  cor.gene <- cor(log2(df_ass2$assemble_count + 1), 
                  sequin.gene.exp[geneid, ] %>% 
                    as.matrix() %>%
                    c())
  
  res <- data.frame(cor.tx = cor.tx, cor.gene = cor.gene)
  
}) %>% 
  rbindlist() %>%
  mutate(data = names(all_quant2))

bambu_dir <- list.dirs('../bambu/')
bambu_dir <- bambu_dir[str_detect(bambu_dir, 'dge/R$')]

cor_seq_bambu2 <- lapply(bambu_dir, function(x){

  tx <- readRDS(paste0(x, '/y.rds'))
  txid <- rownames(tx$counts)[(str_detect(rownames(tx$counts), '^R'))] %>% unique()
  print(length(txid))
  tx <- log2(tx$counts + 1)[txid, ] %>% 
    as.matrix() %>% c()
  cor.tx <- cor(tx, 
                sequin.tx.exp[txid,] %>% 
                  as.matrix() %>%
                  c())  
  
  gene <- readRDS(paste0(x, '/x.rds'))
  geneid <- rownames(gene$counts)[(str_detect(rownames(gene$counts), '^R'))] %>% unique()
  print(length(geneid))
  gene <- log2(gene$counts + 1)[geneid, ] %>% 
    as.matrix() %>% c()
  cor.gene <- cor(gene, 
                  sequin.gene.exp[geneid,] %>% 
                    as.matrix() %>%
                    c())
  
  res <- data.frame(cor.tx = cor.tx, cor.gene = cor.gene)
}) %>%
  rbindlist() %>%
  mutate(data = c('bambu.10m.sec','bambu.2m.sec','bambu.5m.sec'))

cor_seq2 <- rbind(cor_seq_bambu2, cor_seq2)

saveRDS(cor_noseq2, 'plot/cor_noseq2.rds')
saveRDS(cor_seq2, 'plot/cor_seq2.rds')
