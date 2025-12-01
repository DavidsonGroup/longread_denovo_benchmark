# this is based on salmon quant for subset data used for assembly

library(tidyverse)
library(ggpubr)
library(tximport)
library(data.table)
library(cowplot)

source("~/lab_davidson/yan.a/software/scripts_denovo/R/quantile_overlap.R")

sqanti_summary <- readRDS('plot/sqanti_summary.rds')
# names(sqanti_summary) <- c('bambu', 'isonform','rattle','rnabloom2')
requested_counts <- readRDS('plot/requested_counts.rds')

args <- c('../bambu/',
          '../rattle/', 
          '../rnabloom2/',
          '../isonform/'
)

salmon_files <- lapply(args, function(x) {
  list.files(x, 'quant.sf', recursive = T, full.names = T)
}) %>% unlist()

salmon_files <- salmon_files[str_detect(salmon_files, 'dge') & !str_detect(salmon_files, 'bk')]

df <- data.frame(path = salmon_files) %>%
  separate(path, c(NA, 'assembler',NA, NA, "sample",NA), sep = '/', remove = F) %>%
  mutate(assembler = str_remove(assembler, 'stranded')) %>%
  # separate(dep1, c(NA, "depth", NA), sep = '_') %>%
  mutate(quant = ifelse(str_detect(sample, 'map|onts'), 'sec', 'pri')#,
         #depth = ifelse(depth == 'ilu', '10m', depth)
  ) %>%
  unite("sqanti_summary", assembler, #depth, 
        sep = '_', remove = F) 

df_list <- split(df, list( df$assembler, #df$depth, 
                           df$quant), drop = T)
# df_list <- Filter(function(x) nrow(x) > 0, df_list)

all_quant2 <- lapply(df_list, function(df){
  
  quant_summary <- lapply(1:3, function(x){
    
    counts <- tximport(df[x,'path'], 
                       type = "salmon", 
                       txOut = T)
    
    counts <- data.frame(counts = counts$counts[,1], 
                         isoform = rownames(counts$counts))
    
    quant <- quantile_overlap(sqanti_summary[[df[x,'sqanti_summary']]], 
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

saveRDS(all_quant2, 'plot/all_quant2.rds')

# these will be false as each data is filtered by there existing assembly
length(unique(lapply(all_quant2, function(x){
  x$tx_exp_tiled$tx
}))) == 1
length(unique(lapply(all_quant2, function(x){
  x$gene_exp_tiled$gene
}))) == 1

## this will be true
all(sapply(1:3, function(i) {
  identical(all_quant2[[i]]$tx_exp_tiled$tx, 
            all_quant2[[i + 3]]$tx_exp_tiled$tx)
}))
all(sapply(1:3, function(i) {
  identical(all_quant2[[i]]$gene_exp_tiled$gene, 
            all_quant2[[i + 3]]$gene_exp_tiled$gene)
}))

# load bambu reference
tx_counts_bambu <- read.table('../bambu/PB_merged_dge/default/counts_transcript.txt', 
                              header = T, sep = '\t') 
rownames(tx_counts_bambu) <- tx_counts_bambu$TXNAME

gene_counts_bambu <- tx_counts_bambu[, -1] %>%
  filter(!str_detect(GENEID, '^Bambu')) %>%
  group_by(GENEID) %>%
  summarise_all(sum) %>%
  data.frame()
rownames(gene_counts_bambu) <- gene_counts_bambu$GENEID
gene_counts_bambu <- gene_counts_bambu[,2:4]
dim(gene_counts_bambu)

tx_counts_bambu <- tx_counts_bambu[,3:5]
tx_counts_bambu <- tx_counts_bambu[rowSums(tx_counts_bambu) != 0, ]
dim(tx_counts_bambu)

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

saveRDS(cor_noseq2, 'plot/cor_noseq2.rds')


