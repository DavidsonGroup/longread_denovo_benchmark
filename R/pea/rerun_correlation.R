library(tidyverse)
library(ggpubr)
library(tximport)
library(data.table)
library(cowplot)

source("~/lab_davidson/yan.a/software/scripts_denovo/R/quantile_overlap.R")

sqanti_summary <- readRDS('plot/sqanti_summary.rds')
requested_counts <- readRDS('plot/requested_counts.rds')

args <- c('../bambu/',
          '../rattle/', 
          '../rnabloom2/',
          '../isonform/',
          '../trinitystranded/'
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

all_quant <- lapply(df_list, function(df){
  
  quant_summary <- lapply(1:12, function(x){
    
    counts <- tximport(df[x,'path'], 
                       type = "salmon", 
                       txOut = T)
    
    counts <- data.frame(counts = counts$counts[,1], 
                         isoform = rownames(counts$counts))
    
    quant <- quantile_overlap(sqanti_summary[[df[x,'sqanti_summary']]], 
                              counts, 'isoform', requested_counts)
    
  })
  
  tx_exp_tiled <- lapply(quant_summary, function(x) {
    data.frame(tx = rownames(x$lcpm_tx),
               assemble_cpm.max = x$lcpm_tx$assemble_cpm.max, # take max if multiple de novo tx match to 1 ref tx
               counts.max = x$lcpm_tx$counts.max)
  }) %>% rbindlist()
  
  gene_exp_tiled <- lapply(quant_summary, function(x) {
    data.frame(gene = rownames(x$lcpm_gene),
               assemble_cpm = x$lcpm_gene$assemble_cpm, # sum up all tx with same sqanti3 associated gene 
               assemble_count = x$lcpm_gene$counts)
  }) %>% rbindlist()
  
  res <- list(tx_exp_tiled = tx_exp_tiled,
              gene_exp_tiled = gene_exp_tiled)
})

saveRDS(all_quant, 'plot/all_quant.rds')

# make sure all list have same order
length(unique(lapply(all_quant, function(x){
  x$tx_exp_tiled$tx
}))) == 1
length(unique(lapply(all_quant, function(x){
  x$gene_exp_tiled$gene
}))) == 1

# load bambu reference
tx_counts_bambu <- readRDS('../bambu/pea_dge/R/y.rds')
tx_counts_bambu <- tx_counts_bambu[rowSums(tx_counts_bambu$counts) != 0, ]
dim(tx_counts_bambu)

gene_counts_bambu <- data.frame(tx_counts_bambu$counts, 
                                gene_id = tx_counts_bambu$genes$gene_id) %>%
  filter(!str_detect(gene_id, '^Bambu')) %>%
  group_by(gene_id) %>%
  summarise_all(sum) %>% 
  data.frame()
rownames(gene_counts_bambu) <- gene_counts_bambu$gene_id
# gene_counts_bambu <- readRDS('../sqanti_sim/design/limma_gene_obj.rds')
# gene_counts_bambu <- gene_counts_bambu[rowSums(gene_counts_bambu$counts) != 0, ]
dim(gene_counts_bambu)

# for non sequin
geneorder <- all_quant[[1]]$gene_exp_tiled %>%
  filter(!str_detect(gene, '^R')) %>% pull(gene) %>% unique()

truecpm_gene <- log2(gene_counts_bambu[,-1] + 1)[geneorder, ] %>% 
  as.matrix() %>%
  c()

txorder <- all_quant[[1]]$tx_exp_tiled %>%
  filter(!str_detect(tx, '^R')) %>% pull(tx) %>% unique()

truecpm_tx <- log2(tx_counts_bambu$counts + 1)[txorder, ] %>%
  as.matrix() %>%
  c()

# 
# cor(all_quant$rattle.10m.pri$tx_exp_tiled$assemble_cpm, truecpm_tx)

# calculate correlation with or withou sequin
# gene_sequin <- grep('^R', all_quant[[1]]$gene_exp_tiled$gene)
# tx_sequin <- grep('^R', all_quant[[1]]$tx_exp_tiled$tx)
  
cor_noseq <- lapply(all_quant[-c(3,6)], function(x){
  cor.tx <- cor(x$tx_exp_tiled$assemble_cpm, truecpm_tx)
  cor.gene <- cor(x$gene_exp_tiled$assemble_cpm, truecpm_gene)
  res <- data.frame(cor.tx = cor.tx, cor.gene = cor.gene)
}) %>% rbindlist() %>%
  mutate(data = names(all_quant)[-c(3,6)])

# add trinity quant (order is different) 
# in trinity rep(c('Drog','Dwt','Trog','Twt'), each = 3)}, bambu rep(c('Twt','Trog','Dwt','Drog'), each = 3)
tx_counts_bambu <- readRDS('../bambu/pea_dge/R/y.rds')
tx_counts_bambu <- tx_counts_bambu[rowSums(tx_counts_bambu$counts) != 0, c(10:12, 7:9, 4:6, 1:3)]
dim(tx_counts_bambu)

gene_counts_bambu <- data.frame(tx_counts_bambu$counts, 
                                gene_id = tx_counts_bambu$genes$gene_id) %>%
  filter(!str_detect(gene_id, '^Bambu')) %>%
  group_by(gene_id) %>%
  summarise_all(sum) %>% 
  data.frame()
rownames(gene_counts_bambu) <- gene_counts_bambu$gene_id
# gene_counts_bambu <- readRDS('../sqanti_sim/design/limma_gene_obj.rds')
# gene_counts_bambu <- gene_counts_bambu[rowSums(gene_counts_bambu$counts) != 0, ]
dim(gene_counts_bambu)

# for non sequin
geneorder <- all_quant[[1]]$gene_exp_tiled %>%
  filter(!str_detect(gene, '^R')) %>% pull(gene) %>% unique()

truecpm_gene <- log2(gene_counts_bambu[,-1] + 1)[geneorder, ] %>% 
  as.matrix() %>%
  c()

txorder <- all_quant[[1]]$tx_exp_tiled %>%
  filter(!str_detect(tx, '^R')) %>% pull(tx) %>% unique()

truecpm_tx <- log2(tx_counts_bambu$counts + 1)[txorder, ] %>%
  as.matrix() %>%
  c()

cor_noseq_trinity <- lapply(all_quant[c(3,6)], function(x){
  cor.tx <- cor(x$tx_exp_tiled$assemble_cpm, truecpm_tx)
  cor.gene <- cor(x$gene_exp_tiled$assemble_cpm, truecpm_gene)
  res <- data.frame(cor.tx = cor.tx, cor.gene = cor.gene)
}) %>% rbindlist() %>%
  mutate(data = names(all_quant)[c(3,6)])

cor_noseq <- rbind(cor_noseq, cor_noseq_trinity)

saveRDS(cor_noseq, 'plot/cor_noseq.rds')


