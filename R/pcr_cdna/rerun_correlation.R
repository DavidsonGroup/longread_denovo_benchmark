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
          '../rattle/', 
          '../rnabloom2/',
          '../isonform/',
          '../trinitystranded/'
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
         depth = ifelse(depth == 'ilu', '10m', depth)) %>%
  unite("sqanti_summary", assembler, depth, sep = '_', remove = F) 

df_list <- split(df, list( df$assembler, df$depth, df$quant), drop = T)
# df_list <- Filter(function(x) nrow(x) > 0, df_list)

all_quant <- lapply(df_list, function(df){
  
  quant_summary <- lapply(1:6, function(x){
    
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
geneorder <- all_quant$rattle.10m.pri$gene_exp_tiled %>%
  filter(!str_detect(gene, '^R')) %>% pull(gene) %>% unique()

truecpm_gene <- log2(gene_counts_bambu[,-1] + 1)[geneorder, ] %>% 
  as.matrix() %>%
  c()

txorder <- all_quant$rattle.10m.pri$tx_exp_tiled %>%
  filter(!str_detect(tx, '^R')) %>% pull(tx) %>% unique()

truecpm_tx <- log2(tx_counts_bambu$counts + 1)[txorder, ] %>%
  as.matrix() %>%
  c()

# 
# cor(all_quant$rattle.10m.pri$tx_exp_tiled$assemble_cpm, truecpm_tx)

# calculate correlation with or withou sequin
gene_sequin <- grep('^R', all_quant$rattle.10m.pri$gene_exp_tiled$gene)
tx_sequin <- grep('^R', all_quant$rattle.10m.pri$tx_exp_tiled$tx)
  
cor_noseq <- lapply(all_quant, function(x){
  cor.tx <- cor(x$tx_exp_tiled$assemble_cpm[-tx_sequin], truecpm_tx)
  cor.gene <- cor(x$gene_exp_tiled$assemble_cpm[-gene_sequin], truecpm_gene)
  res <- data.frame(cor.tx = cor.tx, cor.gene = cor.gene)
}) %>% rbindlist() %>%
  mutate(data = names(all_quant))

## sequin truth 

# load sequin de 

sequin.gene.order <- all_quant$rattle.10m.pri$gene_exp_tiled %>% 
  filter(str_detect(gene, '^R')) %>% pull(gene) %>% unique()
sequin.gene.exp <- anno.gene[, c('a','a','a','b','b','b')] 
rownames(sequin.gene.exp) <- anno.gene$NAME
sequin.gene.exp <- sequin.gene.exp[sequin.gene.order, ] %>%
  as.matrix() %>%
  c()

sequin.tx.order <- all_quant$rattle.10m.pri$tx_exp_tiled %>% 
  filter(str_detect(tx, '^R')) %>% pull(tx) %>% unique()
sequin.tx.exp <- anno.tx[, c('a','a','a','b','b','b')] 
rownames(sequin.tx.exp) <- anno.tx$NAME
sequin.tx.exp <- sequin.tx.exp[sequin.tx.order, ] %>%
  as.matrix() %>%
  c()

cor_seq <- lapply(all_quant, function(x){
  cor.tx <- cor(x$tx_exp_tiled$assemble_cpm[tx_sequin], sequin.tx.exp)
  cor.gene <- cor(x$gene_exp_tiled$assemble_cpm[gene_sequin], sequin.gene.exp)
  res <- data.frame(cor.tx = cor.tx, cor.gene = cor.gene)
}) %>% rbindlist() %>%
  mutate(data = names(all_quant))

bambu_dir <- list.dirs('../bambu/')
bambu_dir <- bambu_dir[str_detect(bambu_dir, 'dge/R$')]

cor_seq_bambu <- lapply(bambu_dir, function(x){
  gene <- readRDS(paste0(x, '/x.rds'))
  gene <- log2(gene$counts + 1)[sequin.gene.order, ] %>% 
    as.matrix() %>% c()
  cor.gene <- cor(gene, sequin.gene.exp)
  
  tx <- readRDS(paste0(x, '/y.rds'))
  tx <- log2(tx$counts + 1)[sequin.tx.order, ] %>% 
    as.matrix() %>% c()
  cor.tx <- cor(tx, sequin.tx.exp)  
  res <- data.frame(cor.tx = cor.tx, cor.gene = cor.gene)
}) %>%
  rbindlist() %>%
  mutate(data = c('bambu.10m.sec','bambu.2m.sec','bambu.5m.sec'))

cor_seq <- rbind(cor_seq_bambu, cor_seq)

saveRDS(cor_noseq, 'plot/cor_noseq.rds')
saveRDS(cor_seq, 'plot/cor_seq.rds')

