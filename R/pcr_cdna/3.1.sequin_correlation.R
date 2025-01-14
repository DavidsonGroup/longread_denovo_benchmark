
library(tidyverse)
library(ggpubr)
library(tximport)
library(cowplot)

source("~/lab_davidson/yan.a/software/scripts_denovo/R/quantile_overlap.R")
source('~/lab_davidson/yan.a/software/scripts_denovo/R/calc_sequin_cor.R')

sqanti_summary <- readRDS('plot/sqanti_summary.rds')
requested_counts <- readRDS('plot/requested_counts.rds')

# load sequin de 
anno <- read.table("~/lab_davidson/yan.a/software/LongReadRNA/sequins/annotations/rnasequin_genes_2.4.tsv",
                   header = TRUE, stringsAsFactors = FALSE) %>%
  mutate(logFC_expected = log(MIX_A / MIX_B),
         isChange = logFC_expected != 0
  )

anno2 <- read.table("~/lab_davidson/yan.a/software/LongReadRNA/sequins/annotations/rnasequin_isoforms_2.4.tsv", 
                    header = TRUE, stringsAsFactors = FALSE) %>%
  mutate(logFC_expected = log(MIX_A / MIX_B),
         isChange = logFC_expected != 0
  )

# for sequin we need to use the individual quantification for de novo tools
filename_list <- Sys.glob("../*/*dge/*/quant.sf") 

df <- filename_list %>%
  str_split(pattern = '/|_quant', simplify = T) %>%
  data.frame() %>%
  mutate(filepath = filename_list) %>%
  filter(!(X2 %in% c('bambu_bk','limma','trinity'))) %>%
  separate(X4, c('sample','depth'),extra = 'drop') %>%
  separate(X5, c(NA, c('assembler','quant'))) %>%
  filter(X2 %in% c('isonform', 'rattle','rnabloom2','trinitystranded')) %>%
  mutate(assembler = factor(X2, levels = c('rattle','rnabloom2','isonform', 'trinitystranded'),
                         labels = c('rattle','rnabloom2','isonform', 'trinity')),
         mix = ifelse(sample %in% c(paste0('barcode0', 1:3), 'H1975'), 'mix_a', 
                        ifelse(sample %in% c(paste0('barcode0', 4:6), 'HCC827'), 'mix_b',NA)),
         depth = ifelse(assembler != 'trinity', depth, '10m')) %>%
  select(assembler, quant, depth, mix, filepath)

result <- split(df$filepath, 
                list(df$quant, df$mix, df$depth, df$assembler))
result <- result[sapply(result, length) > 0]

# for de novo
sqanti_summary_tmp <- sqanti_summary[c(rep(6:8, each=4), # rattle 
                                       rep(9:11, each=4), # rnablomm2
                                       rep(4:5, each = 4), # isonform
                                       rep(12,4) # trinity
                                       )]

sequin_cor <- lapply(1:length(result), function(x){
  # assembler <- str_split(names(result)[x], '\\.', simplify = T)[1]
  
  calc_sequin_cor(sqanti = sqanti_summary_tmp[[x]], files = result[[x]], isoform_count = NULL, 
                  requested_counts = requested_counts)
})

names(sequin_cor) <- names(result)


# for bambu quantification

bambu_counts <- list.files('../bambu/', 'counts_transcript.txt', full.names = T, recursive = T)[1:3]
bambu_isoform_count_list <- append(lapply(bambu_counts, function(x){
  df <- read.table(x, header = T, sep = '\t')
  res <- data.frame(isoform = df$TXNAME, counts = rowSums(df[,3:5]))
}),
lapply(bambu_counts, function(x){
  df <- read.table(x, header = T, sep = '\t')
  res <- data.frame(isoform = df$TXNAME, counts = rowSums(df[,6:8]))
}))
names(bambu_isoform_count_list) <- c('count.mix_a.10m.bambu', 'count.mix_a.2m.bambu','count.mix_a.5m.bambu',
                                     'count.mix_b.10m.bambu', 'count.mix_b.2m.bambu','count.mix_b.5m.bambu')

sqanti_summary_tmp <- sqanti_summary[rep(1:3, 2)]

sequin_cor2 <- lapply(1:length(bambu_isoform_count_list), function(x){
  # assembler <- str_split(names(result)[x], '\\.', simplify = T)[1]
  
  calc_sequin_cor(sqanti = sqanti_summary_tmp[[x]], isoform_count = bambu_isoform_count_list[[x]], 
                  requested_counts = requested_counts)
})

names(sequin_cor2) <- names(bambu_isoform_count_list)

saveRDS(append(sequin_cor2, sequin_cor), 'plot/sequin_cor.rds')




