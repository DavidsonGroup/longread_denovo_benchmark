#!/usr/bin/env Rscript

## prepare seurat object for gene and tx
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Matrix))
library(Matrix)

source('/vast/scratch/users/yan.a/vast_scratch/pb/1m_stranded/bambu/collapse_transcripts_to_genes.R')
source('~/lab_davidson/yan.a/software/scripts_denovo/R/get_gene_tx.R')

args <- commandArgs(trailingOnly=TRUE)

# dir <- 'PB_merged_oarfish/'

# args[1] <- '../rnabloom2/'
# args[2] <- "../rnabloom2/PB_merged_corset/PB_merged-clusters_mod.txt"
# args[3] <- 'rnabloom2_corset'
# args[4] <- 'oarfish'

# args <- NULL
# args[1] <- '../rattle/'
# args[2] <- "../rattle/fx2tab.txt"
# args[3] <- 'rattle_rattle'
# args[4] <- 'oarfish'

outputdir <- paste(args[3], args[4], sep = '_')

dir.create(outputdir)

method <- str_split(args[3], pattern = '_', simplify = T)[2]

dirs <- list.dirs(args[1], recursive = F, full.names = T)

dir_oarfish <- dirs[grepl('oarfish$', dirs)]

exp1 <- ReadMtx(
  mtx = paste0(dir_oarfish, "/PBMC_10k_oarfish_quant.count.mtx"), 
  features = paste0(dir_oarfish, "/PBMC_10k_oarfish_quant.features.txt"),
  cells = paste0(dir_oarfish, "/PBMC_10k_oarfish_quant.barcodes.txt"),
  feature.column = 1,
  mtx.transpose = T
)
rownames(exp1) <- str_remove(rownames(exp1), "\\|.*")

# for gencode tx ids
# tx2gene <- data.frame(fullname = rownames(exp1)) %>%
#   tidyr::separate(
#     fullname,
#     into = c("txid","geneid",NA,NA,"txname","genename","length","type"),
#     sep = "\\|",
#     fill = "right",
#     remove = FALSE
#   )

# gencode.tx <- read.table('oarfish/PBMC_10k_oarfish_quant.features.txt', sep = '\t', header = F) %>%
#   tidyr::separate(
#     V1,
#     into = c("txid","geneid",NA,NA,"txname","genename","length","type"),
#     sep = "\\|",
#     fill = "right",
#     remove = FALSE
#   ) 
# gencode.gene <- gencode.tx %>% select(genename, geneid) %>% distinct()

rt_gene <- get_gene_tx(args[2],
                       list.files(dirs[grepl('sqanti3$', dirs)], pattern = 'classification.txt', full.names = T, recursive = T),
                       method = method)
rownames(rt_gene) <- rt_gene$isoform

saveRDS(rt_gene, paste0(outputdir,'/rt_gene.rds'))

# tx_table <- args[2]
# 
# if (method == 'rattle') {
#   
#   ## use fx2tab results
#   ## command: seqkit fx2tab -nlgGH transcriptome.fasta > fx2tab.txt
#   
#   table <- read.table(tx_table, header = F, sep = '\t') %>% 
#     separate(V1, c('isoform','geneid','origin_txid','counts','label'), sep = ' ') %>% 
#     dplyr::mutate(counts = str_remove(counts, 'total_reads=') %>% as.numeric()) 
#   
# } else if (method %in% c('corset','rnabloom2','bambu','gencode')) {
#   
#   ## result after run corset on rnabloom2 assembly
#   
#   table <- read.table(tx_table, 
#                       sep = '\t', header = F) %>%
#     dplyr::rename(isoform = V1, geneid = V2)
#   
# } else if (method %in% c('trinity','rnaspades')) {
#   
#   table <- read.table(tx_table, 
#                       sep = '\t', header = F) %>%
#     dplyr::rename(isoform = V2, geneid = V1)
#   
# } else if (method == 'isonclust') {
#   
#   ## use fx2tab results
#   ## command: seqkit fx2tab -nlgGH transcriptome.fasta > fx2tab.txt
#   
#   table <- read.table(tx_table, header = F, sep = '\t') %>% 
#     separate(V1, c('geneid','batchid','id'), sep = '_', remove = F) %>%  
#     dplyr::rename(isoform = V1)
# }
# rownames(table) <- table$isoform
# 
# saveRDS(table, paste0(outputdir,'/tx2gene.rds'))

# tx2gene <- read.table('PB_merged_corset/', sep = '\t', header = F)
# colnames(tx2gene) <- c('txid','geneid')
# 
# tx2gene2 <- left_join(tx2gene, gencode.tx, by = c('txid', 'geneid')) %>%
#   mutate(txname = ifelse(is.na(txname), txid, txname)) %>%
#   left_join(gencode.gene, by = c('geneid')) %>%
#   mutate(genename = ifelse(is.na(genename.y), geneid, genename.y))
# rownames(tx2gene2) <- tx2gene2$txid

gene_counts1 <- collapse_transcripts_to_genes(exp1, rt_gene[rownames(exp1), ], gene_col = "geneid")
dim(gene_counts1)

exp2 <- ReadMtx(
  mtx = paste0(dir_oarfish,"/PBMC_20k_rep1_oarfish_quant.count.mtx"), 
  features = paste0(dir_oarfish,"/PBMC_20k_rep1_oarfish_quant.features.txt"),
  cells = paste0(dir_oarfish,"/PBMC_20k_rep1_oarfish_quant.barcodes.txt"),
  feature.column = 1,
  mtx.transpose = T
)
rownames(exp2) <- str_remove(rownames(exp2), "\\|.*")
gene_counts2 <- collapse_transcripts_to_genes(exp2, rt_gene[rownames(exp2), ], gene_col = "geneid")

exp3 <- ReadMtx(
  mtx = paste0(dir_oarfish, "/PBMC_20k_rep2_oarfish_quant.count.mtx"), 
  features = paste0(dir_oarfish,"/PBMC_20k_rep2_oarfish_quant.features.txt"),
  cells = paste0(dir_oarfish,"/PBMC_20k_rep2_oarfish_quant.barcodes.txt"),
  feature.column = 1,
  mtx.transpose = T
)
rownames(exp3) <- str_remove(rownames(exp3), "\\|.*")
gene_counts3 <- collapse_transcripts_to_genes(exp3, rt_gene[rownames(exp3), ], gene_col = "geneid")

stopifnot(all.equal(sort(unique(rt_gene$geneid)), 
                    sort(row.names(gene_counts1)),
                    sort(row.names(gene_counts2)),
                    sort(row.names(gene_counts3))))

stopifnot(length(unique(sapply(list(exp1, exp2, exp3), nrow))) == 1)
stopifnot(length(unique(sapply(list(gene_counts1, gene_counts2, gene_counts3), nrow))) == 1)

## gene
# note seurat will convert _ to -, so need to change it back later
se1 <- CreateSeuratObject(counts = gene_counts1)
se2 <- CreateSeuratObject(counts = gene_counts2)
se3 <- CreateSeuratObject(counts = gene_counts3)

merged_gene <- merge(
  x = se1,
  y = list(se2, se3),
  add.cell.ids = c("rep1_10k","rep1_20k","rep2_20k")  # prefixes for unique cell IDs
)

merged_gene
# An object of class Seurat 
# 205239 features across 50968 samples within 1 assay 
# Active assay: RNA (205239 features, 0 variable features)
# 3 layers present: counts.1, counts.2, counts.3

saveRDS(merged_gene, paste0(outputdir,'/merged_gene.rds'))

gene_counts <- list(data.frame(gene = merged_gene[['RNA']]$counts.1 %>% rownames(), 
                               rep1_10k_counts = merged_gene[['RNA']]$counts.1 %>% rowSums()),
                    data.frame(gene = merged_gene[['RNA']]$counts.2 %>% rownames(), 
                               rep1_20k_counts = merged_gene[['RNA']]$counts.2 %>% rowSums()),
                    data.frame(gene = merged_gene[['RNA']]$counts.3 %>% rownames(), 
                               rep2_20k_counts = merged_gene[['RNA']]$counts.3 %>% rowSums())) %>% 
  reduce(left_join, by = 'gene') %>%
  mutate(gene = str_replace_all(gene, '-', '_'))

write.csv(gene_counts, paste0(outputdir,'/gene_counts_join.csv'))

## tx
# note seurat will convert _ to -, so need to change it back later
se1 <- CreateSeuratObject(counts = exp1)
se2 <- CreateSeuratObject(counts = exp2)
se3 <- CreateSeuratObject(counts = exp3)

merged_tx <- merge(
  x = se1,
  y = list(se2, se3),
  add.cell.ids = c("rep1_10k","rep1_20k","rep2_20k")  # prefixes for unique cell IDs
)

merged_tx
# An object of class Seurat 
# 508641 features across 50968 samples within 1 assay 
# Active assay: RNA (83687 features, 0 variable features)
# 3 layers present: counts.1, counts.2, counts.3

saveRDS(merged_tx, paste0(outputdir,'/merged_tx.rds'))

tx_counts <- list(data.frame(tx = merged_tx[['RNA']]$counts.1 %>% rownames(), 
                             rep1_10k_counts = merged_tx[['RNA']]$counts.1 %>% rowSums()),
                  data.frame(tx = merged_tx[['RNA']]$counts.2 %>% rownames(), 
                             rep1_20k_counts = merged_tx[['RNA']]$counts.2 %>% rowSums()),
                  data.frame(tx = merged_tx[['RNA']]$counts.3 %>% rownames(), 
                             rep2_20k_counts = merged_tx[['RNA']]$counts.3 %>% rowSums())) %>% 
  reduce(left_join, by = 'tx') %>%
  mutate(gene = str_replace_all(tx, '-', '_'))

write.csv(tx_counts, paste0(outputdir,'/tx_counts_join.csv'))

# sanity check if barcodes in subset of data is present in full data
merged_se <- readRDS('~/lab_davidson/yan.a/pb/seurat/merged_se.rds')
sum(colnames(merged_gene) %in% str_remove(colnames(merged_se), '-1$'))
sum(colnames(merged_tx) %in% str_remove(colnames(merged_se), '-1$'))

cat('Finished generating tx and gene count tables!\n')
