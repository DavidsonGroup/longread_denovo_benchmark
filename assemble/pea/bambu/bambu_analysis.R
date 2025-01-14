#!/usr/bin/env Rscript

library(bambu)
library(tidyverse)
# library(BiocFileCache)

args = commandArgs(trailingOnly=TRUE)

bam.files <- list.files(args[1], pattern = 'bam$', full.names = TRUE)
bam.files <- normalizePath(bam.files)
fa.file <- "/home/users/allstaff/yan.a/lab_davidson/yan.a/ref/ncbi/pea_zw6/GCF_024323335.1_CAAS_Psat_ZW6_1.0_genomic.fna"
gtf.file <- "/home/users/allstaff/yan.a/lab_davidson/yan.a/ref/ncbi/pea_zw6/GCF_024323335.1_CAAS_Psat_ZW6_1.0_genomic.gtf"

bambuAnnotations <- prepareAnnotations(gtf.file)

## run default
dir.create(paste0(args[1],'/default'))

se1 <- bambu(reads = bam.files, 
             annotations = bambuAnnotations, 
             stranded = T, 
             genome = fa.file, ncore = 32)

data.frame(mcols(se1)) %>% count(novelTranscript, is.na(readCount))

writeBambuOutput(se1, path = paste0(args[1],'/default'))

saveRDS(se1, paste0(args[1],'/default/se1.rds'))



