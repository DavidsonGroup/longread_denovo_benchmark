#!/usr/bin/env Rscript

library(bambu) # 3.3.4
library(tidyverse)
# library(BiocFileCache)

args = commandArgs(trailingOnly=TRUE)

bam.files <- list.files(args[1], pattern = 'bam$', full.names = TRUE)
bam.files <- normalizePath(bam.files)
fa.file <- '/vast/projects/lab_davidson/yan.a/ref/gencode/GRCh38.primary_assembly.genome.fa'
gtf.file <- '/vast/projects/lab_davidson/yan.a/ref/gencode/gencode.v44.annotation.gtf'

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



