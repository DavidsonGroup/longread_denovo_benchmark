#!/usr/bin/env Rscript

library(bambu)
library(tidyverse)
# library(BiocFileCache)

args = commandArgs(trailingOnly=TRUE)
# args="ONT_merged_dge/"
# args="../bambu/ONT_merged_dge/"

bam.files <- list.files(args[1], pattern = 'bam$', full.names = TRUE)
bam.files <- normalizePath(bam.files)
fa.file <- "/vast/projects/lab_davidson/yan.a/Dong_2023/reference/genome_GRCh38_sequin.fa"
# gtf.file <- "../../sqanti_sim/gencode.v44.subset.gtf"

# bambuAnnotations <- prepareAnnotations(gtf.file)

## run default
dir.create(paste0(args[1],'/default'))

se1 <- bambu(reads = bam.files, 
            # annotations = bambuAnnotations, 
             stranded = F, 
             genome = fa.file, ncore = 32,
            NDR = 1
            )

# se2 <- bambu(reads = bam.files, 
#              stranded = F, 
#              genome = fa.file, ncore = 32)

# NDR=0.5, bambu 3.11
# 23074/11704

# NDR=1, bambu 3.11
# did not finish
# 312564/19628

# NDR=defult, bambu 3.11
# 23074/11704

# NDR=0.1, bambu 3.11
# 3333/3153

# NDR=1, bambu 3.6/3.3.4
# 50103/13772


data.frame(mcols(se1)) %>% count(novelTranscript, is.na(readCount))

writeBambuOutput(se1, path = paste0(args[1],'/default'))

saveRDS(se1, paste0(args[1],'/default/se1.rds'))



