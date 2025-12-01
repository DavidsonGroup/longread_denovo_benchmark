
library(tidyverse)
library(data.table)
library(GenomicFeatures)
library(scales)
source('~/lab_davidson/yan.a/software/scripts_denovo/R/get_stats.R')
source('~/lab_davidson/yan.a/software/scripts_denovo/R/plot_transrate.R')

# set plotting parameters
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cols <- cbPalette[c(1,2,6,4,7,8)]
names(cols) <- c('ref','bambu','isonform','rattle','rnabloom2','trinity')
shapes <- c(15,15,16,17,17,17)
names(shapes) <- c('ref','bambu','corset','isonclust','rattle','trinity')

args <- commandArgs(trailingOnly=TRUE)

args <- c('../bambu/',
          '../rattle/', 
          '../rnabloom2/',
          '../trinitystranded/',
          '../isonform/'
)

dirs <- list.dirs(args, recursive = F, full.names = T)
test <- dirs[grepl('transrate$', dirs)]

# bambu.stats <- makeTxDbFromGFF('../bambu/A549_MCF7_directRNA_dge/default/extended_annotations.gtf') %>%
#   transcriptLengths() %>%
#   dplyr::select(tx_name, tx_len) %>%
#   dplyr::rename(contig_name = tx_name, length = tx_len)

testlist <- lapply(test, function(x){
  
  get_stats(x)$seqstats %>% #dplyr::select(1,2,3,6) %>%
    mutate(assembler = dirname(x) %>% str_remove('../') )  
  
}) %>%
  rbindlist() %>% 
  mutate(depth = 10, # placehold for plotting
         assembler = factor(assembler, 
                            levels = c('bambu','rattle','rnabloom2','isonform','trinitystranded'),
                            labels = c('bambu','rattle','rnabloom2','isonform','trinity'))
         )
write.csv(testlist, 'plot/transrate_metric.csv')

plot_transrate(testlist, depth = F)


