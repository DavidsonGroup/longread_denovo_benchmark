
library(tidyverse)
library(data.table)
library(GenomicFeatures)
source('~/lab_davidson/yan.a/software/scripts_denovo/R/get_stats.R')
source('~/lab_davidson/yan.a/software/scripts_denovo/R/plot_transrate.R')

# set plotting parameters
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#c51b8a")
cols <- cbPalette[c(1,2,6,4,7,8, 3,5,9)]
names(cols) <- c('limma','bambu','isonform','rattle','rnabloom2','trinity','bambudenovo','rnaspades','rnabloom2hybrid')
shapes <- c(15,15,16,17,17,17,17,17,17)
names(shapes) <- c('sim','bambu','corset','isonclust','rattle','trinity','bambudenovo','rnaspades','rnabloom2hybrid')

## set folders to work
args <- commandArgs(trailingOnly=TRUE)

args <- c('../bambu/',
          '../bambudenovo/',
          '../rattle/', 
          '../rnabloom2/',
          '../isonform/',
          '../trinitystranded/',
          '../rnaspades/',
          '../rnabloom2hybrid/'
)

# get files
dirs <- list.dirs(args, recursive = F, full.names = T)
test <- dirs[grepl('transrate$', dirs)]

# prepare dataframe for plotting
# if only 1 depth, force all depth to be 1
testlist <- lapply(test, function(x){
  
  get_stats(x)$seqstats %>% # dplyr::select(1,2,3,6) %>%
    mutate(assembler = dirname(x) %>% str_remove('../'),
           depth = basename(x) %>% str_remove('_transrate$') %>% str_remove('^merged_')) 
  
}) %>%
  rbindlist() %>%
  mutate(depth = ifelse(depth %in% c('ilu', 'hybrid_merged'), '10m', depth)) %>% 
  mutate(depth = factor(depth, levels = c('2m','5m','10m'), labels = c(2,5,10)),
         assembler = factor(assembler, 
                            levels = c('bambu','rattle','rnabloom2','isonform','trinitystranded','bambudenovo','rnaspades','rnabloom2hybrid'),
                            labels = c('bambu','rattle','rnabloom2','isonform','trinity','bambudenovo','rnaspades','rnabloom2hybrid')))
write.csv(testlist, 'plot/transrate_metric.csv')

# all different plots
plot_transrate(testlist, depth = TRUE)

# prepare data for plotting
# read.stats <- fread('../sqanti_sim/1m/ONT_merged.seqlength', header = F)
# colnames(read.stats) <- c('contig_name', 'length')
# read.stats$assembler <- 'raw_reads'
#
# bambu.stats <- makeTxDbFromGFF('../bambu/merged_10m_dge/default/extended_annotations.gtf') %>%
#   transcriptLengths() %>%
#   dplyr::select(tx_name, tx_len) %>%
#   dplyr::rename(contig_name = tx_name, length = tx_len)
# 
# rbindlist(list(testlist, 
#                bambu.stats %>% mutate(depth = '10', assembler = 'bambu')), 
#           fill = T) %>% 
#   ggplot(aes(x=assembler, y=length, fill = assembler, alpha = depth)) +
#   geom_boxplot(position = position_dodge2(preserve = "single")) +
#   theme(axis.text.x=element_text(angle=45, hjust=1)) +
#   stat_summary(fun.y = function(x) log10(mean(10^x)), geom = "point",
#                position = position_dodge(width = 0.75)) +  # Add points to plot ,  col = "red"
#   stat_summary(fun.y = function(x) log10(mean(10^x)), geom = "text", # col = "red",     # Add text to plot
#                vjust = 1.5, aes(label = round(..y.., digits = 0)),
#                position = position_dodge(width = 0.75)) +
#   scale_y_log10() +
#   scale_fill_manual(values = cols) + 
#   ylab('Transcript length (bp)')
# ggsave('plot/transcript_length_with_bambu.pdf', width = 8, height = 4)
