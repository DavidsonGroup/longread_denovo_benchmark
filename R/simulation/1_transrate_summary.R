
library(tidyverse)
library(data.table)
source('~/lab_davidson/yan.a/software/scripts_denovo/R/get_stats.R')
source('~/lab_davidson/yan.a/software/scripts_denovo/R/plot_transrate.R')

# set plotting parameters
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cols <- cbPalette[c(1,2,6,4,7,8)]
names(cols) <- c('limma','bambu','isonform','rattle','rnabloom2','trinity')
shapes <- c(15,15,16,17,17,17)
names(shapes) <- c('sim','bambu','corset','isonclust','rattle','trinity')

args <- commandArgs(trailingOnly=TRUE)

args <- c('../simulation_1m/bambu/',
          '../simulation_1m/rattle/', 
          '../simulation_1m/rnabloom2/',
          '../simulation_1m/trinity/',
          '../simulation_1m/isonform/',
          '../sqanti_sim/'
          #'../bambu/',
          #'../isonform/'
)

dirs <- list.dirs(args, recursive = F, full.names = T)
test <- dirs[grepl('transrate$', dirs)]

testlist <- lapply(test, function(x){
  
  get_stats(x)$seqstats %>% # select(1,2,3,6) %>%
    mutate(assembler = dirname(x) %>% str_remove('../') %>% str_remove('simulation_1m/')) 
  
}) %>%
  rbindlist() %>%
  mutate(depth = 1,
         assembler = factor(assembler, 
                            levels = c('sqanti_sim','bambu','rattle','rnabloom2','isonform','trinity'),
                            labels = c('ref','bambu','rattle','rnabloom2','isonform','trinity')))
write.csv(testlist, 'plot/transrate_metric.csv')

plot_transrate(testlist, depth = F)

# read.stats <- fread('../sqanti_sim/1m/ONT_merged.seqlength', header = F)
# colnames(read.stats) <- c('contig_name', 'length')
# read.stats$method <- 'raw_reads'
# 
# rbindlist(list(testlist, read.stats), fill = T) %>% 
#   ggplot(aes(x=method, y=length)) +
#   geom_boxplot(aes(fill = method)) +
#   theme(axis.text.x=element_text(angle=45, hjust=1)) +
#   stat_summary(fun.y = function(x) log10(mean(10^x)), geom = "point") +  # Add points to plot ,  col = "red"
#   stat_summary(fun.y = function(x) log10(mean(10^x)), geom = "text", # col = "red",     # Add text to plot
#                vjust = 1.5, aes(label = round(..y.., digits = 0))) +
#   scale_y_log10() +
#   ylab('Transcript length (bp)')
# ggsave('plot/transcript_length_with_raw.pdf', width = 5, height = 4)


