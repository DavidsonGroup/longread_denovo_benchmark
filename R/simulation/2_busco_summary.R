
library(tidyverse)
library(data.table)
library(ggpubr)
source('~/lab_davidson/yan.a/software/scripts_denovo/R/read_busco.R')

# set plotting parameters
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#c51b8a")
cols <- cbPalette[c(1,2,6,4,7,8, 3,5,9)]
names(cols) <- c('limma','bambu','isonform','rattle','rnabloom2','trinity','bambudenovo','rnaspades','rnabloom2hybrid')
shapes <- c(15,15,16,17,17,17,17,17,17)
names(shapes) <- c('sim','bambu','corset','isonclust','rattle','trinity','bambudenovo','rnaspades','rnabloom2hybrid')

args <- commandArgs(trailingOnly=TRUE)

args <- c('../simulation_1m/bambu/',
          '../simulation_1m/bambudenovo/',
          '../simulation_1m/isonform/',
          '../simulation_1m/rattle/', 
          '../simulation_1m/rnabloom2/',
          '../simulation_1m/trinity/',
          '../simulation_1m/rnaspades/',
          '../simulation_1m/rnabloom2hybrid/',
          '../sqanti_sim/'
)

dirs <- list.dirs(args, recursive = F, full.names = T)

busco_path <- dirs[grepl('_busco$|corrected$', dirs)]

busco_summary <- lapply(busco_path, function(x) {
  read_busco(x) %>%
    mutate(data = x) %>%
    separate(data, c(NA,NA,NA,'method',NA,'correct'),sep = '/|_', extra = 'merge')
}) 

busco_summary <- busco_summary %>% Reduce(rbind, .) %>%
  mutate(method = factor(method, levels = c('','bambu','rattle','rnabloom2','isonform','trinity','bambudenovo','rnaspades','rnabloom2hybrid'), 
                         labels = c('ref','bambu','rattle','rnabloom2','isonform','trinity','bambudenovo','rnaspades','rnabloom2hybrid')),
    iscorrect = ifelse(grepl('corrected', correct), 'corrected','raw'), 
         Class = factor(Class, levels = busco_summary[[1]]$Class %>% as.character() %>% rev())) 

busco_summary %>% pivot_wider(names_from = Class, values_from = Frequency) %>% 
  arrange(method, correct) %>% 
  write.csv('plot/busco_summary.csv')

busco_summary %>% 
  ggplot(aes(x=method, y=Frequency/13780, fill=Class)) + 
  geom_col(position = 'stack') +
  coord_flip() +
  facet_grid(iscorrect~.) +
  scale_fill_manual(values = rev(c("#56B4E9", "#3492C7", "#F0E442", "#F04442")))
ggsave('plot/busco_summary.pdf', width = 8, height = 5)

busco_summary %>%
  filter(iscorrect == 'corrected') %>%
  ggplot(aes(x=method, y=Frequency/13780, fill=Class)) + 
  geom_col(position = 'stack') +
  coord_flip() +
  facet_grid(iscorrect~.) +
  scale_fill_manual(values = rev(c("#56B4E9", "#3492C7", "#F0E442", "#F04442"))) +
  ggtitle('BUSCO scores (corected)')
ggsave('plot/busco_summary_corrected.pdf', width = 6, height = 3)

busco_summary %>%
  filter(iscorrect == 'raw') %>%
  ggplot(aes(x=method, y=Frequency/13780, fill=Class)) + 
  geom_col(position = 'stack') +
  coord_flip() +
  facet_grid(iscorrect~.) +
  scale_fill_manual(values = rev(c("#56B4E9", "#3492C7", "#F0E442", "#F04442")))+
  ggtitle('BUSCO scores (raw)')
ggsave('plot/busco_summary_raw.pdf', width = 8, height = 3)

# busco_summary %>%
#   tidyr::unite('id', method, iscorrect) %>%
#   filter(id %in% c('ref_raw','bambu_raw','rattle_corrected','rnabloom2_corrected','isonform_corrected','trinity_corrected',
#                    'bambudenovo','rnaspades','rnabloom2hybrid')) %>%
#   mutate(id = factor(id, levels = c('ref_raw','bambu_raw','rattle_corrected','rnabloom2_corrected','isonform_corrected','trinity_corrected',
#                                     ),
#                      labels = c('ref','bambu','rattle','rnabloom2','isonform','trinity'))) %>%
#   ggplot(aes(x=id, y=Frequency/13780, fill=Class)) + 
#   geom_col(position = 'stack') +
#   coord_flip() +
#   scale_fill_manual(values = rev(c("#56B4E9", "#3492C7", "#F0E442", "#F04442")))+
#   ggtitle('BUSCO scores combined (raw for ref, corrected for de novo)') + 
#   xlab('Assembler') +
#   ylab('Proportion')
# ggsave('plot/busco_summary_combined.pdf', width = 5, height = 2)

# 00:31:22
# 2-11:10:07
# 00:54:47
# 02:57:05
# isonform 2-11:05:18

# run_summary <- data.frame(time_min=c(31+22/60, 59*60+10+7/60, 54+47/60, 2*60+57+5/60),
#                           memory_GB=c(28.64, 110.33, 51.56, 123.32),
#                           tx_number=c(47895, 35726, 39488, 67789),
#                           method = c('bambu','rattle', 'rnabloom2','trinity'))

# run_summary <- data.frame(time_min=c(25+35/60, 59*60+10+7/60, 54+47/60, 59*60+5+18/60, 2*60+57+5/60),
#                           memory_GB=c(28.58, 110.33, 51.56, 65.51, 123.32),
#                           tx_number=c(38952, 35726, 39488, 4848, 67789),
#                           method = c('bambu','rattle', 'rnabloom2','isonform','trinity'))
# 
# run_summary %>%
#   ggplot(aes(x = memory_GB, y = time_min/60, 
#              color = method, size = tx_number)) +
#   geom_point() +
#   scale_size_area() +
#   scale_color_manual(values = cols) +
#   scale_y_log10() +
#   scale_x_log10() +
#   ylab('Time (hours)') +
#   xlab('RAM (GB)')
# ggsave('plot/ram_time.pdf', width = 4, height = 3)

