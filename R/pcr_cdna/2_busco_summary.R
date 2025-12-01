
library(tidyverse)
library(data.table)
library(ggpubr)
library(scales)
source('~/lab_davidson/yan.a/software/scripts_denovo/R/read_busco.R')

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

dirs <- list.dirs(args, recursive = F, full.names = T)

busco_path <- dirs[grepl('_busco$|corrected$', dirs)]

busco_summary <- lapply(busco_path, function(x) {
  read_busco(x) %>%
    mutate(data = x) %>%
    separate(data, c(NA,'method',NA,NA,'depth','correct'),sep = '/|_', extra = 'merge')
}) 

busco_summary <- busco_summary %>% Reduce(rbind, .) %>%
  mutate(method = factor(method, 
                         levels = c('bambu','rattle','rnabloom2','isonform','trinitystranded','bambudenovo','rnaspades','rnabloom2hybrid'),
                         labels = c('bambu','rattle','rnabloom2','isonform','trinity','bambudenovo','rnaspades','rnabloom2hybrid')),
         depth = factor(str_replace(depth, 'ilu|merged','10m'), levels = c('2m','5m','10m')),
         iscorrect = ifelse(grepl('corrected', correct), 'corrected','raw'), 
         Class = factor(Class, levels = busco_summary[[1]]$Class %>% as.character() %>% rev())) 

busco_summary %>% pivot_wider(names_from = Class, values_from = Frequency) %>% 
  arrange(method, correct, depth) %>%
  write.csv('plot/busco_summary.csv')

busco_summary %>% 
  ggplot(aes(x=depth, y=Frequency/13780, fill=Class)) + 
  geom_col(position = 'stack') +
  coord_flip() +
  facet_grid(iscorrect~method) +
  scale_fill_manual(values = rev(c("#56B4E9", "#3492C7", "#F0E442", "#F04442")))
ggsave('plot/busco_summary.pdf', width = 12, height = 3)

busco_summary %>%
  filter(iscorrect == 'corrected') %>%
  ggplot(aes(x=depth, y=Frequency/13780, fill=Class)) + 
  geom_col(position = 'stack') +
  coord_flip() +
  facet_grid(iscorrect~method) +
  scale_fill_manual(values = rev(c("#56B4E9", "#3492C7", "#F0E442", "#F04442"))) +
  ggtitle('BUSCO scores (corected)')
ggsave('plot/busco_summary_corrected.pdf', width = 12, height = 2)

busco_summary %>%
  filter(iscorrect == 'raw') %>%
  ggplot(aes(x=depth, y=Frequency/13780, fill=Class)) + 
  geom_col(position = 'stack') +
  coord_flip() +
  facet_grid(iscorrect~method) +
  scale_fill_manual(values = rev(c("#56B4E9", "#3492C7", "#F0E442", "#F04442")))+
  ggtitle('BUSCO scores (raw)')
ggsave('plot/busco_summary_raw.pdf', width = 12, height = 2)

busco_summary %>%
  filter(method == 'bambu' | iscorrect == 'corrected') %>%
  unite('Assembler', method, depth, remove = F) %>%
  # mutate(id = factor(id, levels = c('ref_raw','bambu_raw','rattle_corrected','rnabloom2_corrected','trinity_corrected'),
  #                    labels = c('ref','bambu','rattle','rnabloom2','trinity'))) %>%
  ggplot(aes(x=method, y=Frequency/13780, fill=Class)) + 
  geom_col(position = 'stack') +
  coord_flip() +
  scale_fill_manual(values = rev(c("#56B4E9", "#3492C7", "#F0E442", "#F04442")))+
  facet_grid(.~depth) +
  ggtitle('BUSCO scores combined (raw for ref, corrected for de novo)') + 
  xlab('Assembler') +
  ylab('Proportion')
ggsave('plot/busco_summary_combined.pdf', width = 10, height = 2)



# library(readxl)
# 
# run_summary <- read_excel('time_dong.xlsx')
# 
# run_summary %>%
#   mutate(time = as.numeric(difftime(time, as.POSIXct('1899-12-31', tz = "UTC"), units = "mins"))/60,
#          depth = factor(depth, levels = c('2m','5m','10m'))) %>%
#   ggplot(aes(x = ram, y = time, 
#              color = method, group = method)) +
#   geom_line() + 
#   geom_point(aes(shape = depth, size = tx_num)) +
#   scale_size_area() +
#   scale_color_manual(values = cols) +
#   scale_y_log10() +
#   scale_x_log10() +
#   ylab('Time (hours)') +
#   xlab('RAM (GB)')
# ggsave('plot/ram_time.pdf', width = 4, height = 3)

