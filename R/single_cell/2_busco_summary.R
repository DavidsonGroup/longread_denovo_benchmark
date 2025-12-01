
library(tidyverse)
library(data.table)
library(ggpubr)
source('~/lab_davidson/yan.a/software/scripts_denovo/R/read_busco.R')

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
          '../isonform/'
)

dirs <- list.dirs(args, recursive = F, full.names = T)

busco_path <- dirs[grepl('_busco$|corrected$', dirs)]

busco_summary <- lapply(busco_path, function(x) {
  read_busco(x) %>%
    mutate(data = x) %>%
    separate(data, c(NA, 'method1','long',NA,'correct'),sep = '/|_', extra = 'merge') %>%
    unite(method, method1, long, sep = "")
}) 

busco_summary <- busco_summary %>% Reduce(rbind, .) %>%
  mutate(method = factor(method, levels = c('bambu','rattle','rnabloom2','isonform'), 
                         labels = c('bambu','rattle','rnabloom2','isonform')),
        iscorrect = ifelse(grepl('correct', correct), 'corrected','raw'),
        Class =factor(Class, levels = busco_summary[[1]]$Class %>% as.character() %>% rev())
         )

busco_summary %>% pivot_wider(names_from = Class, values_from = Frequency) %>% 
  arrange(method, correct) %>% 
  write.csv('plot/busco_summary.csv')

busco_summary %>%
  ggplot(aes(x=method, y=Frequency/13780, fill=Class)) + 
  geom_col(position = 'stack') +
  coord_flip() +
  facet_grid(iscorrect~.) +
  scale_fill_manual(values = rev(c("#56B4E9", "#3492C7", "#F0E442", "#F04442")))
ggsave('plot/busco_summary.pdf', width = 8, height = 3)

busco_summary %>%
  filter(iscorrect == 'corrected') %>%
  ggplot(aes(x=method, y=Frequency/13780, fill=Class)) + 
  geom_col(position = 'stack') +
  coord_flip() +
  facet_grid(iscorrect~.) +
  scale_fill_manual(values = rev(c("#56B4E9", "#3492C7", "#F0E442", "#F04442"))) +
  ggtitle('BUSCO scores (corected)')
ggsave('plot/busco_summary_corrected.pdf', width = 8, height = 2)

busco_summary %>%
  filter(iscorrect == 'raw') %>%
  ggplot(aes(x=method, y=Frequency/13780, fill=Class)) + 
  geom_col(position = 'stack') +
  coord_flip() +
  facet_grid(iscorrect~.) +
  scale_fill_manual(values = rev(c("#56B4E9", "#3492C7", "#F0E442", "#F04442")))+
  ggtitle('BUSCO scores (raw)')
ggsave('plot/busco_summary_raw.pdf', width = 8, height = 2)

busco_summary %>%
  tidyr::unite('id', method,iscorrect) %>%
  filter(id %in% c('bambu_raw','rattle_corrected','rnabloom2_corrected','isonform_corrected')) %>%
  mutate(id = factor(id, levels = c('bambu_raw','rattle_corrected','rnabloom2_corrected','isonform_corrected'),
                     labels = c('bambu','rattle','rnabloom2','isonform'))) %>%
  ggplot(aes(x=id, y=Frequency/13780, fill=Class)) + 
  geom_col(position = 'stack') +
  coord_flip() +
  scale_fill_manual(values = rev(c("#56B4E9", "#3492C7", "#F0E442", "#F04442")))+
  ggtitle('BUSCO scores combined (raw for ref, corrected for de novo)') + 
  xlab('Assembler') +
  ylab('Proportion')
ggsave('plot/busco_summary_combined.pdf', width = 5, height = 2)


# run_summary <- data.frame(time_min=as.numeric(hms(c('00:56:34','24:80:74','19:56:23', '03:44:32', '05:06:29')), 'minutes') + c(0, 3*60*24,3*60*24,1*60*24,0),
#                           memory_GB=c(39.08, 125, 217,198,123),
#                           tx_number=c(69478, 28130, 54954, 443874, 163971),
#                           method = c('bambu','isonform','rattle', 'rnabloom2','trinity'))

run_summary <- data.frame(time_min=as.numeric(hms(c('3:16:01','12:33:56', '02:27:00')), 'minutes') + 
                            c(3*60*24,4*60*24,0),
                          memory_GB=c( 119.60, 221.23, 67.99),
                          tx_number=c( 242488, 143901, 508641),
                          method = c('isonform','rattle', 'rnabloom2'))

run_summary %>%
  ggplot(aes(x = memory_GB, y = time_min/60, 
             color = method, size = tx_number)) +
  geom_point() +
  scale_size_area() +
  scale_color_manual(values = cols) +
  ylab('Time (hours)') +
  xlab('RAM (GB)')
ggsave('plot/ram_time.pdf',width = 4, height = 3)

