
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
          '../trinitystranded/'
)

dirs <- list.dirs(args, recursive = F, full.names = T)

busco_path <- dirs[grepl('_busco$|corrected$', dirs)]

busco_summary <- lapply(busco_path, function(x) {
  read_busco(x) %>%
    mutate(data = x) %>%
    separate(data, c(NA,'method',NA,NA,'correct'),sep = '/|_', extra = 'merge')
}) 

busco_summary <- busco_summary %>% Reduce(rbind, .) %>%
  mutate(method = factor(method, levels = c('bambu','rattle','rnabloom2','trinitystranded'), 
                         labels = c('bambu','rattle','rnabloom2','trinity')),
        iscorrect = ifelse(grepl('correct', correct), 'corrected','raw'),
        Class =factor(Class, levels = busco_summary[[1]]$Class %>% as.character() %>% rev())
         )
busco_summary %>% pivot_wider(names_from = Class, values_from = Frequency) %>% 
  arrange(method, correct) %>% 
  write.csv('plot/busco_summary.csv')

busco_summary %>%
  ggplot(aes(x=method, y=Frequency/5366, fill=Class)) + 
  geom_col(position = 'stack') +
  coord_flip() +
  facet_grid(iscorrect~.) +
  scale_fill_manual(values = rev(c("#56B4E9", "#3492C7", "#F0E442", "#F04442")))
ggsave('plot/busco_summary.pdf', width = 8, height = 3)

busco_summary %>%
  filter(iscorrect == 'corrected') %>%
  ggplot(aes(x=method, y=Frequency/5366, fill=Class)) + 
  geom_col(position = 'stack') +
  coord_flip() +
  facet_grid(iscorrect~.) +
  scale_fill_manual(values = rev(c("#56B4E9", "#3492C7", "#F0E442", "#F04442"))) +
  ggtitle('BUSCO scores (corected)')
ggsave('plot/busco_summary_corrected.pdf', width = 8, height = 2)

busco_summary %>%
  filter(iscorrect == 'raw') %>%
  ggplot(aes(x=method, y=Frequency/5366, fill=Class)) + 
  geom_col(position = 'stack') +
  coord_flip() +
  facet_grid(iscorrect~.) +
  scale_fill_manual(values = rev(c("#56B4E9", "#3492C7", "#F0E442", "#F04442")))+
  ggtitle('BUSCO scores (raw)')
ggsave('plot/busco_summary_raw.pdf', width = 8, height = 2)

busco_summary %>%
  tidyr::unite('id', method,iscorrect) %>%
  filter(id %in% c('bambu_raw','rattle_corrected','rnabloom2_corrected','isonform_corrected','trinity_corrected')) %>%
  mutate(id = factor(id, levels = c('bambu_raw','rattle_corrected','rnabloom2_corrected','isonform_corrected','trinity_corrected'),
                     labels = c('bambu','rattle','rnabloom2','isonform','trinity'))) %>%
  ggplot(aes(x=id, y=Frequency/5366, fill=Class)) + 
  geom_col(position = 'stack') +
  coord_flip() +
  scale_fill_manual(values = rev(c("#56B4E9", "#3492C7", "#F0E442", "#F04442")))+
  ggtitle('BUSCO scores combined (raw for ref, corrected for de novo)') + 
  xlab('Assembler') +
  ylab('Proportion')
ggsave('plot/busco_summary_combined.pdf', width = 5, height = 2)


# 34.49, 00:46:12
# 190.85, 1-10:49:05
# 50.17, 03:49:45
# 123.32, 03:35:00
run_summary <- data.frame(time_min=as.numeric(hms(c('00:46:12','34:49:05','03:49:45','03:35:00')), 'minutes') ,
                          memory_GB=c(34.49, 190.85, 50.17, 123.32),
                          tx_number=c(31496, 35638, 97198, 66502),
                          method = c('bambu','rattle', 'rnabloom2','trinity'))

run_summary %>%
  ggplot(aes(x = memory_GB, y = time_min/60, 
             color = method, size = tx_number)) +
  geom_point() +
  scale_size_area() +
  scale_color_manual(values = cols) +
  ylab('Time (hours)') +
  xlab('RAM (GB)')
ggsave('plot/ram_time.pdf',width = 4, height = 3)

