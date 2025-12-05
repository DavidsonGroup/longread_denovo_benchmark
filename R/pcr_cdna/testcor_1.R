
# test 
# check what proportion of dominate transcript in each reference transcript group express
library(data.table)

setwd("/vast/projects/lab_davidson/yan.a/Dong_2023/R_revise")
salmon_count <- readRDS('plot/salmon_count.rds')
sqanti_summary <- readRDS('plot/sqanti_summary.rds')

sqanti_summary <- sqanti_summary[c('bambu_10m', 'bambu_2m','bambu_5m',
                                   rep(c('isonform_2m','isonform_5m'), each = 2), 
                                   rep(c('rattle_10m','rattle_2m','rattle_5m'), each = 3),
                                   rep(c('rnabloom2_10m','rnabloom2_2m','rnabloom2_5m'), each = 2), 
                                   rep('trinity_10m', 2))]

summary2 <- lapply(1:24, function(x) {
  data.frame(counts = salmon_count[[x]]$counts,
             isoform = salmon_count[[x]]$isoform) %>% 
    left_join(sqanti_summary[[x]], by = 'isoform') %>% 
    filter(associated_transcript != 'novel') %>%
    group_by(associated_transcript) %>%
    mutate(percentage = counts / sum(counts) * 100) %>%
    slice_max(percentage, n = 1) %>%
    select(percentage) %>% 
    mutate(data = names(salmon_count)[x])
}) %>% rbindlist()

ggplot(summary2 %>%
         filter(str_detect(data, 'onts')), 
       aes(x = data, y = percentage)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
