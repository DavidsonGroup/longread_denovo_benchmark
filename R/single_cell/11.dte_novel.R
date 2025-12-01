library(ggh4x)
library(tidyverse)

setwd("~/lab_davidson/yan.a/pb/5m_stranded/R/")

args <- c('~/vast_scratch/sqanti3_bambu/sqanti3_results/pb_sc/isonform/',
  '~/vast_scratch/sqanti3_bambu/sqanti3_results/pb_sc/rattle/',
  '~/vast_scratch/sqanti3_bambu/sqanti3_results/pb_sc/rnabloom2/'
)

dirs <- list.dirs(args, recursive = F, full.names = T)

sqanti_files <- list.files(dirs, 'classification.txt', 
                           recursive = T, full.names = T)

sqanti_summary2 <- lapply(sqanti_files, function(x) {
  read.table(x, sep = '\t', header = T) %>%
    dplyr::select(c(1:8,15:17,25,29,30,41,46)) %>%
    dplyr::mutate(original = ifelse(str_detect(isoform, '_dup'), str_remove(isoform, '_dup.*'), isoform),
                  category = case_when(
                    # Logic: Is my root in the list of roots that have a duplicate somewhere in the WHOLE table?
                    original %in% original[str_detect(isoform, "_dup")] ~ 'chimera',
                    str_detect(associated_transcript, '^Bambu') ~ 'novel_in_bambu',
                    str_detect(associated_transcript, '^novel') | is.na(associated_transcript) ~ 'novel',
                    TRUE ~ 'reference'
                  ))
})
names(sqanti_summary2) <- c('isonform',
  'rattle',
  'rnabloom2')

dte <- lapply(list.files(pattern = 'dte_tx.rds', full.names = T, recursive = T), 
              function(x){
                readRDS(x) %>% 
                  filter(adj.P.Val<0.05) 
              })
names(dte) <- list.files(pattern = 'dte_tx.rds', full.names = T, recursive = T) %>% dirname() %>% str_remove('./')
dte <- dte[ c('isonform_isonclust_oarfish',
  'rattle_rattle_oarfish', 
  'rnabloom2_corset_oarfish')]

df <- lapply(1:length(dte), function(x) {
  dte[[x]] %>% left_join(sqanti_summary2[[x]], by = 'isoform',
                         suffix = c('.gencode','.bambu')) %>%
    mutate(category = ifelse(is.na(category), 'unmapped', category)) %>% # umapped is still in dte, but not annotated in sqanti3
    count(category)
}) %>% purrr::reduce(left_join, by = 'category') 

colnames(df) <- c('category', names(sqanti_summary2) )
write.csv(df, 'plot/dte_novel_in_bambu.csv')

df %>%
  pivot_longer(-1,
               names_to = "assembler",
               values_to = "count") %>%
  ggplot(aes(x = assembler, y = count, fill = category)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal(base_size = 14) +
  labs(x = NULL, y = "Count", fill = "Category") +
  scale_fill_brewer(palette = "Set2") +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  force_panelsizes(cols = unit(3, 'inch'), 
                   rows = unit(3, 'inch'))
ggsave('plot/dte_novel_in_bambu.pdf', width = 5, height = 5)
