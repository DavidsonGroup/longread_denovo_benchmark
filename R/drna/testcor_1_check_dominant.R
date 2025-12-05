# test 
# check what proportion of dominate transcript in each reference transcript group express
library(data.table)

salmon_count <- readRDS('plot/salmon_count.rds')
sqanti_summary <- readRDS('plot/sqanti_summary.rds')

sqanti_summary <- sqanti_summary[c('isonform','isonform','rattle','rattle','rattle','rnabloom2','rnabloom2','trinity','trinity')]

summary <- lapply(1:9, function(x) {
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

ggplot(summary %>%
         filter(str_detect(data, 'onts')), 
       aes(x = data, y = percentage)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


## NEED CHECK
count <- salmon_count$A549_MCF7_directRNA_quant_rattle_onts
sq <- sqanti_summary$rattle

data.frame(counts = count$counts,
           isoform = rownames(count)) %>% 
  left_join(sq, by = 'isoform') %>% 
  filter(associated_transcript != 'novel') %>%
  group_by(associated_transcript) %>%
  mutate(percentage = counts / sum(counts) * 100) %>%
  slice_max(percentage, n = 1) %>%
  pull(percentage) %>% length 



count <- readRDS('../../Dong_2023/R/plot/salmon_count.rds')
sq <- readRDS('../../Dong_2023/R/plot/sqanti_summary.rds')

count <- count$merged_10m_quant_rnabloom2_onts
sq <- sq$rnabloom2_10m

count <- count$merged_10m_quant_rattle_onts
sq <- sq$rattle_10m

data.frame(counts = count$counts,
           isoform = rownames(count)) %>% 
  left_join(sq, by = 'isoform') %>% 
  filter(associated_transcript != 'novel') %>%
  group_by(associated_transcript) %>%
  mutate(percentage = counts / sum(counts) * 100) %>%
  slice_max(percentage, n = 1) %>%
  pull(percentage) %>% summary 


cor(quant_summary$A549_MCF7_directRNA_quant_rattle_onts$lcpm_tx$assemble_cpm.sum, 
    quant_summary$A549_MCF7_directRNA_quant_rattle_onts$lcpm_tx$true_cpm )

cor(quant_summary$A549_MCF7_directRNA_quant_rattle_onts$lcpm_tx$assemble_cpm.max, 
    quant_summary$A549_MCF7_directRNA_quant_rattle_onts$lcpm_tx$true_cpm )

cor(quant_summary$A549_MCF7_directRNA_quant_rnabloom2_onts$lcpm_tx$assemble_cpm.sum, 
    quant_summary$A549_MCF7_directRNA_quant_rnabloom2_onts$lcpm_tx$true_cpm )

cor(quant_summary$A549_MCF7_directRNA_quant_rnabloom2_onts$lcpm_tx$assemble_cpm.max, 
    quant_summary$A549_MCF7_directRNA_quant_rnabloom2_onts$lcpm_tx$true_cpm )



cor(all_quant$rattle.sec$tx_exp_tiled$counts.max, tx_counts_bambu$counts %>% c)
cor(all_quant$rnabloom2.sec$tx_exp_tiled$counts.max, tx_counts_bambu$counts %>% c)

cor( log2(all_quant$rattle.sec$tx_exp_tiled$counts.max + 1), log2(tx_counts_bambu$counts %>% c +1))
cor( log2(all_quant$rnabloom2.sec$tx_exp_tiled$counts.max + 1), log2(tx_counts_bambu$counts %>% c +1))

cor(all_quant$rattle.sec$tx_exp_tiled$assemble_cpm.max, truecpm_tx)
cor(all_quant$rnabloom2.sec$tx_exp_tiled$assemble_cpm.max, truecpm_tx)

cor( log2(all_quant$rattle.sec$tx_exp_tiled$counts.max + .1), log2(tx_counts_bambu$counts %>% c +.1))
cor( log2(all_quant$rnabloom2.sec$tx_exp_tiled$counts.max + .1), log2(tx_counts_bambu$counts %>% c +.1))

cor(all_quant$rattle.sec$tx_exp_tiled$counts.max, tx_counts_bambu$counts %>% c, method = 'spearman')
cor(all_quant$rnabloom2.sec$tx_exp_tiled$counts.max, tx_counts_bambu$counts %>% c, method = 'spearman')


data.frame(rnabloom2_lcpm = all_quant$rnabloom2.sec$tx_exp_tiled$assemble_cpm, 
           true_lcpm = truecpm_tx) %>% 
  ggplot(aes(x = rnabloom2_lcpm, y = true_lcpm)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1)

