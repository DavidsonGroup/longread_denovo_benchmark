library(tidyverse)
library(ggpubr)
library(tximport)
library(cowplot)

setwd('/vast/projects/lab_davidson/yan.a/Dong_2023/R_revise')
source("~/lab_davidson/yan.a/software/scripts_denovo/R/quantile_overlap.R")

sqanti_summary <- readRDS("/vast/projects/lab_davidson/yan.a/Dong_2023/R_revise/plot/sqanti_summary.rds")
requested_counts <- readRDS("/vast/projects/lab_davidson/yan.a/Dong_2023/R_revise/plot/requested_counts.rds")

tx_counts_bambu <- readRDS('../bambu/merged_10m_dge/R/y.rds')

tx_counts_bambu <- tx_counts_bambu$counts[rowSums(tx_counts_bambu$counts) != 0, ]

dim(tx_counts_bambu)

gtf <- plyranges::read_gff('../reference/annotation_gencodev44_sequin.gtf') %>%
  dplyr::filter(type == 'transcript')
tx2gene <- gtf %>% data.frame() %>%
  select(transcript_id, gene_id, gene_name, gene_type)

# rnabloom2 in sample1
rb_ref <- sqanti_summary$rnabloom2_10m # 362135 de novo transcripts

count_s1_rb <- tximport('../rnabloom2/merged_10m_dge/barcode01_10m_quant_rnabloom2_onts/quant.sf', 
                     type = "salmon", 
                     txOut = T)$counts
count_s1_rb <- data.frame(counts = count_s1_rb[,1]) %>%
  rownames_to_column(var = 'isoform')

quant1_rb <- quantile_overlap(rb_ref, 
                           count_s1_rb, 'isoform', requested_counts)

rb_tx1 <- quant1_rb$lcpm_tx %>% 
  left_join(data.frame(tx_counts_bambu) %>% rownames_to_column(var = 'txid'),
            by = 'txid') %>%
  filter(str_detect(txid, '^R')) %>%
  mutate(geneid = str_remove(txid, '_[0-9]$'),
         counts.max = ifelse(isassemble=='not assembled', NA, counts.max))

rb_tx1 %>%
  left_join(tx2gene, by = c('txid' = 'transcript_id')) %>%
  filter(!is.na(gene_id)) %>%
  group_by(gene_id) %>% 
  mutate(percentage = counts.max / sum(counts.max, na.rm = T) * 100) %>%
  slice_max(percentage, n = 1) %>%
  pull(percentage) %>%
  summary

# rattle s1
rt_ref <- sqanti_summary$rattle_10m 

count_s1_rt <- tximport('../rattle/merged_10m_dge/barcode01_10m_quant_rattle_onts/quant.sf', 
                     type = "salmon", 
                     txOut = T)$counts
count_s1_rt <- data.frame(counts = count_s1_rt[,1]) %>%
  rownames_to_column(var = 'isoform')

quant1_rt <- quantile_overlap(rt_ref, 
                              count_s1_rt, 'isoform', requested_counts)

rt_tx1 <- quant1_rt$lcpm_tx %>% 
  left_join(data.frame(tx_counts_bambu) %>% rownames_to_column(var = 'txid'),
            by = 'txid') %>%
  filter(str_detect(txid, '^R')) %>%
  mutate(geneid = str_remove(txid, '_[0-9]$'),
         counts.max = ifelse(isassemble=='not assembled', NA, counts.max))

rt_tx1 %>%
  left_join(tx2gene, by = c('txid' = 'transcript_id')) %>%
  filter(!is.na(gene_id)) %>%
  group_by(gene_id) %>% 
  mutate(percentage = counts.max / sum(counts.max, na.rm = T) * 100) %>%
  slice_max(percentage, n = 1) %>%
  pull(percentage) %>%
  summary

gid <- rt_tx1 %>%
  left_join(tx2gene, by = c('txid' = 'transcript_id')) %>%
  filter(!is.na(gene_id)) %>%
  group_by(gene_id) %>% 
  mutate(percentage = counts.max / sum(counts.max, na.rm = T) * 100,
         true_pct = true_counts / sum(true_counts, na.rm = T) * 100,
         pct_diff = true_pct - percentage,
         n_tx = n()) %>% 
  filter(n_tx > 1) %>% 
  slice_max(pct_diff, n = 1) %>% 
  ungroup() %>% 
  arrange(desc(pct_diff)) %>%
  pull(gene_id) %>% head(n=10)

rb_tx1 %>% 
  select(counts.max, barcode01_10m_sorted, txid, geneid) %>%
  left_join(rt_tx1 %>%
              select(counts.max, barcode01_10m_sorted, txid, geneid),
            by = c('txid','geneid','barcode01_10m_sorted'), 
            suffix = c('.rnabloom2','.rattle')) %>%
  group_by(geneid) %>%
  mutate(
    has_rattle_na = any(is.na(counts.max.rattle)),
    has_rnabloom2_na = any(is.na(counts.max.rnabloom2)),
    group = case_when(
      has_rattle_na & has_rnabloom2_na ~ "missed",
      !has_rattle_na & !has_rnabloom2_na ~ "detected",
      has_rattle_na ~ "rattlemissing",
      has_rnabloom2_na ~ "rnabloom2missing"
    )
  ) %>%
  select(-has_rattle_na, -has_rnabloom2_na) %>%  # Remove helper columns
  ungroup() %>% 
  pivot_longer(cols = c(counts.max.rnabloom2, counts.max.rattle, barcode01_10m_sorted),
               values_to = 'counts',
               names_to = 'source') %>%
  #filter(group != 'detected') %>%
  filter(geneid %in% gid) %>%
  ggplot(aes(x = txid, y = counts, colour = source)) +
  geom_point(position = position_jitter(width = 0.2, height = 0)) +
  facet_grid( group ~ geneid, scales = 'free', space = 'free_x') +  # Rows by group, columns by geneid
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6)) +
  scale_colour_manual(
    values = c(
      "counts.max.rnabloom2"                      = "#D55E00",
      "counts.max.rattle"                         = "#009E73",
      "barcode01_10m_sorted" = "#999999"
    ),
    name   = "Source",
    labels = c("RNA-Bloom2", "RATTLE", "Truth")
  ) 







