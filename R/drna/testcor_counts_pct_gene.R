library(tidyverse)
library(ggpubr)
library(tximport)
library(cowplot)

setwd('/vast/projects/lab_davidson/yan.a/dRNA_gencode/R/')
source("~/lab_davidson/yan.a/software/scripts_denovo/R/quantile_overlap.R")

sqanti_summary <- readRDS("plot/sqanti_summary.rds")
requested_counts <- readRDS("plot/requested_counts.rds")

tx_counts_bambu <- readRDS('../bambu/A549_MCF7_directRNA_dge/R/y.rds')

tx_counts_bambu <- tx_counts_bambu$counts[rowSums(tx_counts_bambu$counts) != 0, ]

dim(tx_counts_bambu)

gtf <- plyranges::read_gff('~/lab_davidson/yan.a/ref/gencode/gencode.v44.annotation.gtf') %>%
  dplyr::filter(type == 'transcript')
tx2gene <- gtf %>% data.frame() %>%
  select(transcript_id, gene_id, gene_name, gene_type)

# rnabloom2 in sample1
rb_ref <- sqanti_summary$rnabloom2 # 362135 de novo transcripts

count_s1_rb <- tximport('../rnabloom2/A549_MCF7_directRNA_dge/SGNex_A549_directRNA_replicate4_run1_filtered_quant_rnabloom2_onts/quant.sf', 
                        type = "salmon", 
                        txOut = T)$counts
count_s1_rb <- data.frame(counts = count_s1_rb[,1]) %>%
  rownames_to_column(var = 'isoform')

quant1_rb <- quantile_overlap(rb_ref, 
                              count_s1_rb, 'isoform', requested_counts)

rb_tx1 <- quant1_rb$lcpm_tx %>% 
  left_join(data.frame(tx_counts_bambu) %>% rownames_to_column(var = 'txid'),
            by = 'txid') %>%
  # filter(str_detect(txid, '^R')) %>%
  mutate(
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
rt_ref <- sqanti_summary$rattle

count_s1_rt <- tximport('../rattle/A549_MCF7_directRNA_dge/SGNex_A549_directRNA_replicate4_run1_filtered_quant_rattle_onts/quant.sf', 
                        type = "salmon", 
                        txOut = T)$counts
count_s1_rt <- data.frame(counts = count_s1_rt[,1]) %>%
  rownames_to_column(var = 'isoform')

quant1_rt <- quantile_overlap(rt_ref, 
                              count_s1_rt, 'isoform', requested_counts)

rt_tx1 <- quant1_rt$lcpm_tx %>% 
  left_join(data.frame(tx_counts_bambu) %>% rownames_to_column(var = 'txid'),
            by = 'txid') %>%
  # filter(str_detect(txid, '^R')) %>%
  mutate(
         counts.max = ifelse(isassemble=='not assembled', NA, counts.max))

rt_tx1 %>%
  left_join(tx2gene, by = c('txid' = 'transcript_id')) %>%
  filter(!is.na(gene_id)) %>%
  group_by(gene_id) %>% 
  mutate(percentage = counts.max / sum(counts.max, na.rm = T) * 100,
         n_tx = n()) %>% 
  filter(n_tx > 1) %>% 
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



gid <- c('ENSG00000107281.10','ENSG00000247077.8','ENSG00000250317.9')

rb_tx1 %>% left_join(tx2gene, by = c('txid' = 'transcript_id')) %>% 
  select(counts.max, SGNex_A549_directRNA_replicate4_run1_sorted, txid, gene_id) %>%
  left_join(rt_tx1 %>% left_join(tx2gene, by = c('txid' = 'transcript_id')) %>%
              select(counts.max, SGNex_A549_directRNA_replicate4_run1_sorted, txid, gene_id),
            by = c('txid','gene_id','SGNex_A549_directRNA_replicate4_run1_sorted'), suffix = c('.rnabloom2','.rattle')) %>% 
  filter(gene_id %in% gid) %>%
  pivot_longer(cols = c(counts.max.rnabloom2, counts.max.rattle, SGNex_A549_directRNA_replicate4_run1_sorted),
               values_to = 'counts',
               names_to = 'source') %>%
  mutate(source = factor(source, levels = c("counts.max.rnabloom2", "counts.max.rattle" ,"SGNex_A549_directRNA_replicate4_run1_sorted"))) %>%
  ggplot(aes(x = txid, y = counts, colour = source)) +
  geom_point(position = position_jitter(width = 0.2, height = 0)) +
  facet_wrap( . ~ gene_id, scales = 'free', nrow = 1) +  # Rows by group, columns by geneid
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6)) +
  scale_colour_manual(
    values = c(
      "counts.max.rnabloom2"                      = "#D55E00",
      "counts.max.rattle"                         = "#009E73",
      "SGNex_A549_directRNA_replicate4_run1_sorted" = "#999999"
    ),
    name   = "Source",
    labels = c("RNA-Bloom2", "RATTLE", "Truth")
  ) +
  force_panelsizes(rows = unit(2, 'inch'), cols = unit(2, 'inch'))

ggsave('plot/rattle_fail.pdf', width = 10, height = 7)




rt_diff <- rt_tx1 %>%
  left_join(tx2gene, by = c('txid' = 'transcript_id')) %>%
  filter(!is.na(gene_id)) %>%
  group_by(gene_id) %>% 
  mutate(n_tx = n()) %>%
  filter(n_tx > 1, sum(counts.max, na.rm = T) > 0) %>%
  mutate(percentage = counts.max / sum(counts.max, na.rm = T) * 100,
         true_pct = true_counts / sum(true_counts, na.rm = T) * 100
  ) %>% 
  replace_na(list(percentage=0, true_pct=0)) %>%
  mutate(pct_diff = percentage - true_pct) %>%
  slice_max(pct_diff, n = 1, with_ties = F) %>%
  select(percentage, true_pct,pct_diff, n_tx) 

rb_diff <- rb_tx1 %>%
  left_join(tx2gene, by = c('txid' = 'transcript_id')) %>%
  filter(!is.na(gene_id)) %>%
  group_by(gene_id) %>% 
  mutate(n_tx = n()) %>%
  filter(n_tx >1, sum(counts.max, na.rm = T) > 0) %>%
  mutate(percentage = counts.max / sum(counts.max, na.rm = T) * 100,
         true_pct = true_counts / sum(true_counts, na.rm = T) * 100
  ) %>% 
  replace_na(list(percentage=0, true_pct=0)) %>%
  mutate(pct_diff = percentage - true_pct) %>%
  slice_max(pct_diff, n = 1, with_ties = F) %>%
  select(percentage, true_pct,pct_diff, n_tx) 

left_join(rt_diff, rb_diff, by = 'gene_id',
          suffix = c('.rattle','.rnabloom2') ) %>%
  select(pct_diff.rattle, pct_diff.rnabloom2) %>%
  pivot_longer(
    cols = c(pct_diff.rattle, pct_diff.rnabloom2),
    names_to = "method",
    values_to = "pct_diff"
  ) %>%
  mutate(method = dplyr::recode(
    method,
    "pct_diff.rattle"    = "RATTLE",
    "pct_diff.rnabloom2" = "RNA-Bloom2"
  )) %>%
  ggplot( aes(x = pct_diff, fill = method, colour = method)) +
  geom_histogram(alpha = 0.4, position = "identity", bins = 30) +
  scale_fill_manual(
    values = c(
      "RNA-Bloom2" = "#D55E00",
      "RATTLE"     = "#009E73"
    )
  ) +
  scale_color_manual(
    values = c(
      "RNA-Bloom2" = "#D55E00",
      "RATTLE"     = "#009E73"
    )
  ) +
  labs(x = "Percentage difference", y = "Count", fill = "Method", colour = "Method") +
  theme_bw()
ggsave('plot/deviation_from_truthpct.pdf', width = 5, height = 3)



# rb_tx1 %>% filter(isassemble=='not assembled') %>% pull(SGNex_A549_directRNA_replicate4_run1_sorted) %>% sum
# sum(rb_tx1$counts.max, na.rm = T)
# sum(rb_tx1$counts.sum, na.rm = T)
# sum(rb_tx1$SGNex_A549_directRNA_replicate4_run1_sorted, na.rm = T)
# 
# rt_tx1 %>% filter(isassemble=='not assembled') %>% pull(SGNex_A549_directRNA_replicate4_run1_sorted) %>% sum
# sum(rt_tx1$counts.max, na.rm = T)
# sum(rt_tx1$counts.sum, na.rm = T)
# sum(rt_tx1$SGNex_A549_directRNA_replicate4_run1_sorted, na.rm = T)

data.frame(total_counts_rb = count_s1_rb$counts %>% sum,
           missedtx_counts_rb = rb_tx1 %>% filter(isassemble=='not assembled') %>% pull(SGNex_A549_directRNA_replicate4_run1_sorted) %>% sum,
           total_counts_rt = count_s1_rt$counts %>% sum,
           missedtx_counts_rt = rt_tx1 %>% filter(isassemble=='not assembled') %>% pull(SGNex_A549_directRNA_replicate4_run1_sorted) %>% sum 
) %>% pivot_longer(1:4) %>%
  separate(name, c('group',NA, 'assembler')) %>%
  ggplot(aes(x = assembler, y = value, fill = group)) +
  geom_col(position = "identity", width = 0.6) + 
  # Use distinct colors: a neutral one for Total, a bright one for Missed
  scale_fill_manual(values = c("total" = "#e0e0e0", "missedtx" = "#377eb8")) +
  scale_y_continuous(labels = scales::comma) + # Removes scientific notation
  theme_minimal() +
  scale_x_discrete(labels = c("rb" = "RNA-Bloom2", "rt" = "RATTLE")) +
  labs(title = "dRNA",
       y = "Read counts",
       x = "Assembler")
ggsave('plot/missed_pct.pdf', width = 4, height = 3)
