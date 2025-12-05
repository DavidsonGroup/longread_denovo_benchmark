
# correlation stratified by Multiplicity (Number of contigs in reference transcript) 
setwd('~/lab_davidson/yan.a/dRNA_gencode/R/')

sqanti_summary <- readRDS('plot/sqanti_summary.rds')
quant_summary <- readRDS('plot/quant_summary.rds')

rt_ref <- sqanti_summary$rattle %>% 
  filter(associated_transcript != 'novel') %>%
  group_by(associated_transcript) %>%
  summarise(n = n()) # for each reference tx, how many de novo matched to it
  
rt_tx <- quant_summary$A549_MCF7_directRNA_quant_rattle_count$lcpm_tx

cor(cpm(rt_tx$true_counts, log = T), cpm(rt_tx$counts.max, log=T))
cor(cpm(rt_tx$true_counts), cpm(rt_tx$counts.max))
cor(rt_tx$true_counts, rt_tx$counts.max)

rb_ref <- sqanti_summary$rnabloom2 %>% 
  filter(associated_transcript != 'novel') %>%
  group_by(associated_transcript) %>%
  summarise(n = n())

rb_tx <- quant_summary$A549_MCF7_directRNA_quant_rnabloom2_onts$lcpm_tx

cor(cpm(rb_tx$true_counts, log = T), cpm(rb_tx$counts.max, log=T))
cor(cpm(rb_tx$true_counts), cpm(rb_tx$counts.max))
cor(rb_tx$true_counts, rb_tx$counts.max)

library(dplyr)
library(ggplot2)

# 1. Prepare RNA-Bloom 2 data
# We rename n() to 'num_denovo' to be explicit
rb_data <- rb_tx %>% 
  left_join(rb_ref, by = c('txid'='associated_transcript')) %>%
  replace_na(list(n = 0)) %>%
  group_by(n) %>%
  summarise(
    correlation = cor(assemble_cpm.max, true_cpm), 
    num_ref = n()
  ) %>% 
  ungroup() %>%
  # head(n=30) %>%
  # Calculate percentage: (num_denovo / total) * 100
  mutate(
    pct_ref = (num_ref / sum(num_ref)) * 100,
    Tool = "RNA-Bloom2"
  )

# 2. Prepare Rattle data
# We also name the count 'num_denovo' here so they align in the plot
rt_data <- rt_tx %>% 
  left_join(rt_ref, by = c('txid'='associated_transcript')) %>%
  replace_na(list(n = 0)) %>%
  group_by(n) %>%
  summarise(
    correlation = cor(assemble_cpm.max, true_cpm), 
    num_ref = n()
  ) %>% 
  ungroup() %>%
  # head(n=30) %>%
  # Calculate percentage: (num_denovo / total) * 100
  mutate(
    pct_ref = (num_ref / sum(num_ref)) * 100,
    Tool = "RATTLE"
  )

# 3. Combine the two datasets
plot_data <- bind_rows(rb_data, rt_data)

# 4. Plot correlation vs Multiplicity (Number of contigs in reference transcript)
ggplot(plot_data, aes(x = n, y = correlation, color = Tool)) +
  geom_line(linewidth = 1) + 
  geom_point(aes(size = num_ref)) +   # Size is based on the count column
  scale_size_continuous(
    trans = "log10",                     # Log scale for the size
    name = "Number of reference\ntranscripts (log10)"
  ) +
  labs(
    y = "Correlation", 
    x = "Multiplicity (Number of contigs in reference transcript)",
    title = "Correlation in each reference transcript group binned by number of de novo transcript"
  ) +
  xlim(c(0,30)) +
  theme_minimal()

# histogram Plot
ggplot(plot_data, aes(x = n, y = pct_ref, fill = Tool)) +
  geom_col(position = "dodge", width = 0.8) +
  labs(
    x = "Multiplicity (Number of contigs in reference transcript)",
    y = "% of reference transcripts",
    title = "Distribution of de novo transcripts by reference transcript group"
  ) +
  xlim(c(-1,20)) + 
  theme_minimal()

ggplot(plot_data, aes(x = n, y = num_ref, fill = Tool)) +
  geom_col(position = "dodge", width = 0.8) +
  labs(
    x = "Multiplicity (Number of contigs in reference transcript)",
    y = "Number of reference transcripts",
    title = "Distribution of de novo transcripts by reference transcript group"
  ) +
  xlim(c(-1,20)) + 
  theme_minimal()







