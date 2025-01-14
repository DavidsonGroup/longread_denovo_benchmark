library(ggh4x)

load("/vast/projects/davidson_longread/yan.a/simulation_20240501/R/de_data.RData")
## ROC like curve
## DTU-TX
df_dtutx %>%
  filter(quant %in% c('map','onts','count'),
         clustering %in% c('bambu','sim','corset')) %>%
  ggplot(aes(x = rowid, y = csum, group = method, color = assembler)) + 
  geom_line(aes(linetype = assembler %in% c('trinity','bambu', 'limma')),
            linewidth = 1) +
  scale_color_manual(values = cols) +
  ylab('Unique True Positives') +
  xlab('Top Ranked Assembled Transcripts') + 
  ggtitle('DTU-TX') +
  force_panelsizes(rows = unit(3, "in"),
                   cols = unit(3, "in"))
ggsave('plot/ROC_dtutx_simplified.pdf', width = 7, height = 5)

# DTE-TX
df_dtetx %>%
  filter(quant %in% c('map','onts','count'),
         clustering %in% c('bambu','sim','corset')) %>%
  ggplot(aes(x = rowid, y = csum, group = method, color = assembler)) + 
  geom_line(aes(linetype = assembler %in% c('trinity','bambu', 'limma')),
            linewidth = 1) +
  scale_color_manual(values = cols) +
  ylab('Unique True Positives') +
  xlab('Top Ranked Assembled Transcripts') + 
  ggtitle('DTE-Tx') +
  force_panelsizes(rows = unit(3, "in"),
                   cols = unit(3, "in"))
ggsave('plot/ROC_dtetx_simplified.pdf', width = 7, height = 5)

# DTU-gene
df_dtugene %>%
  filter(quant %in% c('map','onts','count'),
         clustering %in% c('bambu','sim','corset')) %>%
  ggplot(aes(x = rowid, y = csum, group = method, color = assembler)) + 
  geom_line(aes(linetype = assembler %in% c('trinity','bambu', 'limma')),
            linewidth = 1) +
  scale_color_manual(values = cols) +
  ylab('Unique True Positives') +
  xlab('Top Ranked Clusters') + 
  ggtitle('DTU-Gene') +
  force_panelsizes(rows = unit(3, "in"),
                   cols = unit(3, "in"))
ggsave('plot/ROC_dtugene_simplified.pdf', width = 7, height = 5)

# DGE-gene
df_dge %>%
  filter(quant %in% c('map','onts','count'),
         clustering %in% c('bambu','sim','corset')) %>%
  ggplot(aes(x = rowid, y = csum, group = method, color = assembler)) + 
  geom_line(aes(linetype = assembler %in% c('trinity','bambu', 'limma')),
            linewidth = 1) +
  scale_color_manual(values = cols) +
  ylab('Unique True Positives') +
  xlab('Top Ranked Clusters') + 
  ggtitle('DGE-Gene') +
  force_panelsizes(rows = unit(3, "in"),
                   cols = unit(3, "in"))
ggsave('plot/ROC_dgegene_simplified.pdf', width = 7, height = 5)
