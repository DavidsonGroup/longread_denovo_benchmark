library(UpSetR)
library(ggpubr)
library(rtracklayer)
library(ggVennDiagram)
library(plyranges)
library(tidyverse)
library(scales)
library(cowplot)
library(ggh4x)
source('~/lab_davidson/yan.a/software/scripts_denovo/R/get_precision_recall.R')
source('~/lab_davidson/yan.a/software/scripts_denovo/R/roc.R')

# set plotting parameters
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#c51b8a")
cols <- cbPalette[c(1,2,6,4,7,8, 3,5,9)]
names(cols) <- c('limma','bambu','isonform','rattle','rnabloom2','trinity','bambudenovo','rnaspades','rnabloom2hybrid')
shapes <- c(15,15,16,17,17,17,17,17,17)
names(shapes) <- c('sim','bambu','corset','isonclust','rattle','trinity','bambudenovo','rnaspades','rnabloom2hybrid')

# load DE from bambu_10m
bambu_anno <- readRDS('bambu_bambu_10m/summary.rds')

true_dte <- readRDS('../bambu/merged_10m_dge/R/dte_tx.rds') %>% 
  filter(adj.P.Val < 0.05 & chrom != 'chrIS') %>% pull(transcript_id) %>%
  na.omit()
true_dge <- readRDS('../bambu/merged_10m_dge/R/dge_gene.rds') %>% filter(adj.P.Val < 0.05) %>%
  left_join(bambu_anno$asm_map, by = c('genes' = 'geneid')) %>% filter(!str_detect(assigned_gene, '^R')) %>%
  pull(assigned_gene) %>%
  na.omit()
true_dtutx <- readRDS('../bambu/merged_10m_dge/R/dtu_tx.rds') %>% 
  filter(FDR < 0.05  & chrom != 'chrIS') %>% pull(transcript_id) %>%
  na.omit()
true_dtugene <- readRDS('../bambu/merged_10m_dge/R/dtu_gene.rds') %>% filter(FDR < 0.05) %>%   
  left_join(bambu_anno$asm_map, by = c('gene_id' = 'geneid')) %>% filter(!str_detect(assigned_gene, '^R')) %>%
  pull(assigned_gene) %>%
  na.omit()

## evaluate DTU at tx and gene level
files.dtu <- list.files('.', pattern = 'dtu_tx.rds', recursive = T)

ntest <- length(files.dtu)

test <- lapply(files.dtu, function(x){
  readRDS(x) %>% filter(chrom != 'chrIS')
  }) 
names(test) <- dirname(files.dtu)

pdf('plot/dtu_tx_upset.pdf', width = 8, height = 10)

dtu.list <- list()
dtu.list[[1]] <- true_dtutx

dtu.list[2:(ntest+1)] <- lapply(1:ntest, function(x){
  method <- dirname(files.dtu)[x]
  #txmethod = str_split(method, '_', simplify = T)[1] 
  txmethod = paste(str_split(method, '_', simplify = T)[c(1,3)] , collapse = '_') 
  
  ids <- test[[x]] %>% 
    rownames_to_column() %>% filter(FDR<0.05) %>%
    mutate(associated_transcript = ifelse(associated_transcript == 'novel', 
                                          # txmethod,
                                          paste(txmethod, rowname, sep = '_'),
                                          associated_transcript)) %>% 
    pull(associated_transcript)
}) 
 
# dtu.list[[12]] <- readRDS('../sqanti_sim/design/limma_dtutx.rds') %>%
#   filter(FDR<0.05) %>% pull(transcript_id)

names(dtu.list) <- c('bambu_bambu_10m_count', dirname(files.dtu))

fromList(dtu.list) %>%
  upset(set_size.show = T, keep.order = T, 
        sets = names(dtu.list),
        nintersects = 50, order.by = 'freq')
dev.off() # actually overlap can be more between assemblies

files <- list.files('.', pattern = 'summary.rds', recursive = T)
files <- files[grepl('/summary', files)]
filelist <- lapply(paste0(str_replace(files.dtu, "^(.*)_[^/]+/dtu_tx\\.rds$", "\\1"), '/summary.rds'), function(x){
  readRDS(x) 
})
names(filelist) <- names(test)

files.dtugene <- list.files('.', pattern = 'dtu_gene.rds', recursive = T)

test2 <- lapply(files.dtugene, function(x){
  readRDS(x) %>% 
    left_join(filelist[[str_remove(x, '/dtu_gene.rds')]]$asm_map, 
              by = 'geneid') %>% 
    filter(!str_detect(assigned_gene, '^R'))
}) 
names(test2) <- dirname(files.dtugene)

pdf('plot/dtu_gene_upset.pdf', width = 8, height = 10)

dtugene.list <- list()
dtugene.list[[1]] <- true_dtugene

dtugene.list[2:(ntest+1)] <- lapply(1:ntest, function(x){
  if (x %in% 1:(ntest)) {
    test2[[x]] %>% 
      rownames_to_column() %>% 
      filter(FDR<0.05) %>%
      # left_join(filelist[[x]]$asm_map, by= 'geneid') %>% 
      pull(assigned_gene)
  } else { 
    test2[[x]] %>% rownames_to_column() %>% 
      filter(FDR<0.05) %>%
      mutate(assigned_gene = gene_id) %>%  # for bambu use geneid
      pull(assigned_gene)
  }
})
# dtugene.list[[12]] <- readRDS('../sqanti_sim/design/limma_dtugene.rds') %>%
#   filter(FDR<0.05) %>% pull(gene_id)

names(dtugene.list) <- c('bambu_bambu_10m_count', dirname(files.dtu))

fromList(dtugene.list) %>%
  upset(set_size.show = T, keep.order = T, 
        sets = names(dtugene.list),
        nintersects = 50, order.by = 'freq')

dev.off()


## evaluate DTE and DGE

files.dte <- list.files('.', pattern = 'dte_tx.rds', recursive = T)
test3 <- lapply(files.dte, function(x){
  readRDS(x) %>% filter(chrom != 'chrIS')
}) 
names(test3) <- dirname(files.dte)

pdf('plot/dte_tx_upset.pdf', width = 8, height = 10)

dtetx.list <- list()
dtetx.list[[1]] <- true_dte

dtetx.list[2:(ntest+1)] <- lapply(1:ntest, function(x){
  
  method <- dirname(files.dtu)[x]
  #txmethod = str_split(method, '_', simplify = T)[1] 
  txmethod = paste(str_split(method, '_', simplify = T)[c(1,3)] , collapse = '_') 

  ids <- test3[[x]] %>% 
    rownames_to_column() %>% 
    filter(adj.P.Val < 0.05) %>%
    mutate(associated_transcript = ifelse(associated_transcript == 'novel', 
                                          paste(txmethod, rowname,  ## different clustering should make no diff
                                                sep = '_'),
                                          associated_transcript)) %>% 
    pull(associated_transcript)
  
})
# dtetx.list[[12]] <- readRDS('../sqanti_sim/design/limma_dtetx.rds') %>%
#   filter(adj.P.Val<0.05) %>% pull(transcript_id)

names(dtetx.list) <- c('bambu_bambu_10m_count', dirname(files.dte))

# DTE is independent of clustering, but dependent on quantification

fromList(dtetx.list[!(lengths(dtetx.list) %>% duplicated())]) %>%
  upset(set_size.show = T, keep.order = T, 
        sets = names(dtetx.list[!(lengths(dtetx.list) %>% duplicated())]),
        nintersects = 50, order.by = 'freq')

dev.off()

files.dge <- list.files('.', pattern = 'dge_gene.rds', recursive = T)
test4 <- lapply(files.dge, function(x){
  readRDS(x) %>% rownames_to_column(var = 'geneid') %>% 
    left_join(filelist[[str_remove(x, '/dge_gene.rds')]]$asm_map, 
              by = 'geneid') %>% 
    filter(!str_detect(assigned_gene, '^R'))
}) 
names(test4) <- dirname(files.dge)

pdf('plot/dge_gene_upset.pdf', width = 8, height = 10)

dge.list <- list()
dge.list[[1]] <- true_dge

dge.list[2:(ntest+1)] <- lapply(1:ntest, function(x){
  if (x %in% (1:ntest)) {
    test4[[x]] %>% 
      # rownames_to_column() %>% 
      filter(adj.P.Val < 0.05) %>%
      # left_join(filelist[[x]]$asm_map, by= c('rowname' = 'geneid')) %>% 
      pull(assigned_gene)
  } else { # for bambu
    test4[[x]] %>% 
      rownames_to_column(var = 'assigned_gene') %>% 
      filter(adj.P.Val < 0.05) %>%
      # left_join(filelist[[x]]$asm_map, by= c('rowname' = 'geneid')) %>% 
      pull(assigned_gene)
  }
})

# dge.list[[12]] <- readRDS('../sqanti_sim/design/limma_dgegene.rds') %>%
#   filter(adj.P.Val<0.05) %>% pull(gene_id)

names(dge.list) <-  c('bambu_bambu_10m_count', dirname(files.dge))

fromList(dge.list) %>%
  upset(set_size.show = T, keep.order = T, 
        sets = names(dge.list),
        nintersects = 50, order.by = 'freq')

dev.off()

plot_pr <- function(delist) {
  
  df.pr <- sapply(1:45, function(x) {
    
    get_precision_recall(delist[[1]],
                         delist[[x+1]])
  }, simplify = T)
  
  colnames(df.pr) <- names(delist)[-1]
  
  p1 <- df.pr %>% t %>%
    data.frame() %>%
    rownames_to_column() %>% 
    separate(rowname, c('assembler','clustering','depth','quant')) %>%
    mutate(depth = factor(depth, levels = c('2m','5m','10m'))) %>%
    filter(quant %in% c('map','onts','sum')) %>%
    ggscatter(x = 'Precision', y = 'Recall', 
              shape = 'clustering', color = 'assembler', 
              # fill = 'clustering', alpha = 'clustering',
              size = 5,
              label = 'clustering', repel = T
              ) +
    theme(legend.position = "top") + 
    scale_color_manual(values = cols) +
    scale_shape_manual(values = shapes) +
    facet_wrap(.~depth) +
    guides(fill=guide_legend(nrow=2,byrow=TRUE)) #+
    # xlim(c(0,1)) +
    # ylim(c(0,1))
  
  p2 <- df.pr %>% t %>%
    data.frame() %>%
    select(4:6) %>%
    rownames_to_column() %>% 
    separate(rowname, c('assembler','clustering','depth','quant')) %>%
    mutate(depth = factor(depth, levels = c('2m','5m','10m'))) %>%
    filter(quant %in% c('map','onts','sum')) %>%
    # dplyr::select(1:6) %>%
    pivot_longer(5:7) %>% 
    ggdotchart(x='name', y = 'value', #facet.by = 'assembler',
               dot.size = 3, sorting = 'none',
               add = "segment", color = 'assembler', shape = 'clustering',
               rotate = TRUE)  +
    scale_color_manual(values = cols) +
    scale_shape_manual(values = shapes) +
    facet_wrap(.~depth)
  
  p <- cowplot::plot_grid(p1, p2)
  
  return(p)
}

plot_pr(dtu.list)
ggsave('plot/dtu_tx_pr.pdf', width = 18, height = 5)

plot_pr(dtugene.list)
ggsave('plot/dtu_gene_pr.pdf', width = 18, height = 5)

plot_pr(dtetx.list)
ggsave('plot/dtx_tx_pr.pdf', width = 18, height = 5)

plot_pr(dge.list)
ggsave('plot/dge_gene_pr.pdf', width = 18, height = 5)

line <- c('2m' = 0.1, '5m' = 0.5, '10m' = 1)

## ROC like curve
## DTU-TX
df_dtutx_all <- lapply(names(test), 
                       function(x) {
                         dt_roc_all(test[[x]], x, truth = true_dtutx)
                       }) %>%
  Reduce(rbind,. ) 

df_dtutx_dedup <- lapply(names(test), 
                         function(x) {
                           dt_roc_dedup(test[[x]], x, truth = true_dtutx)
                         }) %>%
  Reduce(rbind,. ) 

df_dtutx_all %>%
  separate(method, c('assembler','clustering','depth','quant'), remove = F) %>%
  filter(quant %in% c('map','onts','count','sum'),
         depth %in% c('5m', '10m'),
         # For trinity assembler, keep only trinity clustering
         (assembler == "trinity" & clustering == "trinity") |
           # For non-trinity assemblers, keep only bambu, sim, or corset clustering
           (assembler != "trinity" & clustering %in% c("bambu", "sim", "corset"))) %>%
  mutate(lr_denovo = factor(assembler %in% c('trinity','bambu', 'limma'), 
                            labels = c('lr_denovo', 'others'))) %>% 
  ggplot(aes(x = rowid, y = csum, group = method, color = assembler)) + 
  geom_line(aes(linetype = lr_denovo, linewidth = depth)) +
  scale_color_manual(values = cols) +
  scale_linewidth_manual(values = line) +
  ylab('Unique True Positives') +
  xlab('Top Ranked Assembled Transcripts') + 
  ggtitle('DTU-TX') +
  force_panelsizes(rows = unit(3, "in"),
                   cols = unit(3, "in"))
ggsave('plot/ROC_dtutx_simplified.pdf', width = 7, height = 5)

df_dtutx_dedup %>%
  separate(method, c('assembler','clustering','depth', 'quant'), remove = F) %>%
  filter(quant %in% c('map','onts','count','sum'),
         depth %in% c('5m', '10m'),
         # For trinity assembler, keep only trinity clustering
         (assembler == "trinity" & clustering == "trinity") |
           # For non-trinity assemblers, keep only bambu, sim, or corset clustering
           (assembler != "trinity" & clustering %in% c("bambu", "sim", "corset"))) %>%
  mutate(lr_denovo = factor(assembler %in% c('trinity','bambu', 'limma'), 
                            labels = c('lr_denovo', 'others'))) %>% 
  ggplot(aes(x = csum.fp, y = csum.tp, group = method, color = assembler)) + 
  geom_line(aes(linetype = lr_denovo, linewidth = depth)) +
  scale_color_manual(values = cols) +
  scale_linewidth_manual(values = line) +
  ylab('Unique True Positives') +
  xlab('Unique False Postives') + 
  ggtitle('DTU-TX') +
  force_panelsizes(rows = unit(3, "in"),
                   cols = unit(3, "in"))
ggsave('plot/ROC_dtutx_fp_simplified.pdf', width = 7, height = 5)


# DTE-TX
idx <- which(!duplicated(sapply(test3, function(x) dim(x)[1], simplify = T)))

list <- test3[idx]

df_dtetx_all <- lapply(names(list), 
                       function(x) {
                         dt_roc_all(list[[x]], x, truth = true_dte)
                       }) %>%
  Reduce(rbind,. ) 

df_dtetx_dedup <- lapply(names(list), 
                         function(x) {
                           dt_roc_dedup(list[[x]], x, truth = true_dte)
                         }) %>%
  Reduce(rbind,. ) 

df_dtetx_all %>%
  separate(method, c('assembler','clustering','depth', 'quant'), remove = F) %>%
  filter(quant %in% c('map','onts','count','sum'),
         depth %in% c('5m', '10m')) %>% 
  mutate(lr_denovo = factor(assembler %in% c('trinity','bambu', 'limma'), 
                            labels = c('lr_denovo', 'others'))) %>% 
  ggplot(aes(x = rowid, y = csum, group = method, color = assembler)) + 
  geom_line(aes(linetype = lr_denovo, linewidth = depth)) +
  scale_color_manual(values = cols) +
  scale_linewidth_manual(values = line) +
  ylab('Unique True Positives') +
  xlab('Top Ranked Assembled Transcripts') + 
  ggtitle('DTE-TX') +
  force_panelsizes(rows = unit(3, "in"),
                   cols = unit(3, "in"))
ggsave('plot/ROC_dtetx_simplified.pdf', width = 7, height = 5)

df_dtetx_dedup %>%
  separate(method, c('assembler','clustering','depth', 'quant'), remove = F) %>%
  filter(quant %in% c('map','onts','count','sum'),
         depth %in% c('5m', '10m')) %>% 
  mutate(lr_denovo = factor(assembler %in% c('trinity','bambu', 'limma'), 
                            labels = c('lr_denovo', 'others'))) %>% 
  ggplot(aes(x = csum.fp, y = csum.tp, group = method, color = assembler)) + 
  geom_line(aes(linetype = lr_denovo, linewidth = depth)) +
  scale_color_manual(values = cols) +
  scale_linewidth_manual(values = line) +
  ylab('Unique True Positives') +
  xlab('Unique False Postives') + 
  ggtitle('DTE-TX') +
  force_panelsizes(rows = unit(3, "in"),
                   cols = unit(3, "in"))
ggsave('plot/ROC_dtetx_fp_simplified.pdf', width = 7, height = 5)

# DTU-gene
df_dtugene_all <- lapply(names(test2), 
                         function(x) {
                           dg_roc_all(df = test2[[x]], name = x, 
                                      map = filelist[[x]], truth = true_dtugene)
                         }) %>%
  Reduce(rbind,. ) 

df_dtugene_dedup <- lapply(names(test2), 
                           function(x) {
                             dg_roc_dedup(df = test2[[x]], name = x, 
                                          map = filelist[[x]], truth = true_dtugene)
                           }) %>%
  Reduce(rbind,. ) 

df_dtugene_all %>%
  separate(method, c('assembler','clustering','depth', 'quant'), remove = F) %>%
  filter(quant %in% c('map','onts','count','sum'),
         depth %in% c('5m', '10m'),
         # For trinity assembler, keep only trinity clustering
         (assembler == "trinity" & clustering == "trinity") |
           # For non-trinity assemblers, keep only bambu, sim, or corset clustering
           (assembler != "trinity" & clustering %in% c("bambu", "sim", "corset"))) %>% 
  mutate(lr_denovo = factor(assembler %in% c('trinity','bambu', 'limma'), 
                            labels = c('lr_denovo', 'others'))) %>% 
  ggplot(aes(x = rowid, y = csum, group = method, color = assembler)) + 
  geom_line(aes(linetype = lr_denovo, linewidth = depth)) +
  scale_color_manual(values = cols) +
  scale_linewidth_manual(values = line) +
  ylab('Unique True Positives') +
  xlab('Top Ranked Clusters') + 
  ggtitle('DTU-Gene') +
  force_panelsizes(rows = unit(3, "in"),
                   cols = unit(3, "in"))
ggsave('plot/ROC_dtugene_simplified.pdf', width = 7, height = 5)

df_dtugene_dedup %>%
  separate(method, c('assembler','clustering','depth', 'quant'), remove = F) %>%
  filter(quant %in% c('map','onts','count','sum'),
         depth %in% c('5m', '10m'),
         # For trinity assembler, keep only trinity clustering
         (assembler == "trinity" & clustering == "trinity") |
           # For non-trinity assemblers, keep only bambu, sim, or corset clustering
           (assembler != "trinity" & clustering %in% c("bambu", "sim", "corset"))) %>% 
  mutate(lr_denovo = factor(assembler %in% c('trinity','bambu', 'limma'), 
                            labels = c('lr_denovo', 'others'))) %>% 
  ggplot(aes(x = csum.fp, y = csum.tp, group = method, color = assembler)) + 
  geom_line(aes(linetype = lr_denovo, linewidth = depth)) +
  scale_color_manual(values = cols) +
  scale_linewidth_manual(values = line) +
  ylab('Unique True Positives') +
  xlab('Unique False Postives') + 
  ggtitle('DTU-Gene') +
  force_panelsizes(rows = unit(3, "in"),
                   cols = unit(3, "in"))
ggsave('plot/ROC_dtugene_fp_simplified.pdf', width = 7, height = 5)

# DGE
df_dge_all <- lapply(names(test4), 
                     function(x) {
                       dg_roc_all(df = test4[[x]], 
                                  name = x, 
                                  map = filelist[[x]], truth = true_dge)
                     }) %>%
  Reduce(rbind,. ) 

df_dge_dedup <- lapply(names(test4), 
                       function(x) {
                         dg_roc_dedup(df = test4[[x]], 
                                      name = x, 
                                      map = filelist[[x]], truth = true_dge)
                       }) %>%
  Reduce(rbind,. ) 

df_dge_all %>%
  separate(method, c('assembler','clustering','depth', 'quant'), remove = F) %>%
  filter(quant %in% c('map','onts','count','sum'),
         depth %in% c('5m', '10m'),
         # For trinity assembler, keep only trinity clustering
         (assembler == "trinity" & clustering == "trinity") |
           # For non-trinity assemblers, keep only bambu, sim, or corset clustering
           (assembler != "trinity" & clustering %in% c("bambu", "sim", "corset"))) %>% 
  mutate(lr_denovo = factor(assembler %in% c('trinity','bambu', 'limma'), 
                            labels = c('lr_denovo', 'others'))) %>% 
  ggplot(aes(x = rowid, y = csum, group = method, color = assembler)) + 
  geom_line(aes(linetype = lr_denovo, linewidth = depth)) +
  scale_color_manual(values = cols) +
  scale_linewidth_manual(values = line) +
  ylab('Unique True Positives') +
  xlab('Top Ranked Clusters') + 
  ggtitle('DGE-Gene') +
  force_panelsizes(rows = unit(3, "in"),
                   cols = unit(3, "in"))
ggsave('plot/ROC_dge_simplified.pdf', width = 7, height = 5)

df_dge_dedup %>%
  separate(method, c('assembler','clustering','depth', 'quant'), remove = F) %>%
  filter(quant %in% c('map','onts','count','sum'),
         depth %in% c('5m', '10m'),
         # For trinity assembler, keep only trinity clustering
         (assembler == "trinity" & clustering == "trinity") |
           # For non-trinity assemblers, keep only bambu, sim, or corset clustering
           (assembler != "trinity" & clustering %in% c("bambu", "sim", "corset"))) %>% 
  mutate(lr_denovo = factor(assembler %in% c('trinity','bambu', 'limma'), 
                            labels = c('lr_denovo', 'others'))) %>% 
  ggplot(aes(x = csum.fp, y = csum.tp, group = method, color = assembler)) + 
  geom_line(aes(linetype = lr_denovo, linewidth = depth)) +
  scale_color_manual(values = cols) +
  scale_linewidth_manual(values = line) +
  ylab('Unique True Positives') +
  xlab('Unique False Postives') + 
  ggtitle('DGE-Gene') +
  force_panelsizes(rows = unit(3, "in"),
                   cols = unit(3, "in"))
ggsave('plot/ROC_dge_fp_simplified.pdf', width = 7, height = 5)

save.image(file = 'de_data.RData')
