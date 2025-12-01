library(UpSetR)
library(ggpubr)
library(rtracklayer)
library(ggVennDiagram)
library(plyranges)
library(tidyverse)
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

true_dtugene <- readRDS("/vast/projects/davidson_longread/yan.a/simulation_20240501/sqanti_sim/design/true_dtugene.rds")
true_dge <- readRDS("/vast/projects/davidson_longread/yan.a/simulation_20240501/sqanti_sim/design/true_dge.rds")
true_dte <- readRDS("/vast/projects/davidson_longread/yan.a/simulation_20240501/sqanti_sim/design/true_dte.rds")
true_dtutx <- readRDS("/vast/projects/davidson_longread/yan.a/simulation_20240501/sqanti_sim/design/true_dtutx.rds")

# # x is the true set, y is the predicted set
# get_precision_recall <- function(x,y) {
#   tp <- length(intersect(x,y))
#   fn <- length(setdiff(x,y))
#   fp <- length(setdiff(y,x))
#   res <- c(tp, fn, fp)
#   names(res) <- c('tp','fn','fp')
#   return(res)
# }

## evaluate DTU at tx and gene level

files.dtu <- c(list.files('.', pattern = 'dtu_tx.rds', recursive = T),
               list.files('../simulation_1m/bambu/ONT_merged_dge/', pattern = 'dtu_tx.rds', recursive = T, full.names = T))

ntest <- length(files.dtu)

test <- lapply(files.dtu, function(x){
  readRDS(x) 
}) 

#names(test) <- dirname(files.dtu)
names(test) <- c(dirname(files.dtu)[1: (ntest-1)], 
                 'bambu_bambu_count')

pdf('plot/dtu_tx_upset.pdf', width = 8, height = 5)

dtu.list <- list()
dtu.list[[1]] <- true_dtutx

dtu.list[2: (ntest + 1)] <- lapply(1:ntest, function(x){
  method <- dirname(files.dtu)[x]
  txmethod = str_split(method, '_', simplify = T)[1] 
  ids <- test[[x]] %>% 
    rownames_to_column() %>% filter(FDR<0.05) %>%
    mutate(associated_transcript = ifelse(associated_transcript == 'novel', 
                                          paste(txmethod, rowname, sep = '_'),
                                          associated_transcript)) %>% 
    pull(associated_transcript)
}) 

dtu.list[[ntest+2]] <- readRDS('../sqanti_sim/design/limma_dtutx.rds') %>%
  filter(FDR<0.05) %>% pull(transcript_id)

names(dtu.list) <- c('requested', names(test), 'limma_sim_count')

fromList(dtu.list) %>%
  upset(set_size.show = T, keep.order = T, 
        sets = names(dtu.list) ,
        nintersects = 50, order.by = 'freq')
dev.off() # actually overlap can be more between assemblies

files.dtugene <- c(list.files('.', pattern = 'dtu_gene.rds', recursive = T),
                   list.files('../simulation_1m/bambu/ONT_merged_dge/', pattern = 'dtu_gene.rds', recursive = T, full.names = T))

test2 <- lapply(files.dtugene, function(x){
  readRDS(x) 
}) 
names(test2) <- c(dirname(files.dtugene)[1:(ntest-1)], 
                  'bambu_bambu_count')

files <- list.files('.', pattern = 'summary.rds', recursive = T)
files <- files[!grepl('plot', files)]

filelist <- lapply(c(rep(files[2:8], each = 2), 
                     files[9:11], 
                     rep(files[12:13], each = 2), 
                     files[1]), function(x){
  readRDS(x) 
}) # to match the order of DE object
names(filelist) <- names(test2) 

pdf('plot/dtu_gene_upset.pdf', width = 8, height = 5)

dtugene.list <- list()
dtugene.list[[1]] <- true_dtugene

dtugene.list[2:(ntest+1)] <- lapply(1:ntest, function(x){
  if (x %in% 1:(ntest-1)) {
    test2[[x]] %>% 
      rownames_to_column() %>% 
      filter(FDR<0.05) %>%
      left_join(filelist[[x]]$asm_map, by= 'geneid') %>% 
      pull(assigned_gene)
  } else { 
    test2[[x]] %>% rownames_to_column() %>% 
      filter(FDR<0.05) %>%
      mutate(assigned_gene = gene_id) %>%  # for bambu use geneid
      pull(assigned_gene)
  }
})

dtugene.list[[(ntest+2)]] <- readRDS('../sqanti_sim/design/limma_dtugene.rds') %>%
  filter(FDR<0.05) %>% pull(gene_id)

names(dtugene.list) <- c('requested', names(test2) ,'limma_sim_count')

fromList(dtugene.list) %>%
  upset(set_size.show = T, keep.order = T, 
        sets = names(dtugene.list) ,
        nintersects = 50, order.by = 'freq')

dev.off()

## evaluate DTE and DGE

files.dte <- c(list.files('.', pattern = 'dte_tx.rds', recursive = T),
               list.files('../simulation_1m/bambu/ONT_merged_dge/', pattern = 'dte_tx.rds', recursive = T, full.names = T))

test3 <- lapply(files.dte, function(x){
  readRDS(x) 
}) 
names(test3) <- c(dirname(files.dte)[1:(ntest-1)], 
                  'bambu_bambu_count')

pdf('plot/dte_tx_upset.pdf', width = 8, height = 5)

dtetx.list <- list()
dtetx.list[[1]] <- true_dte

dtetx.list[2:(ntest+1)] <- lapply(1:ntest, function(x){
  
  method = dirname(files.dte)[x]
  txmethod = str_split(method, '_', simplify = T)[1] 
  ids <- test3[[x]] %>% 
    rownames_to_column() %>% 
    filter(adj.P.Val < 0.05) %>%
    mutate(associated_transcript = ifelse(associated_transcript == 'novel', 
                                          paste(txmethod, ## different clustering should make no diff
                                                rowname, sep = '_'),
                                          associated_transcript)) %>% 
    pull(associated_transcript)
  
})
dtetx.list[[(ntest+2)]] <- readRDS('../sqanti_sim/design/limma_dtetx.rds') %>%
  filter(adj.P.Val<0.05) %>% pull(transcript_id)

names(dtetx.list) <- c('requested', names(test3) ,'limma_sim_count')

# DTE is independent of clustering, but dependent on quantification
fromList(dtetx.list[!(lengths(dtetx.list) %>% duplicated())]) %>%
  upset(set_size.show = T, keep.order = T, 
        sets = names(dtetx.list[!(lengths(dtetx.list) %>% duplicated())]),
        nintersects = 50, order.by = 'freq')

dev.off()

files.dge <- c(list.files('.', pattern = 'dge_gene.rds', recursive = T),
               list.files('../simulation_1m/bambu/ONT_merged_dge/', pattern = 'dge_gene.rds', recursive = T, full.names = T))

test4 <- lapply(files.dge, function(x){
  readRDS(x)
}) 
names(test4) <- c(dirname(files.dge)[1:ntest-1], 
                  'bambu_bambu_count')

pdf('plot/dge_gene_upset.pdf', width = 8, height = 5)

dge.list <- list()
dge.list[[1]] <- true_dge

dge.list[2:(ntest+1)] <- lapply(1:ntest, function(x){
  if (x %in% (1:ntest-1)) {
    test4[[x]] %>% 
      rownames_to_column() %>% 
      filter(adj.P.Val < 0.05) %>%
      left_join(filelist[[x]]$asm_map, by= c('rowname' = 'geneid')) %>% 
      pull(assigned_gene)
  } else { # for bambu
    test4[[x]] %>% 
      rownames_to_column(var = 'assigned_gene') %>% 
      filter(adj.P.Val < 0.05) %>%
      # left_join(filelist[[x]]$asm_map, by= c('rowname' = 'geneid')) %>% 
      pull(assigned_gene)
  }
})

dge.list[[ntest+2]] <- readRDS('../sqanti_sim/design/limma_dgegene.rds') %>%
  filter(adj.P.Val<0.05) %>% pull(gene_id)

names(dge.list) <- c('requested', names(test4),'limma_sim_count')

fromList(dge.list) %>%
  upset(set_size.show = T, keep.order = T, 
        sets = names(dge.list),
        nintersects = 50, order.by = 'freq')

dev.off()

plot_pr <- function(delist) {
  
  df.pr <- sapply(1:23, function(x) {
    
    get_precision_recall(delist[[1]],
                         delist[[x+1]])
  }, simplify = T)
  
  colnames(df.pr) <- names(delist)[-1]
  
  p1 <- df.pr %>% t %>%
    data.frame() %>%
    rownames_to_column() %>% 
    separate(rowname, c('assembler','clustering','quant')) %>%
    ggscatter(x = 'Precision', y = 'Recall', 
              shape = 'clustering', color = 'assembler', 
              fill = 'quant', #alpha = 'quant',
              size = 5,
              label = 'quant', repel = T) +
    theme(legend.position = "top") + 
    scale_color_manual(values = cols) +
    scale_shape_manual(values = shapes) +
    guides(fill=guide_legend(nrow=2,byrow=TRUE)) #+
  #xlim(c(0,1)) +
  #ylim(c(0,1))
  
  p2 <- df.pr %>% t %>%
    data.frame() %>%
    select(4:6) %>%
    rownames_to_column() %>% 
    separate(rowname, c('assembler', 'clustering', 'quant')) %>%
    # dplyr::select(1:6) %>%
    pivot_longer(4:6) %>% 
    ggdotchart(x='name', y = 'value', #facet.by = 'assembler',
               dot.size = 3, sorting = 'none',
               add = "segment", color = 'assembler', shape = 'clustering',
               rotate = TRUE)  +
    scale_color_manual(values = cols) +
    scale_shape_manual(values = shapes)
  
  p <- cowplot::plot_grid(p1, p2)
  return(p)
}

plot_pr(dtu.list)
ggsave('plot/dtu_tx_pr.pdf', width = 12, height = 5)

plot_pr(dtugene.list)
ggsave('plot/dtu_gene_pr.pdf', width = 12, height = 5)

plot_pr(dtetx.list)
ggsave('plot/dtx_tx_pr.pdf', width = 12, height = 5)

plot_pr(dge.list)
ggsave('plot/dge_gene_pr.pdf', width = 12, height = 5)


## ROC like curve
## DTU-TX
list <- append(test, 
               list(limma_sim_count = readRDS('../sqanti_sim/design/limma_dtutx.rds') %>% 
                      mutate(associated_transcript = transcript_id)))

df_dtutx_all <- lapply(names(list), 
                   function(x) {
                     dt_roc_all(list[[x]], x, truth = true_dtutx)
                   }) %>%
  Reduce(rbind,. ) 

df_dtutx_dedup <- lapply(names(list), 
                   function(x) {
                     dt_roc_dedup(list[[x]], x, truth = true_dtutx)
                   }) %>%
  Reduce(rbind,. ) 

df_dtutx_all %>%
  separate(method, c('assembler','clustering','quant'), remove = F) %>%
  filter(quant %in% c('map','onts','count','sum'),
         # For trinity assembler, keep only trinity clustering
         (assembler == "trinity" & clustering == "trinity") |
           # For non-trinity assemblers, keep only bambu, sim, or corset clustering
           (assembler != "trinity" & clustering %in% c("bambu", "sim", "corset"))) %>%
  mutate(lr_denovo = factor(assembler %in% c('trinity','bambu', 'limma'), 
                            labels = c('lr_denovo', 'others'))) %>% 
  ggplot(aes(x = rowid, y = csum, group = method, color = assembler)) + 
  geom_line(aes(linetype = lr_denovo), linewidth = 1) +
  scale_color_manual(values = cols) +
  ylab('Unique True Positives') +
  xlab('Top Ranked Assembled Transcripts') + 
  ggtitle('DTU-TX') +
  force_panelsizes(rows = unit(3, "in"),
                   cols = unit(3, "in"))
ggsave('plot/ROC_dtutx_simplified.pdf', width = 7, height = 5)

df_dtutx_dedup %>%
  separate(method, c('assembler','clustering','quant'), remove = F) %>%
  filter(quant %in% c('map','onts','count','sum'),
         # For trinity assembler, keep only trinity clustering
         (assembler == "trinity" & clustering == "trinity") |
           # For non-trinity assemblers, keep only bambu, sim, or corset clustering
           (assembler != "trinity" & clustering %in% c("bambu", "sim", "corset"))) %>%
  mutate(lr_denovo = factor(assembler %in% c('trinity','bambu', 'limma'), 
                            labels = c('lr_denovo', 'others'))) %>% 
  ggplot(aes(x = csum.fp, y = csum.tp, group = method, color = assembler)) + 
  geom_line(aes(linetype = lr_denovo), linewidth = 1) +
  scale_color_manual(values = cols) +
  ylab('Unique True Positives') +
  xlab('Unique False Postives') + 
  ggtitle('DTU-TX') +
  force_panelsizes(rows = unit(3, "in"),
                   cols = unit(3, "in"))
ggsave('plot/ROC_dtutx_fp_simplified.pdf', width = 7, height = 5)


# DTE-TX
idx <- which(!duplicated(sapply(test3, function(x) dim(x)[1], simplify = T)))

list <- append(test3[idx], 
               list(limma_sim_count = readRDS('../sqanti_sim/design/limma_dtetx.rds') %>% 
                      mutate(associated_transcript = transcript_id)))

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
  separate(method, c('assembler','clustering','quant'), remove = F) %>%
  filter(quant %in% c('map','onts','count','sum')) %>% 
  mutate(lr_denovo = factor(assembler %in% c('trinity','bambu', 'limma'), 
                            labels = c('lr_denovo', 'others'))) %>% 
  ggplot(aes(x = rowid, y = csum, group = method, color = assembler)) + 
  geom_line(aes(linetype = lr_denovo), linewidth = 1) +
  scale_color_manual(values = cols) +
  ylab('Unique True Positives') +
  xlab('Top Ranked Assembled Transcripts') + 
  ggtitle('DTE-TX') +
  force_panelsizes(rows = unit(3, "in"),
                   cols = unit(3, "in"))
ggsave('plot/ROC_dtetx_simplified.pdf', width = 7, height = 5)

df_dtetx_dedup %>%
  separate(method, c('assembler','clustering','quant'), remove = F) %>%
  filter(quant %in% c('map','onts','count','sum')) %>% 
  mutate(lr_denovo = factor(assembler %in% c('trinity','bambu', 'limma'), 
                            labels = c('lr_denovo', 'others'))) %>% 
  ggplot(aes(x = csum.fp, y = csum.tp, group = method, color = assembler)) + 
  geom_line(aes(linetype = lr_denovo), linewidth = 1) +
  scale_color_manual(values = cols) +
  ylab('Unique True Positives') +
  xlab('Unique False Postives') + 
  ggtitle('DTE-TX') +
  force_panelsizes(rows = unit(3, "in"),
                   cols = unit(3, "in"))
ggsave('plot/ROC_dtetx_fp_simplified.pdf', width = 7, height = 5)

# DTU-gene
list <- append(test2, 
               list(limma_sim_count = readRDS('../sqanti_sim/design/limma_dtugene.rds') ))

df_dtugene_all <- lapply(names(list), 
                       function(x) {
                         dg_roc_all(df = list[[x]], name = x, 
                                    map = filelist[[x]], truth = true_dtugene)
                       }) %>%
  Reduce(rbind,. ) 

df_dtugene_dedup <- lapply(names(list), 
                         function(x) {
                           dg_roc_dedup(df = list[[x]], name = x, 
                                        map = filelist[[x]], truth = true_dtugene)
                         }) %>%
  Reduce(rbind,. ) 

df_dtugene_all %>%
  separate(method, c('assembler','clustering','quant'), remove = F) %>%
  filter(quant %in% c('map','onts','count','sum'),
         # For trinity assembler, keep only trinity clustering
         (assembler == "trinity" & clustering == "trinity") |
           # For non-trinity assemblers, keep only bambu, sim, or corset clustering
           (assembler != "trinity" & clustering %in% c("bambu", "sim", "corset"))) %>% 
  mutate(lr_denovo = factor(assembler %in% c('trinity','bambu', 'limma'), 
                            labels = c('lr_denovo', 'others'))) %>% 
  ggplot(aes(x = rowid, y = csum, group = method, color = assembler)) + 
  geom_line(aes(linetype = lr_denovo), linewidth = 1) +
  scale_color_manual(values = cols) +
  ylab('Unique True Positives') +
  xlab('Top Ranked Clusters') + 
  ggtitle('DTU-Gene') +
  force_panelsizes(rows = unit(3, "in"),
                   cols = unit(3, "in"))
ggsave('plot/ROC_dtugene_simplified.pdf', width = 7, height = 5)

df_dtugene_dedup %>%
  separate(method, c('assembler','clustering','quant'), remove = F) %>%
  filter(quant %in% c('map','onts','count','sum'),
         # For trinity assembler, keep only trinity clustering
         (assembler == "trinity" & clustering == "trinity") |
           # For non-trinity assemblers, keep only bambu, sim, or corset clustering
           (assembler != "trinity" & clustering %in% c("bambu", "sim", "corset"))) %>% 
  mutate(lr_denovo = factor(assembler %in% c('trinity','bambu', 'limma'), 
                            labels = c('lr_denovo', 'others'))) %>% 
  ggplot(aes(x = csum.fp, y = csum.tp, group = method, color = assembler)) + 
  geom_line(aes(linetype = lr_denovo), linewidth = 1) +
  scale_color_manual(values = cols) +
  ylab('Unique True Positives') +
  xlab('Unique False Postives') + 
  ggtitle('DTU-Gene') +
  force_panelsizes(rows = unit(3, "in"),
                   cols = unit(3, "in"))
ggsave('plot/ROC_dtugene_fp_simplified.pdf', width = 7, height = 5)

# DGE
list <- append(test4, 
               list(limma_sim_count = readRDS('../sqanti_sim/design/limma_dgegene.rds') ))

df_dge_all <- lapply(names(list), 
                         function(x) {
                           dg_roc_all(df = list[[x]] %>% rownames_to_column(var = 'geneid'), 
                                      name = x, 
                                      map = filelist[[x]], truth = true_dge)
                         }) %>%
  Reduce(rbind,. ) 

df_dge_dedup <- lapply(names(list), 
                           function(x) {
                             dg_roc_dedup(df = list[[x]] %>% rownames_to_column(var = 'geneid'), 
                                          name = x, 
                                          map = filelist[[x]], truth = true_dge)
                           }) %>%
  Reduce(rbind,. ) 

df_dge_all %>%
  separate(method, c('assembler','clustering','quant'), remove = F) %>%
  filter(quant %in% c('map','onts','count','sum'),
         # For trinity assembler, keep only trinity clustering
         (assembler == "trinity" & clustering == "trinity") |
           # For non-trinity assemblers, keep only bambu, sim, or corset clustering
           (assembler != "trinity" & clustering %in% c("bambu", "sim", "corset"))) %>% 
  mutate(lr_denovo = factor(assembler %in% c('trinity','bambu', 'limma'), 
                            labels = c('lr_denovo', 'others'))) %>% 
  ggplot(aes(x = rowid, y = csum, group = method, color = assembler)) + 
  geom_line(aes(linetype = lr_denovo), linewidth = 1) +
  scale_color_manual(values = cols) +
  ylab('Unique True Positives') +
  xlab('Top Ranked Clusters') + 
  ggtitle('DGE-Gene') +
  force_panelsizes(rows = unit(3, "in"),
                   cols = unit(3, "in"))
ggsave('plot/ROC_dge_simplified.pdf', width = 7, height = 5)

df_dge_dedup %>%
  separate(method, c('assembler','clustering','quant'), remove = F) %>%
  filter(quant %in% c('map','onts','count','sum'),
         # For trinity assembler, keep only trinity clustering
         (assembler == "trinity" & clustering == "trinity") |
           # For non-trinity assemblers, keep only bambu, sim, or corset clustering
           (assembler != "trinity" & clustering %in% c("bambu", "sim", "corset"))) %>% 
  mutate(lr_denovo = factor(assembler %in% c('trinity','bambu', 'limma'), 
                            labels = c('lr_denovo', 'others'))) %>% 
  ggplot(aes(x = csum.fp, y = csum.tp, group = method, color = assembler)) + 
  geom_line(aes(linetype = lr_denovo), linewidth = 1) +
  scale_color_manual(values = cols) +
  ylab('Unique True Positives') +
  xlab('Unique False Postives') + 
  ggtitle('DGE-Gene') +
  force_panelsizes(rows = unit(3, "in"),
                   cols = unit(3, "in"))
ggsave('plot/ROC_dge_fp_simplified.pdf', width = 7, height = 5)

save.image(file = 'de_data.RData')




