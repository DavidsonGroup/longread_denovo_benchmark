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

# load sequin de 
anno <- read.table("~/lab_davidson/yan.a/software/LongReadRNA/sequins/annotations/rnasequin_genes_2.4.tsv",
                   header = TRUE, stringsAsFactors = FALSE) %>%
  mutate(logFC_expected = log(MIX_A / MIX_B),
         isChange = logFC_expected != 0
  )

anno2 <- read.table("~/lab_davidson/yan.a/software/LongReadRNA/sequins/annotations/rnasequin_isoforms_2.4.tsv", 
                    header = TRUE, stringsAsFactors = FALSE) %>%
  mutate(logFC_expected = log(MIX_A / MIX_B),
         isChange = logFC_expected != 0
  )

proportion <- anno2 %>% 
  separate(NAME, c('R','g','t'), '_', remove = F) %>%
  tidyr::unite(gene, c('R','g'), sep = '_') %>%
  left_join(anno, by = c('gene' = 'NAME')) %>%
  group_by(gene) %>% 
  summarise(pct_A = MIX_A.x / sum(MIX_A.x), 
            pct_B = MIX_B.x / sum(MIX_B.x),
            txID = NAME) %>%
  ungroup() %>%
  mutate(FC = pct_A / pct_B,
         isChange = abs(FC-1) > 0.0001)

true_dte <- anno2 %>% filter(isChange) %>% pull(NAME) # 95
true_dge <- anno %>% filter(isChange) %>% pull(NAME) # 44
true_dtutx <- proportion %>% filter(isChange) %>% pull(txID) # 56
true_dtugene <- proportion %>% filter(isChange) %>% pull(gene) %>% unique() #28


## evaluate DTU at tx and gene level
files.dtu <- c(list.files('.', pattern = 'dtu_tx.rds', recursive = T, full.names = T),
               list.files('../bambu/', pattern = 'dtu_tx.rds', recursive = T, full.names = T)
               )

ntest <- length(files.dtu)

test <- lapply(files.dtu, function(x){
  readRDS(x) %>% filter(chrom == 'chrIS')
}) 
# names(test) <- dirname(files.dtu)
names(test) <- c(basename(dirname(files.dtu)[1:(ntest-3)]), 
                 'bambu_bambu_10m_count', 'bambu_bambu_2m_count', 'bambu_bambu_5m_count')

pdf('plot/sequin_dtu_tx_upset.pdf', width = 8, height = 10)

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

# dtu.list[[24]] <- readRDS('../bambu/') %>%
#   filter(FDR<0.05) %>% pull(transcript_id)

names(dtu.list) <- c('sequin', names(test))

fromList(dtu.list) %>%
  upset(set_size.show = T, keep.order = T, 
        sets = names(dtu.list),
        nintersects = 50, order.by = 'freq')
dev.off() # actually overlap can be more between assemblies

files <- list.files('.', pattern = 'summary.rds', recursive = T)
files <- files[grepl('/summary', files)]

filelist <- lapply(c(paste0(str_replace(files.dtu[1:45], "^(.*)_[^/]+/dtu_tx\\.rds$", "\\1"), '/summary.rds'),
                   files[1:3]), function(x){
  readRDS(x) 
})
names(filelist) <- names(test)

files.dtugene <- c(list.files('.', pattern = 'dtu_gene.rds', recursive = T, full.names = T),
                   list.files('../bambu/', pattern = 'dtu_gene.rds', recursive = T, full.names = T)
                   )
test2 <- lapply(1:ntest, function(x){
  
  if (x <= ntest-3) {
    readRDS(files.dtugene[x]) %>% 
      left_join(filelist[[x]]$asm_map, 
                by = 'geneid') %>% 
      filter(str_detect(assigned_gene, '^R'))
  } else {
    readRDS(files.dtugene[x]) %>% 
      left_join(filelist[[x]]$asm_map, 
                by = c('gene_id'='geneid')) %>% 
      filter(str_detect(assigned_gene, '^R'))
  }

}) 
#names(test2) <- dirname(files.dtugene)
names(test2) <- c(basename(dirname(files.dtugene)[1:(ntest-3)]), 
                 'bambu_bambu_10m_count', 'bambu_bambu_2m_count', 'bambu_bambu_5m_count')

# add gene chrom map
# gff <- read_gff('../reference/annotation_gencodev44_sequin.gtf') %>% 
#   filter(type == 'gene') %>%
#   data.frame() %>%
#   select(seqnames, gene_id)

pdf('plot/sequin_dtu_gene_upset.pdf', width = 8, height = 10)

dtugene.list <- list()
dtugene.list[[1]] <- true_dtugene

dtugene.list[2:(ntest+1)] <- lapply(1:ntest, function(x){
  if (x %in% 1:(ntest-3)) {
    test2[[x]] %>% #rownames_to_column() %>% 
      filter(FDR<0.05) %>%
      #left_join(filelist[[x]]$asm_map, by= 'geneid') %>% 
      #left_join(gff, by = c('assigned_gene' = 'gene_id')) %>% ## match and add seqname
      #filter(seqnames == 'chrIS') %>%
      pull(assigned_gene)
  } else {
    test2[[x]] %>% rownames_to_column() %>% 
      filter(FDR<0.05) %>%
      # mutate(assigned_gene = gene_id) %>%  # for bambu use geneid
      # left_join(filelist[[x]]$asm_map, by= c('gene_id'='geneid')) %>% 
      # left_join(gff, by = c('gene_id')) %>% ## match and add seqname
      # filter(seqnames == 'chrIS') %>%
      pull(assigned_gene)
  }
})

# dtugene.list[[12]] <- readRDS('../sqanti_sim/design/limma_dtugene.rds') %>%
#   filter(FDR<0.05) %>% pull(gene_id)

names(dtugene.list) <- c('sequin', names(test2))

fromList(dtugene.list) %>%
  upset(set_size.show = T, keep.order = T, 
        sets = names(dtugene.list),
        nintersects = 50, order.by = 'freq')

dev.off()


## evaluate DTE and DGE

files.dte <- c(list.files('.', pattern = 'dte_tx.rds', recursive = T, full.names = T),
               list.files('../bambu/', pattern = 'dte_tx.rds', recursive = T, full.names = T)
)

test3 <- lapply(files.dte, function(x){
  readRDS(x) %>% filter(chrom == 'chrIS')
}) 
#names(test3) <- dirname(files.dte)
names(test3) <- c(basename(dirname(files.dte)[1:(ntest-3)]), 
                  'bambu_bambu_10m_count', 'bambu_bambu_2m_count', 'bambu_bambu_5m_count')

pdf('plot/sequin_dte_tx_upset.pdf', width = 8, height = 10)

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
                                          paste(txmethod, rowname, ## different clustering should make no diff
                                                sep = '_'),
                                          associated_transcript)) %>% 
    pull(associated_transcript)
  
})
# dtetx.list[[12]] <- readRDS('../sqanti_sim/design/limma_dtetx.rds') %>%
#   filter(adj.P.Val<0.05) %>% pull(transcript_id)

names(dtetx.list) <- c('sequin', names(test3))

# DTE is independent of clustering, but dependent on quantification

fromList(dtetx.list[!(lengths(dtetx.list) %>% duplicated())]) %>%
  upset(set_size.show = T, keep.order = T, 
        sets = names(dtetx.list[!(lengths(dtetx.list) %>% duplicated())]),
        nintersects = 50, order.by = 'freq')

dev.off()

files.dge <- c(list.files('.', pattern = 'dge_gene.rds', recursive = T, full.names = T),
               list.files('../bambu/', pattern = 'dge_gene.rds', recursive = T, full.names = T)
)
test4 <- lapply(1:ntest, function(x){
  readRDS(files.dge[x]) %>% 
    rownames_to_column(var = 'geneid') %>% 
    left_join(filelist[[x]]$asm_map, 
              by = 'geneid') %>% 
    filter(str_detect(assigned_gene, '^R'))
}) 
#names(test4) <- dirname(files.dge)
names(test4) <- c(basename(dirname(files.dge)[1:(ntest-3)]), 
                  'bambu_bambu_10m_count', 'bambu_bambu_2m_count', 'bambu_bambu_5m_count')

pdf('plot/sequin_dge_gene_upset.pdf', width = 8, height = 10)

dge.list <- list()
dge.list[[1]] <- true_dge

dge.list[2:(ntest+1)] <- lapply(1:ntest, function(x){
  if (x %in% 1:(ntest-3)) {
    test4[[x]] %>% # rownames_to_column() %>% 
      filter(adj.P.Val < 0.05) %>%
      # left_join(filelist[[x]]$asm_map, by= c('rowname' = 'geneid')) %>% 
      # left_join(gff, by = c('assigned_gene' = 'gene_id')) %>% ## match and add seqname
      # filter(seqnames == 'chrIS') %>%
      pull(assigned_gene)
  } else {
    test4[[x]] %>% 
      # rownames_to_column(var = 'assigned_gene') %>%
      filter(adj.P.Val < 0.05) %>%
      # mutate(assigned_gene = rowname) %>%  # for bambu use geneid
      # left_join(filelist[[x]]$asm_map, by= c('rowname' = 'geneid')) %>% 
      # left_join(gff, by = c('genes' = 'gene_id')) %>% ## match and add seqname
      # filter(seqnames == 'chrIS') %>%
      pull(assigned_gene)
  }
})

names(dge.list) <-  c('sequin', names(test4))

fromList(dge.list) %>%
  upset(set_size.show = T, keep.order = T, 
        sets = names(dge.list) ,
        nintersects = 50, order.by = 'freq')

dev.off()

plot_pr <- function(delist) {
  
  df.pr <- sapply(1:48, function(x) {
    
    get_precision_recall(delist[[1]],
                         delist[[x+1]])
  }, simplify = T)
  
  colnames(df.pr) <- names(delist)[-1]
  
  p1 <- df.pr %>% t %>%
    data.frame() %>%
    rownames_to_column() %>% 
    separate(rowname, c('assembler','clustering','depth','quant')) %>%
    mutate(depth = factor(depth, levels = c('2m','5m','10m'))) %>%
    filter(quant %in% c('map','onts', 'count','sum')) %>%
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
    filter(quant %in% c('map','onts','count','sum')) %>%
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
ggsave('plot/sequin_dtu_tx_pr.pdf', width = 18, height = 5)

plot_pr(dtugene.list)
ggsave('plot/sequin_dtu_gene_pr.pdf', width = 18, height = 5)

plot_pr(dtetx.list)
ggsave('plot/sequin_dtx_tx_pr.pdf', width = 18, height = 5)

plot_pr(dge.list)
ggsave('plot/sequin_dge_gene_pr.pdf', width = 18, height = 5)


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
ggsave('plot/sequin_ROC_dtutx_simplified.pdf', width = 7, height = 5)

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
ggsave('plot/sequin_ROC_dtutx_fp_simplified.pdf', width = 7, height = 5)


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
ggsave('plot/sequin_ROC_dtetx_simplified.pdf', width = 7, height = 5)

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
ggsave('plot/sequin_ROC_dtetx_fp_simplified.pdf', width = 7, height = 5)

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
ggsave('plot/sequin_ROC_dtugene_simplified.pdf', width = 7, height = 5)

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
ggsave('plot/sequin_ROC_dtugene_fp_simplified.pdf', width = 7, height = 5)

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
ggsave('plot/sequin_ROC_dge_simplified.pdf', width = 7, height = 5)

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
ggsave('plot/sequin_ROC_dge_fp_simplified.pdf', width = 7, height = 5)

save.image(file = 'de_data_sequin.RData')
