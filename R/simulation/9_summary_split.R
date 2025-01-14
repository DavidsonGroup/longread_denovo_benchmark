library(UpSetR)
library(ggpubr)
library(rtracklayer)
library(ggVennDiagram)
library(plyranges)
library(tidyverse)

# set plotting parameters
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cols <- cbPalette[c(1,2,6,4,7,8)]
names(cols) <- c('limma','bambu','isonform','rattle','rnabloom2','trinity')
shapes <- c(15,15,16,17,17,17)
names(shapes) <- c('sim','bambu','corset','isonclust','rattle','trinity')

load("/vast/projects/davidson_longread/yan.a/simulation_20240501/R/de_data.RData")

## Gene level
# 1, DGEs with concordant tx expression changes
true_dge_whole <- readRDS("/vast/projects/davidson_longread/yan.a/simulation_20240501/sqanti_sim/design/dgeid.rds")
# 2, 1-tx caused DTU also DGE
true_dtugene_dtx <- setdiff(readRDS("/vast/projects/davidson_longread/yan.a/simulation_20240501/sqanti_sim/design/true_dtugene.rds"),
                            readRDS("/vast/projects/davidson_longread/yan.a/simulation_20240501/sqanti_sim/design/dtuid.rds")) # or just use comid.rds
# 3, isoform switching only DTUs
true_dtugene_switch <- readRDS("/vast/projects/davidson_longread/yan.a/simulation_20240501/sqanti_sim/design/dtuid.rds")

## Transcript level
# 4, DTU switching, 2000
true_dtutx_switch <- readRDS("/vast/projects/davidson_longread/yan.a/simulation_20240501/sqanti_sim/design/dtutxid.rds")
# 5, dte with single tx changing, 1000
true_dte_single <- readRDS("/vast/projects/davidson_longread/yan.a/simulation_20240501/sqanti_sim/design/combotxid.rds")
# 6, all DTE - single TX changing - switching DTU, 3933
true_dte_dge <- setdiff(readRDS("/vast/projects/davidson_longread/yan.a/simulation_20240501/sqanti_sim/design/true_dte.rds"),
                        c(true_dtutx_switch, true_dte_single))
# 7, all DTU - switching DTU, 3927
true_dtutx_dte <- setdiff(readRDS("/vast/projects/davidson_longread/yan.a/simulation_20240501/sqanti_sim/design/true_dtutx.rds"),
                        c(true_dtutx_switch))

gene1 <- true_dge_whole
gene2 <- true_dtugene_dtx
gene3 <- true_dtugene_switch 

tx1 <- true_dtutx_switch
tx2 <- true_dte_single
tx3 <- true_dte_dge
tx4 <- true_dtutx_dte

## true DGE is gene1 + gene2, 2000
## true DTU-gene gene2 + gene3, 2000
## true DTE is tx1 + tx2 + tx3, 6933
## true DTU-tx is tx1 + tx4, 5927

pdf('plot/ROC_de_summary_split.pdf')
## ROC like curve
## DTU-TX
plot_dtutx <- function(delist = test, true_set = gene1) {
  df_dtutx <- lapply(1:ntest, function(x){
    delist[[x]] %>% 
      rowid_to_column() %>%
      mutate(associated_transcript = ifelse(associated_transcript == 'novel', 
                                            paste(names(delist)[x], rowid, sep = '_'),
                                            associated_transcript),
             is.de = as.numeric(associated_transcript %in% true_set)) %>% 
      select(associated_transcript, is.de) %>%
      unique() %>% 
      rowid_to_column() %>%
      mutate(csum = cumsum(is.de)) %>% 
      select(rowid, csum) %>% 
      mutate(method = names(delist)[x]) %>%
      separate(method, c('assembler','clustering','quant'), remove = F) 
  }) %>% Reduce(rbind,. ) %>%
    rbind(readRDS('../sqanti_sim/design/limma_dtutx.rds') %>% 
            mutate(is.de = as.numeric(transcript_id %in% true_set)) %>% 
            mutate(csum = cumsum(is.de)) %>% 
            rowid_to_column() %>% select(rowid, csum) %>% 
            mutate(assembler = 'limma', clustering = 'sim', quant = 'count', method = 'limma') ) 
  
  p <- df_dtutx %>%
    filter(quant %in% c('map','onts','count')) %>%
    ggplot(aes(x = rowid, y = csum, group = method, color = assembler)) + 
    geom_line(aes(linetype = clustering == 'corset' #, alpha = quant
    )) +
    scale_color_manual(values = cols) +
    ylab('Unique True Positives') +
    xlab('Top Ranked Assembled Transcripts') + 
    ggtitle('DTU-TX')
  return(p)
}

plot_dtutx(delist = test, true_set = tx1) + 
  ggtitle('DTU-Tx due to switching') # DTU due to isoform switching
plot_dtutx(delist = test, true_set = tx4) +
  ggtitle('DTU-Tx due to 1 tx from the gene been changed') # DTU due to 1-tx been changed
plot_dtutx(delist = test, true_set = c(tx1,tx4))

# ggsave('plot/ROC_dtutx.pdf', width = 8, height = 5)

# DTE-TX
plot_dtetx <- function(delist = test3, true_set = gene1) {
  idx <- which(!duplicated(sapply(delist, function(x) dim(x)[1], simplify = T)))
  df_dtetx <- lapply(idx, function(x){
    delist[[x]] %>% 
      rowid_to_column() %>%
      mutate(associated_transcript = ifelse(associated_transcript == 'novel', 
                                            paste(names(delist)[x], rowid, sep = '_'),
                                            associated_transcript),
             is.de = as.numeric(associated_transcript %in% true_set)) %>% 
      select(associated_transcript, is.de) %>%
      unique() %>% 
      rowid_to_column() %>%
      mutate(csum = cumsum(is.de)) %>% 
      select(rowid, csum) %>% 
      mutate(method = names(delist)[x]) %>%
      separate(method, c('assembler','clustering','quant'), remove = F) 
  }) %>% Reduce(rbind,. ) %>%
    rbind(readRDS('../sqanti_sim/design/limma_dtetx.rds') %>% 
            mutate(is.de = as.numeric(transcript_id %in% true_set)) %>% 
            mutate(csum = cumsum(is.de)) %>% 
            rowid_to_column() %>% select(rowid, csum) %>% 
            mutate(assembler = 'limma', clustering = 'sim', quant = 'count', method = 'limma') ) 
  
  p <- df_dtetx %>%
    filter(quant %in% c('map','onts','count')) %>%
    ggplot(aes(x = rowid, y = csum, group = method, color = assembler)) + 
    geom_line() +
    scale_color_manual(values = cols) +
    ylab('Unique True Positives') +
    xlab('Top Ranked Assembled Transcripts') + 
    ggtitle('DTE-Tx')
  return(p)
}

plot_dtetx(delist = test3, true_set = tx1) +
  ggtitle('DTE-Tx due to switching') # DTE due to isoform switching
plot_dtetx(delist = test3, true_set = tx2) + 
  ggtitle('DTE-Tx due to 1-tx been changed') # DTE due to 1-tx been changed
plot_dtetx(delist = test3, true_set = tx3)+ 
  ggtitle('DTE-Tx due to gene been changed') # DTE due to gene been changed

plot_dtetx(delist = test3, true_set = c(tx1,tx2,tx3))

# ggsave('plot/ROC_dtetx.pdf', width = 8, height = 5)


# DTU-gene
plot_dtugene <- function(delist = test2, mapping = filelist, true_set = gene1){
  df_dtugene <- lapply(1:ntest, function(x){
    if (x %in% 1:(ntest-1)) {
      delist[[x]] %>% 
        # rowid_to_column() %>%
        left_join(mapping[[x]]$asm_map, by= 'geneid') %>% 
        mutate(is.de = assigned_gene %in% true_set) %>%
        select(assigned_gene, is.de) %>%
        unique() %>% 
        rowid_to_column() %>%
        mutate(csum = cumsum(is.de)) %>% 
        select(rowid, csum) %>% 
        mutate(method = names(delist)[x]) %>%
        separate(method, c('assembler','clustering','quant'), remove = F) 
    } else {
      delist[[x]] %>% 
        # rowid_to_column() %>%
        # left_join(filelist[[x]]$asm_map, by= c('gene_id'='geneid')) %>% 
        mutate(assigned_gene = gene_id) %>%
        mutate(is.de = assigned_gene %in% true_set) %>%
        select(assigned_gene, is.de) %>%
        unique() %>% 
        rowid_to_column() %>%
        mutate(csum = cumsum(is.de)) %>% 
        select(rowid, csum) %>% 
        mutate(method = names(delist)[x]) %>%
        separate(method, c('assembler','clustering','quant'), remove = F) 
    }
  }) %>% Reduce(rbind,. ) %>%
    rbind(readRDS('../sqanti_sim/design/limma_dtugene.rds') %>% 
            mutate(is.de = as.numeric(gene_id %in% true_set)) %>% 
            mutate(csum = cumsum(is.de)) %>% 
            rowid_to_column() %>% select(rowid, csum) %>% 
            mutate(assembler = 'limma', clustering = 'sim', quant = 'count', method = 'limma') ) 
  
  p <- df_dtugene %>%
    filter(quant %in% c('map','onts','count')) %>%
    ggplot(aes(x = rowid, y = csum, group = method, color = assembler)) + 
    geom_line(aes(linetype = clustering == 'corset' #, alpha = quant
    )) +
    scale_color_manual(values = cols) +
    ylab('Unique True Positives') +
    xlab('Top Ranked Clusters') + 
    ggtitle('DTU-Gene')
  return(p)
}
# ggsave('plot/ROC_dtugene.pdf', width = 8, height = 5)

plot_dtugene(delist = test2, true_set = gene3) + 
  ggtitle('DTU-gene due to switching') # DTU-gene due to isoform switching
plot_dtugene(delist = test2, true_set = gene2) + 
  ggtitle('DTU-gene due to 1 tx from the gene been changed') # DTU-gene due to 1-tx been changed
plot_dtugene(delist = test2, true_set = c(gene2, gene3))

# DGE-gene
plot_dgegene <- function(delist = test4, mapping = filelist, true_set = gene1){
  df_dge <- lapply(1:ntest, function(x){
    if (x %in% 1:(ntest-1)) {
      delist[[x]] %>% 
        rownames_to_column() %>%
        # rowid_to_column() %>%
        left_join(filelist[[x]]$asm_map, by= c('rowname' = 'geneid')) %>% 
        mutate(is.de = assigned_gene %in% true_set) %>%
        select(assigned_gene, is.de) %>%
        unique() %>% 
        rowid_to_column() %>%
        mutate(csum = cumsum(is.de)) %>% 
        select(rowid, csum) %>% 
        mutate(method = names(delist)[x]) %>%
        separate(method, c('assembler','clustering','quant'), remove = F) 
    } else {
      delist[[x]] %>% 
        rownames_to_column(var = 'assigned_gene') %>%
        # rowid_to_column() %>%
        # left_join(filelist[[x]]$asm_map, by= c('rowname' = 'geneid')) %>% 
        mutate(is.de = assigned_gene %in% true_set) %>%
        select(assigned_gene, is.de) %>%
        unique() %>% 
        rowid_to_column() %>%
        mutate(csum = cumsum(is.de)) %>% 
        select(rowid, csum) %>% 
        mutate(method = names(delist)[x]) %>%
        separate(method, c('assembler','clustering','quant'), remove = F) 
    } 
  }) %>% Reduce(rbind,. ) %>%
    rbind(readRDS('../sqanti_sim/design/limma_dgegene.rds') %>% 
            mutate(is.de = as.numeric(gene_id %in% true_set)) %>% 
            mutate(csum = cumsum(is.de)) %>% 
            rowid_to_column() %>% select(rowid, csum) %>% 
            mutate(assembler = 'limma', clustering = 'sim', quant = 'count', method = 'limma') ) 
  
  p <- df_dge %>%
    filter(quant %in% c('map','onts','count')) %>%
    ggplot(aes(x = rowid, y = csum, group = method, color = assembler)) + 
    geom_line(aes(linetype = clustering == 'corset' #, alpha = quant)
    )) +
    scale_color_manual(values = cols) +
    ylab('Unique True Positives') +
    xlab('Top Ranked Clusters') + 
    ggtitle('DGE-Gene')
  return(p)
}
plot_dgegene(delist = test4, true_set = gene1) +
  ggtitle('DGE-gene due to gene change')# DGE-gene due to gene change
plot_dgegene(delist = test4, true_set = gene2) +
  ggtitle('DGE-gene due to 1-tx change') # DTU-gene due to 1-tx been changed
plot_dgegene(delist = test4, true_set = c(gene2, gene1))

dev.off()
# ggsave('plot/ROC_dgegene.pdf', width = 8, height = 5)

# save.image(file = 'de_data.RData')

