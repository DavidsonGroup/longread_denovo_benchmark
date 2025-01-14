library(data.table)
library(ggrepel)
library(scales)
library(tidyverse)
library(ggpubr)
library(fastmatch)
library(UpSetR)
source("~/lab_davidson/yan.a/software/scripts_denovo/R/upset_sqanti.R")
source("~/lab_davidson/yan.a/software/scripts_denovo/R/get_union_f1.R")
source("~/lab_davidson/yan.a/software/scripts_denovo/R/get_cluster_exp.R")
source("~/lab_davidson/yan.a/software/scripts_denovo/R/calculate_denovo_cluster.R")

# set plotting parameters
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cols <- cbPalette[c(1,2,6,4,7,8)]
names(cols) <- c('ref','bambu','isonform','rattle','rnabloom2','trinity')
shapes <- c(15,15,16,17,17,17)
names(shapes) <- c('ref','bambu','corset','isonclust','rattle','trinity')

files <- list.files('.', pattern = 'summary.rds', recursive = T) 
files <- files[grepl('m/', files)]

filelist <- lapply(files, function(x){
  readRDS(x)
})
names(filelist) <- dirname(files)

lapply(names(filelist), function(x){
  filelist[[x]]$pure_denovo %>% 
    mutate(num_bin = cut(Unique_Count, breaks = c(0,1,2,5,10,Inf))) %>%
    group_by(num_bin) %>%
    summarise(n = sum(n)) %>%
    mutate(method = x)
}) %>% rbindlist() %>%
  separate(method, c('assembler','clustering','depth'), remove = F) %>%
  mutate(depth = factor(depth, levels = c('2m','5m','10m')),
         clustering = ifelse(clustering == 'corset', 'corset', 'native')) %>%
  ggplot(aes(x = depth, y = n, fill = num_bin)) + 
  geom_bar(stat = 'identity', position = 'stack') +
  geom_text(aes(label = n), colour = "white", size = 3.5,
            position=position_stack(vjust=0.5)) + 
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_grid(clustering~assembler) + 
  ylab('Number of de novo clusters') +
  guides(fill=guide_legend(nrow=2, byrow=TRUE)) 
ggsave('plot/pct_denovocluster_with_genes.pdf', width = 8, height = 8)

lapply(names(filelist), function(x){
  filelist[[x]]$pure_known %>% 
    mutate(num_bin = cut(Unique_Count, breaks = c(0,1,2,5,10,Inf))) %>%
    group_by(num_bin) %>%
    summarise(n = sum(n)) %>%
    mutate(method = x)
}) %>% rbindlist() %>%
  separate(method, c('assembler','clustering','depth'), remove = F) %>%
  mutate(depth = factor(depth, levels = c('2m','5m','10m')),
         clustering = ifelse(clustering == 'corset', 'corset', 'native')) %>%
  ggplot(aes(x = depth, y = n, fill = num_bin)) + 
  geom_bar(stat = 'identity', position = 'stack') +
  geom_text(aes(label = n), colour = "white", size = 3.5,
            position=position_stack(vjust=0.5)) + 
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_grid(clustering~assembler) + 
  ylab('Number of genes') +
  guides(fill=guide_legend(nrow=2, byrow=TRUE)) 
ggsave('plot/pct_genes_with_denovoclusters.pdf', width = 8, height = 8)

## would not change even if we changed the minimap2 paramters, because previously novel transcript (without txid) are still in the correct and same gene
# note the evaluation set is not identical here, and can not use intersection to get a shared evaluation set
lapply(filelist, function(x) {
  data.frame(precision = x$Precision, 
             recall = x$Recall,
             Num_eval = x$known_denovotx)
}) %>% rbindlist() %>%
  mutate(method = dirname(files)) %>%
  separate(method, c('assembler', 'clustering', 'depth'), remove = F) %>%
  mutate(depth = factor(depth, levels = c('2m','5m','10m'))) %>%
  ggplot(aes(x = precision, y = recall, color = assembler, shape = clustering, label = clustering)) + 
  geom_point(aes(size = depth)) +  # geom_point(aes(size = Num_eval)) 
  geom_text_repel() +
  scale_colour_manual(values = cols) +  
  scale_shape_manual(values = shapes) +
  xlim(c(0,1)) +
  ylim(c(0,1)) 
ggsave('plot/cluster_precision_recall.pdf', width = 6, height = 5)


lapply(filelist, function(x) {
  data.frame(precision = x$Precision, 
             recall = x$Recall,
             Num_eval = x$known_denovotx)
}) %>% rbindlist() %>%
  mutate(method = dirname(files)) %>%
  separate(method, c('assembler', 'clustering', 'depth'), remove = F) %>%
  mutate(depth = factor(depth, levels = c('2m','5m','10m'))) %>%
  ggplot(aes(x = precision, y = recall, color = assembler, shape = clustering, label = clustering)) + 
  geom_point(aes(size = Num_eval)) + 
  geom_text_repel() +
  scale_colour_manual(values = cols) +
  scale_shape_manual(values = shapes) +
  facet_grid(assembler~depth) + 
  xlim(c(0,1)) +
  ylim(c(0,1)) 
ggsave('plot/cluster_precision_recall_facet.pdf', width = 10, height = 9)


# ignore novel assigned gene
# lapply(filelist, function(x) {
#   data.frame(one2one = (x$asm_map %>% filter(!duplicated(assigned_gene) & !grepl('novel', assigned_gene)) %>% nrow()) / nrow(x$asm_map),
#              total_cluster = x$total_denovocluster,
#              one2one_num = nrow((x$asm_map %>% filter(!duplicated(assigned_gene) & !grepl('novel', assigned_gene)))))
# }) %>% rbindlist() %>%
#   mutate(method = dirname(files)) %>%
#   separate(method, c('assembler', 'clustering'), remove = F) %>%
#   ggplot(aes(x = total_cluster, y = one2one, color = clustering, shape = assembler)) + 
#   geom_point(size = 3) + 
#   ylim(c(0,1)) 

lapply(filelist, function(x) {
  data.frame(one2one = (x$asm_map %>% filter(!duplicated(assigned_gene) & !grepl('novel', assigned_gene)) %>% nrow()) / nrow(x$asm_map),
             total_cluster = x$total_denovocluster,
             one2one_num = nrow((x$asm_map %>% filter(!duplicated(assigned_gene) & !grepl('novel', assigned_gene)))))
}) %>% rbindlist() %>%
  mutate(method = dirname(files)) %>%
  separate(method, c('assembler', 'clustering', 'depth'), remove = F) %>%
  mutate(depth = factor(depth, levels = c('2m','5m','10m'))) %>%
  ggplot(aes(x = total_cluster, y = one2one_num, color = assembler, shape = clustering, label = depth)) + 
  geom_point(size = 5) + 
  geom_text_repel() +
  scale_colour_manual(values = cols) +
  scale_shape_manual(values = shapes) +
  geom_abline(slope = 1)
ggsave('plot/cluster_number_unique_cluster_gene.pdf', width = 6, height = 5)

lapply(filelist, function(x) {
  data.frame(pure_denovo = x$pure_denovo$n[1]/(sum(x$pure_denovo$n)),
             pure_known = x$pure_known$n[1]/(sum(x$pure_known$n)))
}) %>% rbindlist() %>%
  mutate(method = dirname(files)) %>%
  separate(method, c('assembler', 'clustering', 'depth'), remove = F) %>%
  mutate(depth = factor(depth, levels = c('2m','5m','10m'))) %>%
  ggplot(aes(x = pure_known, y = pure_denovo, color = assembler, shape = clustering, label = depth)) + 
  geom_point(size = 5) + 
  geom_text_repel() + 
  scale_colour_manual(values = cols) +
  scale_shape_manual(values = shapes) +
  ylim(c(0,1)) +
  xlim(c(0,1)) 
ggsave('plot/cluster_pct_uniq_cluster_gene.pdf', width = 6, height = 5)

lapply(filelist, function(x) {
  data.frame(precision = x$Precision, 
             recall = x$Recall,
             F1 = x$F1,
             ARI= x$eval["adj_rand_index"],
             homogeneity = x$eval["homogeneity"], 
             completeness = x$eval["completeness"])
}) %>% rbindlist() %>%
  mutate(method = dirname(files)) %>%
  separate(method, c('assembler', 'clustering', 'depth'), remove = F) %>%
  mutate(depth = factor(depth, levels = c('2m','5m','10m')),
         clustering = ifelse(clustering =='corset', 'corset', 'native')) %>% 
  pivot_longer(1:6) %>% 
  ggdotchart(x='name', y = 'value', facet.by = 'depth',
             dot.size = 3, sorting = 'none',
             add = "segment", color = 'assembler', shape = 'clustering',
             rotate = TRUE)  +
  scale_color_manual(values = cols)
ggsave('plot/cluster_all_metrics.pdf', width = 10, height = 5)

lapply(filelist, function(x) {
  data.frame(precision = x$Precision, 
             recall = x$Recall,
             F1 = x$F1)
}) %>% rbindlist() %>%
  mutate(method = dirname(files)) %>%
  separate(method, c('assembler', 'clustering', 'depth'), remove = F) %>%
  mutate(depth = factor(depth, levels = c('2m','5m','10m')),
         clustering = ifelse(clustering =='corset', 'corset', 'native')) %>% 
  # dplyr::select(1:6) %>%
  pivot_longer(1:3) %>% 
  ggplot(aes(x=clustering, y=value, group = assembler, colour = assembler)) + 
  geom_point() +
  geom_line() + 
  facet_grid(name~depth) +
  ylim(c(0,1)) + 
  scale_color_manual(values = cols)
ggsave('plot/cluster_f1_metrics.pdf', width = 8, height = 5)

lapply(filelist, function(x) {
  data.frame(pure_denovo = x$pure_denovo$n[1]/(sum(x$pure_denovo$n)),
             pure_known = x$pure_known$n[1]/(sum(x$pure_known$n))) %>% 
    cbind(c(x$eval, x[c(2:5,7,11:14)]) %>% unlist() %>% data.frame() %>% t())
}) %>% rbindlist() %>%
  mutate(method = dirname(files)) %>%
  separate(method, c('assembler', 'clustering', 'depth'), remove = F) %>%
  mutate(depth = factor(depth, levels = c('2m','5m','10m')),
         clustering = ifelse(clustering =='corset', 'corset', 'native')) %>% 
  write.csv('plot/cluster_all_metrics.csv')

## get cluster level gene expression

# remove bambu clsuter results from list
# filelist_denovo <- filelist[-c(1:3)]

salmon_count <- readRDS('plot/salmon_count.rds')
requested_counts <- readRDS('plot/requested_counts.rds')

# map the quant and clustering
df1 <- data.frame(cluster_file = names(filelist)) %>%
  separate(cluster_file, c('assembler','clustering','depth'), remove = F) %>%
  unite(dataset, 'assembler','depth', remove = F)

df2 <- data.frame(quant_file = names(salmon_count)) %>%
  separate(quant_file, c(NA, 'depth', NA, 'assembler', 'quant'), remove = F) %>%
  filter(quant %in% c('onts','map') | assembler == 'bambu') %>%
  unite(dataset, 'assembler','depth')
df2$dataset[12] <- 'trinity_10m'

df1 <- df1 %>%
  left_join(df2, 'dataset') %>%
  mutate(quant_file2 = str_replace(quant_file, 'onts', 'ontp')) %>%
  mutate(quant_file2 = str_replace(quant_file2, 'map', 'align')) 

# run get_cluster_exp in multiple combinations

# add up counts by cluster (secondary vs primary for ont, in the case of trinity, it will be map vs align, map is mostly superior) 
# and then use asm_map to assign cluster to gene
# if multiple cluster match to same gene, count cluster by adding up or taking only the max

# for no sequin
all <- lapply(1:nrow(df1), function(x) {
  res1 <- get_cluster_exp(salmon_count[[df1$quant_file[x]]], 
                          filelist[[df1$cluster_file[x]]],
                          requested_counts = requested_counts %>% filter(!str_detect(transcript_id, '^R')))
  
  res2 <- get_cluster_exp(salmon_count[[df1$quant_file2[x]]], 
                          filelist[[df1$cluster_file[x]]],
                          requested_counts = requested_counts %>% filter(!str_detect(transcript_id, '^R')))
  
  list(sec=res1, pri=res2)
})
names(all) <- df1$cluster_file

cor.sec.sum <- sapply(all, function(x) {
  cor(x$sec[, c('logbambu', 'logassigned_sum')])[1,2]
}, simplify = T) 
cor.sec.max <- sapply(all, function(x) {
  cor(x$sec[, c('logbambu', 'logassigned_max')])[1,2]
}, simplify = T) 
cor.pri.sum <- sapply(all, function(x) {
  cor(x$pri[, c('logbambu', 'logassigned_sum')])[1,2]
}, simplify = T) 
cor.pri.max <- sapply(all, function(x) {
  cor(x$pri[, c('logbambu', 'logassigned_max')])[1,2]
}, simplify = T) 

cor <- data.frame(cluster_file = names(all),
                  cor.sec.sum = cor.sec.sum, cor.sec.max = cor.sec.max, 
                  cor.pri.sum = cor.pri.sum, cor.pri.max = cor.pri.max) %>%
  separate(cluster_file, c('assembler','clustering','depth'), remove = F) %>%
  unite(dataset, 'assembler','depth', remove = F)
write.csv(cor, 'plot/cor_cluster.csv')

cor %>%
  pivot_longer(6:9, values_to = 'cor') %>% 
  filter(name == 'cor.sec.max') %>%
  mutate(corset = ifelse(clustering == 'corset', 'corset', 'native'),
         depth = factor(depth, levels = c('2m','5m', '10m'))) %>%
  ggplot(aes(x = assembler, y = cor, colour = assembler)) +
  geom_point(aes(shape = depth, size = 1)) + 
  facet_grid(corset~.) +
  scale_color_manual(values = cols) +
  ylim(c(0,1)) +
  ylab('Correlation (cluster vs Ref gene)')
ggsave('plot/cluster_exp_cor.pdf', width = 10, height = 4)

# all$isonform_corset_2m$sec %>%
#   ggscatter(x='logbambu', 'logassigned_sum') +
#   geom_abline()

cor.corset <- cor %>% filter(clustering == 'corset')
cor.native <- cor %>% filter(clustering != 'corset')

cor1 <- full_join(cor.corset, cor.native,
                  by = c('dataset', 'assembler', 'depth'), 
                  suffix = c('.corset', '.native'))
  
cor2 <- read.csv('plot/cor.csv')
cor2 <- left_join(cor2 %>% 
                unite(dataset, 'assembler','depth', remove = F) %>%
                filter(quant %in% c('onts', 'map') | assembler == 'bambu') %>%
                select(c(2,3,4,6)),
              cor2 %>%   
                unite(dataset, 'assembler','depth', remove = F) %>%
                filter(quant %in% c('ontp', 'align') | assembler == 'bambu') %>% 
                select(c(2,3,4,6)), by = 'dataset',
              suffix = c('.sec', '.pri')) 

cor.final <- left_join(cor1, cor2, 
                       by = 'dataset')
write.csv(cor.final, 'plot/cor_finalcombined.csv')

cor.final %>%
  select(c(2,3,5,6:9,12:17,19:20)) %>% # 12 columns of correlation
  pivot_longer(-c(1:3)) %>% 
  filter(str_detect(name, 'sec') & !str_detect(name, 'sum') & !str_detect(name, 'tx')) %>%
  mutate(method = ifelse(name == 'cor.sec.max.corset', 'corset_cluster',
                         ifelse(name == 'cor.sec.max.native', 'native_cluster',
                                ifelse(name == 'gene_cor.sec', 'sqanti3_gene', 'sqanti3_transcript'))),
         depth = factor(depth, levels = c('2m','5m', '10m'))) %>%
  ggplot(aes(x = assembler, y = value, 
             colour = assembler)) +
  geom_point(aes(shape = method), size = 3) + 
  facet_grid(~depth, scales = 'free_x') +
  scale_color_manual(values = cols) +
  ylim(c(0,1)) +
  ylab('Correlation (De novo vs Bambu)')
ggsave('plot/cor_exp_gene.pdf', width = 10, height = 3)

cor.final %>%
  select(c(2,3,5,6:9,12:17,19:20)) %>% # 12 columns of correlation
  pivot_longer(-c(1:3)) %>% 
  filter(str_detect(name, 'tx_cor.sec')) %>%
  mutate(method = ifelse(name == 'cor.sec.max.corset', 'corset_cluster',
                         ifelse(name == 'cor.sec.max.native', 'native_cluster',
                                ifelse(name == 'gene_cor.sec', 'sqanti3_gene', 'sqanti3_transcript'))),
         depth = factor(depth, levels = c('2m','5m', '10m'))) %>%
  ggplot(aes(x = assembler, y = value, 
             colour = assembler)) +
  geom_point(aes(shape = method), size = 3) + 
  facet_grid(~depth, scales = 'free_x') +
  scale_color_manual(values = cols) +
  ylim(c(0,1)) +
  ylab('Correlation (De novo vs Bambu)')
ggsave('plot/cor_exp_tx.pdf', width = 10, height = 3)

# for sequin cluster correlation

all.sequin <- lapply(1:nrow(df1), function(x) {
  res1 <- get_cluster_exp(salmon_count[[df1$quant_file[x]]], 
                          filelist[[df1$cluster_file[x]]],
                          requested_counts = requested_counts %>% filter(str_detect(transcript_id, '^R')))
  
  res2 <- get_cluster_exp(salmon_count[[df1$quant_file2[x]]], 
                          filelist[[df1$cluster_file[x]]],
                          requested_counts = requested_counts %>% filter(str_detect(transcript_id, '^R')))
  
  list(sec=res1, pri=res2)
})
names(all.sequin) <- df1$cluster_file

cor.sec.sum <- sapply(all.sequin, function(x) {
  cor(x$sec[, c('logbambu', 'logassigned_sum')])[1,2]
}, simplify = T) 
cor.sec.max <- sapply(all.sequin, function(x) {
  cor(x$sec[, c('logbambu', 'logassigned_max')])[1,2]
}, simplify = T) 
cor.pri.sum <- sapply(all.sequin, function(x) {
  cor(x$pri[, c('logbambu', 'logassigned_sum')])[1,2]
}, simplify = T) 
cor.pri.max <- sapply(all.sequin, function(x) {
  cor(x$pri[, c('logbambu', 'logassigned_max')])[1,2]
}, simplify = T) 

cor.sequin <- data.frame(cluster_file = names(all),
                  cor.sec.sum = cor.sec.sum, cor.sec.max = cor.sec.max, 
                  cor.pri.sum = cor.pri.sum, cor.pri.max = cor.pri.max) %>%
  separate(cluster_file, c('assembler','clustering','depth'), remove = F) %>%
  unite(dataset, 'assembler','depth', remove = F)

write.csv(cor.sequin, 'plot/cor_cluster_sequin.csv')

cor.corset <- cor.sequin %>% filter(clustering == 'corset')
cor.native <- cor.sequin %>% filter(clustering != 'corset')

cor1 <- full_join(cor.corset, cor.native,
                  by = c('dataset', 'assembler', 'depth'), 
                  suffix = c('.corset', '.native'))

cor2 <- readRDS('plot/sequin_cor.rds')

cor2 <- data.frame(gene_cor_a = sapply(cor2[grep('mix_a', names(cor2))], function(x) x$cor.gene['logA'], simplify = T),
           gene_cor_b = sapply(cor2[grep('mix_b', names(cor2))], function(x) x$cor.gene['logB'], simplify = T),
           tx_cor_a = sapply(cor2[grep('mix_a', names(cor2))], function(x) x$cor.tx['logA'], simplify = T),
           tx_cor_b = sapply(cor2[grep('mix_b', names(cor2))], function(x) x$cor.tx['logB'], simplify = T)) %>%
  rownames_to_column() %>% 
  separate(rowname, c('quant', 'mix', 'depth', 'assembler', 'log'), sep = '\\.')

cor2 <- left_join(cor2 %>% 
                    unite(dataset, 'assembler','depth', remove = F) %>%
                    filter(quant %in% c('onts', 'map') | assembler == 'bambu') %>%
                    select(-c(2,6)),
                  cor2 %>%   
                    unite(dataset, 'assembler','depth', remove = F) %>%
                    filter(quant %in% c('ontp', 'align') | assembler == 'bambu') %>% 
                    select(-c(2,6)), by = 'dataset',
                  suffix = c('.sec', '.pri')) 

cor.final.sequin <- left_join(cor1, cor2, 
                       by = 'dataset')
write.csv(cor.final.sequin, 'plot/cor_finalcombined_sequin.csv')


## calculating f1 on shared tx onnly

sqanti_summary <- readRDS("plot/sqanti_summary.rds")

# no sequin
tx_union <- lapply(sqanti_summary, function(x){
  x[, c('associated_transcript', 'associated_gene')]
}) %>% rbindlist() %>%
  filter(!str_detect(associated_transcript, '^novel') & !is.na(associated_transcript) & !str_detect(associated_gene, '_') & !str_detect(associated_transcript, '^R')) %>%
  unique()

feature_count <- upset_sqanti(sqanti_summary[c(1,5,6,9,12)], featuretype = 'transcript', 
                              subcat = 'all', names = names(sqanti_summary)[c(1,5,6,9,12)])
ncol <- length(sqanti_summary[c(1,5,6,9,12)])

subset_gt4 <- tx_union %>% 
  filter(associated_transcript %in% rownames(feature_count)[rowSums(feature_count)>=(ncol-1)])
subset_shared <- tx_union %>% 
  filter(associated_transcript %in% rownames(feature_count)[rowSums(feature_count)==ncol])

lapply(names(filelist), function(x) {
  print(x)
  subset_known <- filelist[[x]]$sqanti_anno %>%
    mutate(gene_category = ifelse(structural_category %in% c('antisense','fusion','intergenic','genic_intron', NA), 
                                  structural_category,
                                  'known_gene')) %>%
    mutate(gene_category = ifelse(is.na(gene_category), 'not_mapped', gene_category)) %>%
    dplyr::filter(gene_category=='known_gene' ) %>% ## any de novo transcript mapped to known gene (FSM/ISM, NNC/NIC, genic)
    select(associated_transcript, associated_gene)
  
  df <- rbind(# get_union_f1(sqanti_anno = filelist[[x]]$sqanti_anno, truth = tx_union),
    get_union_f1(sqanti_anno = filelist[[x]]$sqanti_anno, truth = subset_gt4, remove.nomatch = F),
    get_union_f1(sqanti_anno = filelist[[x]]$sqanti_anno, truth = subset_shared, remove.nomatch = F),
    get_union_f1(sqanti_anno = filelist[[x]]$sqanti_anno, truth = subset_gt4),
    get_union_f1(sqanti_anno = filelist[[x]]$sqanti_anno, truth = subset_shared),
    get_union_f1(sqanti_anno = filelist[[x]]$sqanti_anno, truth = subset_known, keep.all.truth = F)) %>%
    mutate(clustering = x)
  df$method = c('all_but1_w_nomatch', 'shared_w_nomatch', 
                'all_but1_wo_nomatch', 'shared_wo_nomatch', 'de_novo')
  df
}) %>% rbindlist() %>% 
  write.csv('plot/new_cluster_metric.csv')

# sequin
tx_union <- lapply(sqanti_summary, function(x){
  x[, c('associated_transcript', 'associated_gene')]
}) %>% rbindlist() %>%
  filter(str_detect(associated_transcript, '^R')) %>%
  unique()

feature_count <- upset_sqanti(sqanti_summary[c(1,5,6,9,12)], featuretype = 'transcript', 
                              subcat = 'all', names = names(sqanti_summary)[c(1,5,6,9,12)])
ncol <- length(sqanti_summary[c(1,5,6,9,12)])

subset_gt4 <- tx_union %>% 
  filter(associated_transcript %in% rownames(feature_count)[rowSums(feature_count)>=(ncol-1)])
subset_shared <- tx_union %>% 
  filter(associated_transcript %in% rownames(feature_count)[rowSums(feature_count)==ncol])

lapply(names(filelist), function(x) {
  print(x)
  subset_known <- filelist[[x]]$sqanti_anno %>%
    mutate(gene_category = ifelse(structural_category %in% c('antisense','fusion','intergenic','genic_intron', NA), 
                                  structural_category,
                                  'known_gene')) %>%
    mutate(gene_category = ifelse(is.na(gene_category), 'not_mapped', gene_category)) %>%
    dplyr::filter(gene_category=='known_gene' ) %>% ## any de novo transcript mapped to known gene (FSM/ISM, NNC/NIC, genic)
    select(associated_transcript, associated_gene)
  
  df <- rbind(# get_union_f1(sqanti_anno = filelist[[x]]$sqanti_anno, truth = tx_union),
    get_union_f1(sqanti_anno = filelist[[x]]$sqanti_anno, truth = subset_gt4, remove.nomatch = F),
    get_union_f1(sqanti_anno = filelist[[x]]$sqanti_anno, truth = subset_shared, remove.nomatch = F),
    get_union_f1(sqanti_anno = filelist[[x]]$sqanti_anno, truth = subset_gt4),
    get_union_f1(sqanti_anno = filelist[[x]]$sqanti_anno, truth = subset_shared),
    get_union_f1(sqanti_anno = filelist[[x]]$sqanti_anno, truth = subset_known, keep.all.truth = F)) %>%
    mutate(clustering = x)
  df$method = c('all_but1_w_nomatch', 'shared_w_nomatch', 
                'all_but1_wo_nomatch', 'shared_wo_nomatch', 'de_novo')
  df
}) %>% rbindlist() %>% 
  write.csv('plot/new_cluster_metric_sequin.csv')

# calculate redundancy of de novo clusters

sqanti_summary_tmp <- sqanti_summary[c(1:3, rep(4:8, each = 2), 9:11, 12,12)]

cnt.list <- lapply(1:length(filelist), function(x) {
  calculate_denovo_cluster(cluster_summary = filelist[[x]], 
                           sequin = F)
}) 

lapply(cnt.list, function(x) x$counts) %>% rbindlist() %>%
  mutate(dataset = rep(names(filelist), each = 4)) %>%
  write.csv('plot/cluster_count_redundancy.csv')

lapply(1:length(filelist), 
       function(x) {
         length <- filelist[[x]]$sqanti_anno %>%
           mutate(m_rate = nm/length) %>%
           group_by(geneid) %>%
           summarise(mean_length = mean(length, na.rm = T), 
                     mean_rate = mean(m_rate, na.rm = T))
         
         cnt.list[[x]]$asm %>% 
           left_join(length, by ='geneid') %>%
           mutate(method = names(filelist)[x])}) %>% 
  rbindlist() %>% 
  filter(category=='match') %>% 
  pivot_longer(c('num_denovo_tx','mean_length','mean_rate')) %>% 
  ggplot(aes(x = match_type, y = value)) + 
  geom_boxplot() +
  facet_wrap(method ~ name, scales = 'free_y', ncol = 6) +
  stat_compare_means()
ggsave('plot/cluster_redundant_match_feature.pdf', width = 10, height = 20)


cnt.list.sequin <- lapply(1:length(filelist), function(x) {
  calculate_denovo_cluster(cluster_summary = filelist[[x]], 
                           sequin = T)
}) 

lapply(cnt.list.sequin, function(x) x$counts) %>% rbindlist() %>%
  mutate(dataset = rep(names(filelist), each = 4)) %>%
  write.csv('plot/cluster_count_redundancy_sequin.csv')

lapply(1:length(filelist), 
       function(x) {
         length <- filelist[[x]]$sqanti_anno %>%
           mutate(m_rate = nm/length) %>%
           group_by(geneid) %>%
           summarise(mean_length = mean(length, na.rm = T), 
                     mean_rate = mean(m_rate, na.rm = T))
         
         cnt.list.sequin[[x]]$asm %>% 
           left_join(length, by ='geneid') %>%
           mutate(method = names(filelist)[x])}) %>% 
  rbindlist() %>% 
  filter(category=='match') %>% 
  pivot_longer(c('num_denovo_tx','mean_length','mean_rate')) %>% 
  ggplot(aes(x = match_type, y = value)) + 
  geom_boxplot() +
  facet_wrap(method ~ name, scales = 'free_y', ncol = 6) +
  stat_compare_means()
ggsave('plot/cluster_redundant_match_feature_sequin.pdf', width = 10, height = 20)
