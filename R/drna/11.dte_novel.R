library(ggh4x)
library(tidyverse)

setwd("~/lab_davidson/yan.a/dRNA_gencode/R")

args <- c('~/vast_scratch/sqanti3_bambu/sqanti3_results/drna/isonform/',
          '~/vast_scratch/sqanti3_bambu/sqanti3_results/drna/rnabloom2/',
          '~/vast_scratch/sqanti3_bambu/sqanti3_results/drna/rattle/',
          '~/vast_scratch/sqanti3_bambu/sqanti3_results/drna/trinity/'
)

dirs <- list.dirs(args, recursive = F, full.names = T)

sqanti_files <- list.files(dirs, 'classification.txt', 
                           recursive = T, full.names = T)

sqanti_summary2 <- lapply(sqanti_files, function(x) {
  read.table(x, sep = '\t', header = T) %>%
    dplyr::select(c(1:8,15:17,25,29,30,41,46)) %>%
    dplyr::mutate(original = ifelse(str_detect(isoform, '_dup'), str_remove(isoform, '_dup.*'), isoform),
                  category = case_when(
                    # Logic: Is my root in the list of roots that have a duplicate somewhere in the WHOLE table?
                    original %in% original[str_detect(isoform, "_dup")] ~ 'chimera',
                    str_detect(associated_transcript, '^Bambu') ~ 'novel_in_bambu',
                    str_detect(associated_transcript, '^novel') | is.na(associated_transcript) ~ 'novel',
                    TRUE ~ 'reference'
                  ))
})
names(sqanti_summary2) <- c('isonform','rattle','rnabloom2','trinity')

dte <- lapply(list.files(pattern = 'dte_tx.rds', full.names = T, recursive = T), 
              function(x){
  readRDS(x) %>% 
    filter(adj.P.Val<0.05) 
})
names(dte) <- list.files(pattern = 'dte_tx.rds', full.names = T, recursive = T) %>% dirname() %>% str_remove('./')
dte <- dte[c('isonform_isonclust_onts','rattle_rattle_onts','rnabloom2_corset_onts', 'trinity_trinity_map')]

df <- lapply(1:4, function(x) {
  dte[[x]] %>% left_join(sqanti_summary2[[x]], by = 'isoform',
                    suffix = c('.gencode','.bambu')) %>%
    mutate(category = ifelse(is.na(category), 'unmapped', category)) %>% # umapped is still in dte, but not annotated in sqanti3
    count(category)
}) %>% purrr::reduce(left_join, by = 'category') 

colnames(df) <- c('category', names(sqanti_summary2) )
write.csv(df, 'plot/dte_novel_in_bambu.csv')

df %>%
  pivot_longer(2:5,
               names_to = "assembler",
               values_to = "count") %>%
  ggplot(aes(x = assembler, y = count, fill = category)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal(base_size = 14) +
  labs(x = NULL, y = "Count", fill = "Category") +
  scale_fill_brewer(palette = "Set2") +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  force_panelsizes(cols = unit(3, 'inch'), 
                   rows = unit(3, 'inch'))
ggsave('plot/dte_novel_in_bambu.pdf', width = 5, height = 5)

# venn

library(data.table)

fusion_tp <- c(str_replace(read.csv('../jaffal/A549 fusions.csv')$Fusion.Name, pattern = '--', replacement = ':'),
               str_replace(read.csv('../jaffal/MCF7 fusions.csv')$Fusion.Name, pattern = '--', replacement = ':'))

fusion.rnabloom2 <- fread('../jaffal/rnabloom2/transcriptome.reformat/transcriptome.reformat.txt') %>%
  mutate(is.A549 = V4 %in% str_replace(read.csv('../jaffal/A549 fusions.csv')$Fusion.Name, pattern = '--', replacement = ':'), 
         is.MCF7 = V4 %in% str_replace(read.csv('../jaffal/MCF7 fusions.csv')$Fusion.Name, pattern = '--', replacement = ':')) %>%
  mutate(fusion.type = ifelse(is.A549, 'A549', ifelse(is.MCF7, 'MCF7', 'FP'))) 

colnames(fusion.rnabloom2)[c(1,4)] <- c('denovo_txid','fusion')

fusion.rattle <- fread('../jaffal/rattle/transcriptome.reformat/transcriptome.reformat.txt') %>%
  mutate(is.A549 = V4 %in% str_replace(read.csv('../jaffal/A549 fusions.csv')$Fusion.Name, pattern = '--', replacement = ':'), 
         is.MCF7 = V4 %in% str_replace(read.csv('../jaffal/MCF7 fusions.csv')$Fusion.Name, pattern = '--', replacement = ':')) %>%
  mutate(fusion.type = ifelse(is.A549, 'A549', ifelse(is.MCF7, 'MCF7', 'FP'))) 

colnames(fusion.rattle)[c(1,4)] <- c('denovo_txid','fusion')

list(rnabloom2 = dte[[3]] %>% 
       left_join(sqanti_summary2[[3]], by = 'isoform',
                 suffix = c('.gencode','.bambu')) %>%
       filter(category != 'reference') %>%
       left_join(fusion.rnabloom2, by = c('isoform'='denovo_txid')) %>%
       mutate(fusionid = ifelse(is.na(fusion) , isoform, fusion)) %>%
       pull(fusionid), 
     rattle = dte[[2]] %>% 
       left_join(sqanti_summary2[[2]], by = 'isoform',
                 suffix = c('.gencode','.bambu')) %>%
       filter(category != 'reference') %>%
       left_join(fusion.rattle, by = c('isoform'='denovo_txid')) %>%
       mutate(fusionid = ifelse(is.na(fusion) , isoform, fusion)) %>%
       pull(fusionid),
     fusion = fusion_tp) %>%
  ggvenn::ggvenn(. ,
                 show_percentage = FALSE,
                 fill_color = c("#D55E00", "#009E73","#999999"),
                 stroke_size = 0.5, set_name_size = 4, fill_alpha = 0.8)
ggsave('plot/venn_dte_novel_ccle.pdf', width = 5, height = 5)


# same figure 6 layout
pdf('plot/dte_novel_ccle_split.pdf', width = 5, height = 5)
ggvenn::ggvenn(list(novel = sqanti_summary2$rnabloom2 %>% 
                      filter(category != 'reference') %>%
                      pull(isoform),
                    fusion_TP = fusion.rnabloom2 %>% 
                      filter(fusion.type != 'FP') %>% 
                      pull(denovo_txid),
                    dte_denovo = dte$rnabloom2_corset_onts$isoform),
               show_percentage = FALSE,
               fill_color = c("#999999", "#56B4E9", "#D55E00"),
               stroke_size = 0.5, set_name_size = 4, fill_alpha = 0.8) +
  ggtitle('RNAbloom2 transcripts')

ggvenn::ggvenn(list(novel = sqanti_summary2$rattle %>% 
                      filter(category != 'reference') %>%
                      pull(isoform),
                    fusion_TP = fusion.rattle %>% 
                      filter(fusion.type != 'FP') %>% pull(denovo_txid),
                    dte_denovo = dte$rattle_rattle_onts$isoform),
               show_percentage = FALSE,
               fill_color = c("#999999", "#56B4E9", "#D55E00"),
               stroke_size = 0.5, set_name_size = 4, fill_alpha = 0.8) +
  ggtitle('RATTLE transcripts')
dev.off()

# dRNA

chim <- read.table('../rnabloom2_chimeric.txt', header = F)
nm <- fread('../rnabloom2_nm.txt')
colnames(nm) <- c('denovo_txid', 'sam_flag', 'chr', 'start_genome', 'mapq', 'nm', 'start_clip_bp','end_clip_bp','sequence')

noalign <- nm %>% 
  mutate(is.chimeric = denovo_txid %in% chim$V1,
         mismatch = str_remove(nm, 'NM:i:') %>% as.numeric()) %>%
  filter(start_clip_bp > 400 | end_clip_bp > 400 | is.chimeric | mismatch > 200 | is.na(mismatch)) %>% 
  pull(denovo_txid) %>% unique()

ggvenn::ggvenn(list(noalign = noalign,
                    fusion = fusion.rnabloom2 %>% 
                      filter(fusion.type != 'FP') %>% pull(denovo_txid),
                    dte_denovo = dte$rnabloom2_corset_onts$isoform),
               show_percentage = FALSE,
               fill_color = c("#999999", "#56B4E9", "#D55E00"),
               stroke_size = 0.5, set_name_size = 4, fill_alpha = 0.8)

ggsave('plot/venn_noalign.pdf', width = 3, height = 3)

drna.combined <- left_join(nm %>% 
                             mutate(is.chimeric = denovo_txid %in% chim$V1,
                                    mismatch = str_remove(nm, 'NM:i:') %>% as.numeric())  , 
                           dte$rnabloom2_corset_onts, by = c('denovo_txid' = 'isoform')) %>%
  left_join(fusion.rnabloom2, by = 'denovo_txid')

drna.combined %>% 
  filter(adj.P.Val < 0.05, denovo_txid %in% noalign) %>% 
  arrange(adj.P.Val, P.Value) %>% 
  write.csv('plot/rnabloom2_noalign_dte.csv')
# "rb_372057" "rb_357990" "rb_366312" "rb_358307" "rb_361170"

ovlp.drna <- Reduce(intersect, list(noalign, fusion.rnabloom2$denovo_txid, dte$rnabloom2_corset_onts$isoform))

drna.combined %>% 
  filter(denovo_txid %in% ovlp.drna) %>% 
  arrange(denovo_txid, sam_flag) %>% 
  write.csv('plot/rnabloom2_noalign_fusion.csv')

# tx for BCAS4:BCAS3
drna.combined %>% 
  filter(denovo_txid %in% ovlp.drna) %>% 
  arrange(denovo_txid, sam_flag) %>% 
  filter(fusion=='BCAS4:BCAS3') %>% 
  pull(denovo_txid) %>% 
  unique()

# rattle
chim <- read.table('../rattle_chimeric.txt', header = F)
nm <- fread('../rattle_nm.txt')
colnames(nm) <- c('denovo_txid', 'sam_flag', 'chr', 'start_genome', 'mapq', 'nm', 'start_clip_bp','end_clip_bp','sequence')

noalign <- nm %>% 
  mutate(is.chimeric = denovo_txid %in% chim$V1,
         mismatch = str_remove(nm, 'NM:i:') %>% as.numeric()) %>%
  filter(start_clip_bp > 400 | end_clip_bp > 400 | is.chimeric | mismatch > 200 | is.na(mismatch)) %>% 
  pull(denovo_txid) %>% unique()

ggvenn::ggvenn(list(noalign = noalign,
                    fusion = fusion.rattle %>% 
                      filter(fusion.type != 'FP') %>% pull(denovo_txid),
                    dte_denovo = dte$rattle_rattle_onts$isoform),
               show_percentage = FALSE,
               fill_color = c("#999999", "#56B4E9", "#009E73"),
               stroke_size = 0.5, set_name_size = 4, fill_alpha = 0.8)

ggsave('plot/venn_noalign2.pdf', width = 3, height = 3)

drna.combined2 <- left_join(nm %>% 
                             mutate(is.chimeric = denovo_txid %in% chim$V1,
                                    mismatch = str_remove(nm, 'NM:i:') %>% as.numeric())  , 
                           dte$rattle_rattle_onts, by = c('denovo_txid' = 'isoform')) %>%
  left_join(fusion.rattle, by = 'denovo_txid')

drna.combined2 %>% 
  filter(adj.P.Val < 0.05, denovo_txid %in% noalign) %>% 
  arrange(adj.P.Val, P.Value) %>% # pull(denovo_txid) %>% unique() %>% length
  write.csv('plot/rattle_noalign_dte.csv')

ovlp.drna2 <- Reduce(intersect, list(noalign, fusion.rattle$denovo_txid, dte$rattle_rattle_onts$isoform))

drna.combined2 %>% 
  filter(denovo_txid %in% ovlp.drna2) %>% 
  arrange(denovo_txid, sam_flag) %>% #filter(fusion.type != 'FP') %>% pull(denovo_txid) %>% unique() %>% length
  write.csv('plot/rattle_noalign_fusion.csv')



ggvenn::ggvenn(list(rattle_fusion_gene = drna.combined2 %>% 
                      filter(denovo_txid %in% ovlp.drna2) %>% 
                      arrange(denovo_txid, sam_flag) %>% filter(fusion.type != 'FP') %>% 
                      pull(fusion) %>% unique() ,
                    rnabloom2_fusion_gene = drna.combined %>% 
                      filter(denovo_txid %in% ovlp.drna) %>% 
                      arrange(denovo_txid, sam_flag) %>% filter(fusion.type != 'FP') %>% 
                      pull(fusion) %>% unique()),
               show_percentage = FALSE,
               fill_color = c( "#009E73", "#D55E00"),
               stroke_size = 0.5, set_name_size = 4, fill_alpha = 0.8)
ggsave('plot/venn_fusion_gene.pdf', width = 3, height = 3)



# extra effort ignore now
# to make the venn diagrm show duplicated elements rather than only compare unique 
setwd('~/lab_davidson/yan.a/pea_fastq_merged/')

fusion_tp <- c(str_replace(read.csv('jaffal/A549 fusions.csv')$Fusion.Name, pattern = '--', replacement = ':'),
               str_replace(read.csv('jaffal/MCF7 fusions.csv')$Fusion.Name, pattern = '--', replacement = ':'))

fusion.rnabloom2 <- fread('jaffal/rnabloom2/transcriptome.reformat/transcriptome.reformat.txt') %>%
  mutate(is.A549 = V4 %in% str_replace(read.csv('jaffal/A549 fusions.csv')$Fusion.Name, pattern = '--', replacement = ':'), 
         is.MCF7 = V4 %in% str_replace(read.csv('jaffal/MCF7 fusions.csv')$Fusion.Name, pattern = '--', replacement = ':')) %>%
  mutate(fusion.type = ifelse(is.A549, 'A549', ifelse(is.MCF7, 'MCF7', 'FP'))) 

colnames(fusion.rnabloom2)[c(1,4)] <- c('denovo_txid','fusion')

dte.rnabloom2 <- readRDS('R/rnabloom2_corset_onts/dte_tx.rds') %>% 
  filter(adj.P.Val<0.05)

chim.rnabloom2 <- read.table('rnabloom2_chimeric.txt', header = F)
nm.rnabloom2 <- fread('rnabloom2_nm.txt')
colnames(nm.rnabloom2) <- c('denovo_txid', 'sam_flag', 'chr', 'start_genome', 'mapq', 'nm', 'start_clip_bp','end_clip_bp','sequence')

noalign.rnabloom2 <- nm.rnabloom2 %>% 
  mutate(is.chimeric = denovo_txid %in% chim.rnabloom2$V1,
         mismatch = str_remove(nm, 'NM:i:') %>% as.numeric()) %>%
  filter(start_clip_bp > 400 | end_clip_bp > 400 | is.chimeric | mismatch > 200 | is.na(mismatch)) %>% 
  pull(denovo_txid) %>% unique()

# dte_novel.rnabloom2 <- intersect(noalign.rnabloom2, dte.rnabloom2)
dte_novel.rnabloom2 <- left_join(dte.rnabloom2, fusion.rnabloom2, by = c('isoform'='denovo_txid')) %>%
  filter(isoform %in% noalign.rnabloom2) %>%
  mutate(fusionid = ifelse(is.na(fusion) | fusion.type == 'FP', isoform, fusion))


fusion.rattle <- fread('jaffal/rattle/transcriptome.reformat/transcriptome.reformat.txt') %>%
  mutate(is.A549 = V4 %in% str_replace(read.csv('jaffal/A549 fusions.csv')$Fusion.Name, pattern = '--', replacement = ':'), 
         is.MCF7 = V4 %in% str_replace(read.csv('jaffal/MCF7 fusions.csv')$Fusion.Name, pattern = '--', replacement = ':')) %>%
  mutate(fusion.type = ifelse(is.A549, 'A549', ifelse(is.MCF7, 'MCF7', 'FP'))) 

colnames(fusion.rattle)[c(1,4)] <- c('denovo_txid','fusion')

dte.rattle <- readRDS('R/rattle_rattle_onts/dte_tx.rds') %>% 
  filter(adj.P.Val<0.05) 

chim.rattle <- read.table('rattle_chimeric.txt', header = F)
nm.rattle <- fread('rattle_nm.txt')
colnames(nm.rattle) <- c('denovo_txid', 'sam_flag', 'chr', 'start_genome', 'mapq', 'nm', 'start_clip_bp','end_clip_bp','sequence')

noalign.rattle <- nm.rattle %>% 
  mutate(is.chimeric = denovo_txid %in% chim.rattle$V1,
         mismatch = str_remove(nm, 'NM:i:') %>% as.numeric()) %>%
  filter(start_clip_bp > 400 | end_clip_bp > 400 | is.chimeric | mismatch > 200 | is.na(mismatch)) %>% 
  pull(denovo_txid) %>% unique()

# dte_novel.rattle <- intersect(noalign.rattle, dte.rattle)

dte_novel.rattle <- left_join(dte.rattle, fusion.rattle, by = c('isoform'='denovo_txid')) %>%
  filter(isoform %in% noalign.rattle) %>%
  mutate(fusionid = ifelse(is.na(fusion) | fusion.type == 'FP', isoform, fusion))

