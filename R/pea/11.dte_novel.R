library(ggh4x)
library(tidyverse)

setwd("~/lab_davidson/yan.a/pea_fastq_merged/R/")

args <- c(
          '~/vast_scratch/sqanti3_bambu/sqanti3_results/pea/rattle/',
          '~/vast_scratch/sqanti3_bambu/sqanti3_results/pea/rnabloom2/',
          '~/vast_scratch/sqanti3_bambu/sqanti3_results/pea/trinity/'
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
names(sqanti_summary2) <- c(
                            'rattle',
                            'rnabloom2', 'trinity')

dte <- lapply(list.files(pattern = 'dte_tx.rds', full.names = T, recursive = T), 
              function(x){
                readRDS(x) %>% 
                  filter(adj.P.Val<0.05) 
              })
names(dte) <- list.files(pattern = 'dte_tx.rds', full.names = T, recursive = T) %>% dirname() %>% str_remove('./')
dte <- dte[ c(
              'rattle_rattle_onts', 
              'rnabloom2_corset_onts', 'trinity_trinity_map')]

df <- lapply(1:length(dte), function(x) {
  dte[[x]] %>% left_join(sqanti_summary2[[x]], by = 'isoform',
                         suffix = c('.gencode','.bambu')) %>%
    mutate(category = ifelse(is.na(category), 'unmapped', category)) %>% # umapped is still in dte, but not annotated in sqanti3
    count(category)
}) %>% purrr::reduce(left_join, by = 'category') 

colnames(df) <- c('category', names(sqanti_summary2) )
write.csv(df, 'plot/dte_novel_in_bambu.csv')

df %>%
  pivot_longer(-1,
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

tx_busco <- read.table('../busco_id_full.txt', sep = '\t', header = F)
colnames(tx_busco) <- 'buscoid'

tx_busco.rnabloom2 <- read.table('../rnabloom2/pea_busco_corrected/pea_busco_corrected/run_fabales_odb10/full_table.tsv', 
                                 sep = '\t',
                                 header = TRUE,       # Use the first line we read as the header
                                 skip = 2,            # Skip the first 2 metadata comment lines
                                 comment.char = "",   # Treat '#' as part of the header, not a comment
                                 fill = TRUE,
                                 quote = "") %>%
  separate(Sequence, c('denovo_txid','start_tx','end_tx'), remove = F, sep = ':|-')
colnames(tx_busco.rnabloom2) <- c('busco_id','category','Sequence', 'denovo_txid','start_tx','end_tx', 'score', 'Length',  'OrthoDB_url',     'Description')

tx_busco.rattle <- read.table('../rattle/pea_busco_corrected/pea_busco_corrected/run_fabales_odb10/full_table.tsv', 
                              sep = '\t',
                              header = TRUE,       # Use the first line we read as the header
                              skip = 2,            # Skip the first 2 metadata comment lines
                              comment.char = "",   # Treat '#' as part of the header, not a comment
                              fill = TRUE,
                              quote = "") %>%
  separate(Sequence, c('denovo_txid','start_tx','end_tx'), remove = F, sep = ':|-')
colnames(tx_busco.rattle) <- c('busco_id','category','Sequence', 'denovo_txid','start_tx','end_tx', 'score', 'Length',  'OrthoDB_url',     'Description')

list(rnabloom2 = dte[[2]] %>% 
       left_join(sqanti_summary2[[2]], by = 'isoform',
                 suffix = c('.gencode','.bambu')) %>%
       filter(category != 'reference') %>%
       left_join(tx_busco.rnabloom2, by = c('isoform'='denovo_txid')) %>%
       mutate(buscoid = ifelse(is.na(busco_id), isoform, busco_id)) %>%
       pull(buscoid), 
     rattle = dte[[1]] %>% 
       left_join(sqanti_summary2[[1]], by = 'isoform',
                 suffix = c('.gencode','.bambu')) %>%
       filter(category != 'reference') %>%
       left_join(tx_busco.rattle, by = c('isoform'='denovo_txid')) %>%
       mutate(buscoid = ifelse(is.na(busco_id), isoform, busco_id)) %>%
       pull(buscoid),
     BUSCO = tx_busco$buscoid) %>%
  ggvenn::ggvenn(. ,
                 show_percentage = FALSE,
                 fill_color = c("#D55E00", "#009E73","#999999"),
                 stroke_size = 0.5, set_name_size = 4, fill_alpha = 0.8)
ggsave('plot/venn_dte_novel_busco.pdf', width = 5, height = 5)

# same figure 6 layout
pdf('plot/dte_novel_busco_split.pdf', width = 5, height = 5)
ggvenn::ggvenn(list(novel = sqanti_summary2$rnabloom2 %>% 
                      filter(category != 'reference') %>%
                      pull(isoform),
                    busco = tx_busco.rnabloom2$denovo_txid,
                    dte_denovo = dte$rnabloom2_corset_onts$isoform),
               show_percentage = FALSE,
               fill_color = c("#999999", "#56B4E9", "#D55E00"),
               stroke_size = 0.5, set_name_size = 4, fill_alpha = 0.8) +
  ggtitle('RNAbloom2 transcripts')

ggvenn::ggvenn(list(novel = sqanti_summary2$rattle %>% 
                      filter(category != 'reference') %>%
                      pull(isoform),
                    busco = tx_busco.rattle$denovo_txid,
                    dte_denovo = dte$rattle_rattle_onts$isoform),
               show_percentage = FALSE,
               fill_color = c("#999999", "#56B4E9", "#D55E00"),
               stroke_size = 0.5, set_name_size = 4, fill_alpha = 0.8) +
  ggtitle('RATTLE transcripts')
dev.off()

# old analysis version
# pea

chim <- read.table('../rnabloom2_chimeric.tx', header = F)
nm <- fread('../rnabloom2_nm.txt')
colnames(nm) <- c('denovo_txid', 'sam_flag', 'chr', 'start_genome', 'mapq', 'nm', 'start_clip_bp','end_clip_bp','sequence')

noalign <- nm %>% 
  mutate(is.chimeric = denovo_txid %in% chim$V1,
         mismatch = str_remove(nm, 'NM:i:') %>% as.numeric()) %>%
  filter(start_clip_bp > 400 | end_clip_bp > 400 | is.chimeric | mismatch > 100 | is.na(mismatch)) %>% 
  pull(denovo_txid) %>% unique()

ggvenn::ggvenn(list(noalign = noalign,
                    busco = tx_busco.rnabloom2$denovo_txid,
                    dte_denovo = dte$rnabloom2_corset_onts$isoform),
               show_percentage = FALSE,
               fill_color = c("#999999", "#56B4E9", "#D55E00"),
               stroke_size = 0.5, set_name_size = 4, fill_alpha = 0.8)
ggsave('plot/venn_noalign_corrected.pdf', width = 3, height = 3)

pea.combined <- left_join(nm %>% 
                            mutate(is.chimeric = denovo_txid %in% chim$V1,
                                   mismatch = str_remove(nm, 'NM:i:') %>% as.numeric()) , 
                          dte$rnabloom2_corset_onts, by = c('denovo_txid' = 'isoform')) %>%
  left_join(tx_busco.rnabloom2, by = 'denovo_txid')

pea.combined %>% 
  filter(adj.P.Val < 0.05, denovo_txid %in% noalign) %>% 
  arrange(adj.P.Val, P.Value, denovo_txid) %>% 
  write.csv('plot/rnabloom2_noalign_corrected_dte.csv')

pea.combined %>% 
  filter(denovo_txid %in% Reduce(intersect, list(noalign = noalign,
                                                 fusion = tx_busco.rnabloom2$denovo_txid,
                                                 dte_denovo = dte$rnabloom2_corset_onts$isoform))) %>% 
  write.csv('plot/rnabloom2_noalign_corrected_busco.csv')

chim <- read.table('../rattle_chimeric.txt', header = F)
nm <- fread('../rattle_nm.txt')
colnames(nm) <- c('denovo_txid', 'sam_flag', 'chr', 'start_genome', 'mapq', 'nm', 'start_clip_bp','end_clip_bp','sequence')

noalign2 <- nm %>% 
  mutate(is.chimeric = denovo_txid %in% chim$V1,
         mismatch = str_remove(nm, 'NM:i:') %>% as.numeric()) %>%
  filter(start_clip_bp > 400 | end_clip_bp > 400 | is.chimeric | mismatch > 100 | is.na(mismatch)) %>% 
  pull(denovo_txid) %>% unique()

ggvenn::ggvenn(list(noalign = noalign2,
                    fusion = tx_busco.rattle$denovo_txid,
                    dte_denovo = dte$rattle_rattle_onts$isoform),
               show_percentage = FALSE,
               fill_color = c("#999999", "#56B4E9", "#009E73"),
               stroke_size = 0.5, set_name_size = 4, fill_alpha = 0.8)
ggsave('plot/venn_noalign2.pdf', width = 3, height = 3)

pea.combined2 <- left_join(nm %>% 
                            mutate(is.chimeric = denovo_txid %in% chim$V1,
                                   mismatch = str_remove(nm, 'NM:i:') %>% as.numeric()) , 
                          dte$rattle_rattle_onts, by = c('denovo_txid' = 'isoform')) %>%
  left_join(tx_busco.rattle, by = 'denovo_txid')

pea.combined2 %>% 
  filter(adj.P.Val < 0.05, denovo_txid %in% noalign2) %>% 
  arrange(adj.P.Val, P.Value, denovo_txid) %>% 
  write.csv('plot/rattle_noalign_corrected_dte.csv')

pea.combined2 %>% 
  filter(denovo_txid %in% Reduce(intersect, list(noalign = noalign2,
                                                 fusion = tx_busco.rattle$denovo_txid,
                                                 dte_denovo = dte$rattle_rattle_onts$isoform))) %>% 
  write.csv('plot/rattle_noalign_corrected_busco.csv')



ggvenn::ggvenn(list(rattle_busco_gene = pea.combined2 %>% 
                      filter(denovo_txid %in% 
                               Reduce(intersect, list(noalign = noalign2,
                                                      fusion = tx_busco.rattle$denovo_txid,
                                                      dte_denovo = dte$rattle_rattle_onts$isoform))) %>% 
                      arrange(denovo_txid, sam_flag) %>% 
                      pull(busco_id) %>% unique() ,
                    rnabloom2_busco_gene = pea.combined %>% 
                      filter(denovo_txid %in% 
                               Reduce(intersect, list(noalign = noalign,
                                                      fusion = tx_busco.rnabloom2$denovo_txid,
                                                      dte_denovo = dte$rnabloom2_corset_onts$isoform))) %>%
                      arrange(denovo_txid, sam_flag) %>% 
                      pull(busco_id) %>% unique()),
               show_percentage = FALSE,
               fill_color = c( "#009E73", "#D55E00"),
               stroke_size = 0.5, set_name_size = 4, fill_alpha = 0.8)
ggsave('plot/venn_busco_gene.pdf', width = 3, height = 3)






# extra effort ignore now
# to make the venn diagrm show duplicated elements rather than only compare unique 
setwd('~/lab_davidson/yan.a/pea_fastq_merged/')

tx_busco <- read.table('busco_id_full.txt', sep = '\t', header = F)
colnames(tx_busco) <- 'buscoid'

tx_busco.rnabloom2 <- read.table('rnabloom2/pea_busco_corrected/pea_busco_corrected/run_fabales_odb10/full_table.tsv', 
                                 sep = '\t',
                                 header = TRUE,       # Use the first line we read as the header
                                 skip = 2,            # Skip the first 2 metadata comment lines
                                 comment.char = "",   # Treat '#' as part of the header, not a comment
                                 fill = TRUE,
                                 quote = "") %>%
  separate(Sequence, c('denovo_txid','start_tx','end_tx'), remove = F, sep = ':|-')
colnames(tx_busco.rnabloom2) <- c('busco_id','category','Sequence', 'denovo_txid','start_tx','end_tx', 'score', 'Length',  'OrthoDB_url',     'Description')

dte.rnabloom2 <- readRDS('R/rnabloom2_corset_onts/dte_tx.rds') %>% 
  filter(adj.P.Val<0.05)

chim.rnabloom2 <- read.table('rnabloom2_chimeric.txt', header = F)
nm.rnabloom2 <- fread('rnabloom2_nm.txt')
colnames(nm.rnabloom2) <- c('denovo_txid', 'sam_flag', 'chr', 'start_genome', 'mapq', 'nm', 'start_clip_bp','end_clip_bp','sequence')

noalign.rnabloom2 <- nm.rnabloom2 %>% 
  mutate(is.chimeric = denovo_txid %in% chim.rnabloom2$V1,
         mismatch = str_remove(nm, 'NM:i:') %>% as.numeric()) %>%
  filter(start_clip_bp > 400 | end_clip_bp > 400 | is.chimeric | mismatch > 100 | is.na(mismatch)) %>% 
  pull(denovo_txid) %>% unique()

# dte_novel.rnabloom2 <- intersect(noalign.rnabloom2, dte.rnabloom2)
dte_novel.rnabloom2 <- left_join(dte.rnabloom2, 
                                 tx_busco.rnabloom2,
                                 by = c('isoform'='denovo_txid')) %>%
  filter(isoform %in% noalign.rnabloom2) %>%
  mutate(buscoid = ifelse(is.na(busco_id) , isoform, busco_id))


tx_busco.rattle <- read.table('rattle/pea_busco_corrected/pea_busco_corrected/run_fabales_odb10/full_table.tsv', 
                              sep = '\t',
                              header = TRUE,       # Use the first line we read as the header
                              skip = 2,            # Skip the first 2 metadata comment lines
                              comment.char = "",   # Treat '#' as part of the header, not a comment
                              fill = TRUE,
                              quote = "") %>%
  separate(Sequence, c('denovo_txid','start_tx','end_tx'), remove = F, sep = ':|-')
colnames(tx_busco.rattle) <- c('busco_id','category','Sequence', 'denovo_txid','start_tx','end_tx', 'score', 'Length',  'OrthoDB_url',     'Description')

dte.rattle <- readRDS('R/rattle_rattle_onts/dte_tx.rds') %>% 
  filter(adj.P.Val<0.05) 

chim.rattle <- read.table('rattle_chimeric.txt', header = F)
nm.rattle <- fread('rattle_nm.txt')
colnames(nm.rattle) <- c('denovo_txid', 'sam_flag', 'chr', 'start_genome', 'mapq', 'nm', 'start_clip_bp','end_clip_bp','sequence')

noalign.rattle <- nm.rattle %>% 
  mutate(is.chimeric = denovo_txid %in% chim.rattle$V1,
         mismatch = str_remove(nm, 'NM:i:') %>% as.numeric()) %>%
  filter(start_clip_bp > 400 | end_clip_bp > 400 | is.chimeric | mismatch > 100 | is.na(mismatch)) %>% 
  pull(denovo_txid) %>% unique()

# dte_novel.rattle <- intersect(noalign.rattle, dte.rattle)

dte_novel.rattle <- left_join(dte.rattle, 
                              tx_busco.rattle ,
                              by = c('isoform'='denovo_txid')) %>%
  filter(isoform %in% noalign.rattle) %>%
  mutate(buscoid = ifelse(is.na(busco_id), isoform, busco_id))

