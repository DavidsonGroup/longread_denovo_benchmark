# fusion
library(ggvenn)
library(RColorBrewer)
load("de_data.RData")

files <- list.files('../jaffal/', 'csv', recursive = T, full.names = T)
files <- files[!grepl(' ', files)]
files <- files[-c(1,4)] # remove isonform empty file and old trinity

fusion_list <- lapply(files, function(x){
  read.csv(x)$fusion.genes
})
names(fusion_list) <- str_remove(dirname(files), '../jaffal//')

pdf('plot/fusion_upset.pdf', width = 5, height = 5)
fromList(fusion_list) %>% upset()
dev.off()

union <- Reduce(c, fusion_list) %>% unique()

fusion_all <- data.frame(Fusion.Name = union,
           rattle = union %in% fusion_list$rattle,
           rnabloom2 = union %in% fusion_list$rnabloom2,
           trinity = union %in% fusion_list$trinity) %>% 
  mutate(is.A549 = Fusion.Name %in% str_replace(read.csv('../jaffal/A549 fusions.csv')$Fusion.Name, pattern = '--', replacement = ':'), 
         is.MCF7 = Fusion.Name %in% str_replace(read.csv('../jaffal/MCF7 fusions.csv')$Fusion.Name, pattern = '--', replacement = ':')) %>%
  mutate(fusion = ifelse(is.A549, 'A549', ifelse(is.MCF7, 'MCF7', 'FP'))) %>%
  mutate(fusion = factor(fusion, levels = c('FP','A549','MCF7')),
         num = rattle + rnabloom2 + trinity) %>%
  tidyr::unite(combo, rattle:trinity, remove = F)
write.csv(fusion_all, 'plot/fusion_all.csv')

data.table::rbindlist(list(count(fusion_all, rattle, fusion) %>% mutate(method = 'rattle'),
      count(fusion_all, rnabloom2, fusion) %>% mutate(method = 'rnabloom2'),
      count(fusion_all, trinity, fusion) %>% mutate(method = 'trinity'))) %>%
  filter(rattle == T) %>%
  ggplot(aes(x = method, y = n, fill= fusion)) +
  geom_bar(stat = 'identity') +
  xlab('Number of fusion genes')
ggsave('plot/num_fusion_genes.pdf', width = 5, height = 4)

ggvenn(fusion_list, show_elements = T, text_size = 2, label_sep = "\n", 
       fill_color = brewer.pal(name="Set2",n=3)) 
ggsave('plot/fusion_venn.pdf', width = 8, height = 8)

count(fusion_all, combo, fusion) %>%
  mutate(combo = factor(combo, labels = c('trinity','rnabloom2','rnabloom2:trinity','rattle','rattle:trinity','rattle:ranbloom2','all'))) %>%
  ggplot(aes(x = combo, y = n, fill= fusion)) +
  geom_bar(stat = 'identity') +
  xlab('Number of fusion genes') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave('plot/num_fusion_genes_venn.pdf', width = 5, height = 4)
  

for (i in c('rattle','rnabloom2','trinitystranded')) {
  
  contigs <- read.table(paste0('../jaffal/',i,'/transcriptome.reformat/transcriptome.reformat.txt'), 
                        sep = '\t', header = F) %>%
    filter(V4 %in% fusion_list[[i]])
  
  idx <- paste0(i, '_corset_map') %>% str_remove(., 'stranded')
  
  if (i != 'trinitystranded') {
    idx <- paste0(i, '_corset_onts')
  } 
  
  fusions <- test3[[idx]] %>%
    rowid_to_column() %>%
    mutate(associated_transcript = ifelse(associated_transcript == 'novel', 
                                          paste(idx, rowid, sep = '_'),
                                          associated_transcript),
           is.de =  adj.P.Val < 0.05, 
           is.bambu = as.numeric(associated_transcript %in% true_dte),
           is.fusion = as.numeric(isoform %in% contigs$V1)) %>% 
    left_join(contigs, by = c("isoform"='V1')) %>% 
    left_join(fusion_all, by = c('V4' = 'Fusion.Name')) %>%
    select(isoform, associated_transcript, is.de, is.bambu, is.fusion, V4, fusion) %>%
    # filter(is.de==1& !is.na(V4) & fusion != 'FP') %>% 
    filter(is.de==1& !is.na(V4) ) %>% 
    pull(V4) %>% unique()
  
  print(length(fusions))
  print(fusions)
  
}



