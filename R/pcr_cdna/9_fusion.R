# fusion
library(tidyverse)
library(UpSetR)
load("de_data.RData")

files <- list.files('../jaffal/', 'csv', recursive = T, full.names = T)
files <- files[grepl('_', files)]
files <- files[-9]

fusion_list <- lapply(files, function(x){
  read.csv(x)$fusion.genes
})
names(fusion_list) <- str_remove(dirname(files), '../jaffal//')

pdf('plot/fusion_upset.pdf', width = 5, height = 5)
fromList(fusion_list) %>% upset(nsets = 9, sets = names(fusion_list), keep.order = T)

fromList(fusion_list) %>% upset(nsets = 2, sets = names(fusion_list)[1:2], keep.order = T)

fromList(fusion_list) %>% upset(nsets = 3, sets = names(fusion_list)[c(4,5,3)], keep.order = T)

fromList(fusion_list) %>% upset(nsets = 3, sets = names(fusion_list)[c(7,8,6)], keep.order = T)

dev.off()


union <- Reduce(c, fusion_list) %>% unique()

fusion_all <- data.frame(Fusion.Name = union,
                         sapply(fusion_list, function(x) union %in% x) ) %>% 
  mutate(is.NCIH1975 = Fusion.Name %in% str_replace(read.csv('../jaffal/NCIH1975 fusions.csv')$Fusion.Name, pattern = '--', replacement = ':'), 
         is.HCC827 = Fusion.Name %in% str_replace(read.csv('../jaffal/HCC827 fusions.csv')$Fusion.Name, pattern = '--', replacement = ':')) %>%
  mutate(fusion = ifelse(is.NCIH1975, 'NCIH1975', ifelse(is.HCC827, 'HCC827', 'FP'))) %>%
  mutate(fusion = factor(fusion, levels = c('FP','NCIH1975','HCC827')),
         num = rowSums(.[2:10])) %>%
  tidyr::unite(combo, 2:10, remove = F)

write.csv(fusion_all, 'plot/fusion_all.csv')

# data.table::rbindlist(list(count(fusion_all, rattle, fusion) %>% mutate(method = 'rattle'),
#                            count(fusion_all, rnabloom2, fusion) %>% mutate(method = 'rnabloom2'),
#                            count(fusion_all, trinity, fusion) %>% mutate(method = 'trinity'))) %>%
#   filter(rattle == T) %>%
#   ggplot(aes(x = method, y = n, fill= fusion)) +
#   geom_bar(stat = 'identity') +
#   xlab('Number of fusion genes')
# ggsave('plot/num_fusion_genes.pdf', width = 5, height = 4)

# ggvenn(fusion_list, show_elements = T, text_size = 2, label_sep = "\n", 
#        fill_color = brewer.pal(name="Set2",n=3)) 
# ggsave('plot/fusion_venn.pdf', width = 8, height = 8)

count(fusion_all, combo, fusion) %>%
  # mutate(combo = factor(combo, labels = c('trinity','rnabloom2','rnabloom2:trinity','rattle','rattle:trinity','rattle:ranbloom2','all'))) %>%
  ggplot(aes(x = combo, y = n, fill= fusion)) +
  geom_bar(stat = 'identity') +
  xlab('Number of fusion genes') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave('plot/num_fusion_genes_venn.pdf', width = 5, height = 8)


for (i in 1:9) {
  
  name <- names(fusion_list)[i]
  print(name)
  files <- list.files(paste0('../jaffal/',name), pattern = 'ed.txt', recursive = T, full.names = T)
  contigs <- read.table(files, sep = '\t', header = F) #%>%
    #filter(V4 %in% fusion_list[[i]])
  
  idx <- paste0(str_split(name, '_', simplify = T)[1] %>% str_remove(., 'stranded'), 
                '_corset_10m_map')
  
  if (!str_detect(name, 'trinitystranded')) {
    idx <- paste0(str_split(i, '_', simplify = T)[1], '_corset_',str_split(i, '_', simplify = T)[2],'_onts')
  } 
  
  idx <- str_remove(idx, '../jaffal//')
  
  dtefile <- files.dte[grep('onts|map', files.dte)] 
    
  pattern <- paste(c(str_split(name, pattern = '_', simplify = T)[1], 
                     'corset', 
                     str_split(name, pattern = '_', simplify = T)[2]), 
                   collapse = '_')
  
  dtefile <- dtefile[grep(pattern, dtefile)]
  
  # check if DTE overlap with fusion
  fusions <- readRDS(dtefile) %>% 
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
    filter(is.de==1 & !is.na(V4) #& fusion != 'FP'
           ) %>% 
    pull(V4) %>% unique()
  
  print(length(fusions))
  print(fusions)
  
}



