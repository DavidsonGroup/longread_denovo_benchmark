library(tximport)
library(bambu)
source('~/lab_davidson/yan.a/software/scripts_denovo/R/get_gene_tx.R')

get_exp_length_cor <- function(args){
  method <- str_split(args[3], pattern = '_', simplify = T)[2]
  
  dirs <- list.dirs(args[1], recursive = F, full.names = T)
  rt_gene <- get_gene_tx(args[2], 
                         list.files(dirs, pattern = 'classification.txt', full.names = T, recursive = T), 
                         method = method)
  rownames(rt_gene) <- rt_gene$isoform
  
  files <- list.files(str_replace(args[1], 'sqanti3','salmon'), 'quant.sf', recursive = T, full.names = T)[2]  
  salmon_quant_list <-tximport(files, 
                               type = "salmon", 
                               txOut = T)
  
  cor <- rt_gene %>% 
    left_join(data.frame(counts= salmon_quant_list$counts, isoform = rownames(salmon_quant_list$counts))) %>% 
    mutate(exp = log2(counts + 1)) %>% 
    group_by(geneid) %>% 
    select(isoform, geneid, exp, length) %>% 
    dplyr::summarise(cor = cor(length, exp#, method = 'spearman'
    ))  %>% pull(cor) 
  return(cor)
}

args <- c('../rnabloom2/merged_5m_sqanti3/', 
          '../rnabloom2/merged_5m_corset/merged_5m-clusters_mod.txt', 
          'rnabloom2_corset')
cor_rb <- get_exp_length_cor(args)

args <- c('../isonform/merged_5m_sqanti3/',
          '../isonform/merged_5m_corset/merged_5m-clusters_mod.txt',
          'isonform_corset')
cor_is <- get_exp_length_cor(args)

args <- c('../trinitystranded/merged_ilu_sqanti3/',
          '../trinitystranded/merged_ilu_corset/merged_ilu-clusters_mod.txt',
          'trinity_corset')
cor_tr <- get_exp_length_cor(args)

args <- c('../rattle/merged_5m_sqanti3/', 
          '../rattle/merged_5m_corset/merged_5m-clusters_mod.txt', 
          'rattle_corset')
cor_rt <- get_exp_length_cor(args)

# requested_counts <- read.table('../sqanti_sim/1m/sqanti-sim_index_meta.tsv', header = T, sep = '\t')
bambu_counts <- readRDS('../bambu/merged_5m_dge/default/se1.rds')

# bambu_counts <- read.table('../bambu/pea_dge/default/counts_transcript.txt', header = T, sep = '\t') 
requested_counts <- data.frame(sim_counts = rowSums(assays(bambu_counts)$counts),
                               transcript_id = rowData(bambu_counts)$TXNAME, 
                               gene_id = rowData(bambu_counts)$GENEID,
                               length = sapply(width(rowRanges(bambu_counts)), sum, simplify = T)) 
requested_counts <- requested_counts %>% filter(sim_counts != 0 )

cor_sim <- requested_counts %>% 
  mutate(exp = log2(sim_counts + 1)) %>% 
  group_by(gene_id) %>% 
  select(transcript_id, gene_id, exp, length) %>% 
  summarise(cor = cor(length, exp, #method = 'spearman'
  ))  %>% pull(cor) 

rbind(data.frame(method = 'rattle_corset', cor = cor_rt),
      data.frame(method = 'isonform_corset', cor = cor_is),
      data.frame(method = 'rnabloom2_corset', cor = cor_rb),
      data.frame(method = 'trinity_corset', cor = cor_tr),
      data.frame(method = 'sim', cor = cor_sim)) %>%
  ggplot(aes(x = method, y = cor)) +
  geom_violin() +
  geom_boxplot(width=0.1, color="grey", alpha=0.2)
ggsave('plot/cor_txexpression_length.pdf', width = 5, height = 3)
