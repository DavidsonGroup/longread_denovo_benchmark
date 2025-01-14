suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(clevr))
source('~/lab_davidson/yan.a/software/scripts_denovo/R/get_gene_tx.R')

# args[1] <- c('../simulation_1m/rattle/') # parent folder to get method and sqanti3 files
# args[2] <- '../simulation_1m/rattle/fx2tab.txt' # gene tx mapping
# args[3] <- 'rattle_rattle' # output prefix

# args[1] <- c('../trinity/') # parent folder to get method and sqanti3 files
# args[2] <- '../trinity/A549_MCF7_Illimina_merged/A549_MCF7_Illimina_merged_trinity.Trinity.fasta.gene_trans_map' # gene tx mapping
# args[3] <- 'trinity_trinity' # output prefix

args = commandArgs(trailingOnly=TRUE)

dir.create(args[3])
dir <- args[3]

method <- str_split(args[3], pattern = '_', simplify = T)[2]

dirs <- args[1]
# dirs <- list.dirs(args[1], recursive = F, full.names = T)
rt_gene <- get_gene_tx(args[2], 
                       list.files(dirs, pattern = 'classification.txt', full.names = T, recursive = T), 
                       method = method,
                       extra = list.files(dirs, pattern = 'read_score', full.names = T, recursive = F))
rownames(rt_gene) <- rt_gene$isoform

# known gene: "full-splice_match","incomplete-splice_match", "novel_in_catalog","novel_not_in_catalog",
# genic also has gene name, so define as known gene
# "fusion", "antisense","intergenic","genic_intron" either novel or ambiguous
rt_gene <- rt_gene %>%
  mutate(gene_category = ifelse(structural_category %in% c('antisense','fusion','intergenic','genic_intron', NA), 
                                structural_category,
                                'known_gene')) %>%
  mutate(gene_category = ifelse(is.na(gene_category), 'not_mapped', gene_category)) 

rt_gene$geneid %>% unique() %>% length()

# with NA
# eval_report_clusters(rt_gene$geneid, rt_gene$associated_gene) %>% unlist()

## remove NA and the result is slightly different
print('Metrics for all annotated assembled transcripts')

eval_report_clusters(rt_gene %>%
                       dplyr::filter(!is.na(associated_gene)) %>%
                       pull(geneid), 
                     rt_gene %>%
                       dplyr::filter(!is.na(associated_gene)) %>%
                       pull(associated_gene)) %>% unlist()

subset <- rt_gene %>% dplyr::filter(gene_category=='known_gene') 

# only known genes and the result is slightly different
eval <- eval_report_clusters(subset %>%
                       pull(geneid), 
                     subset %>%
                       pull(associated_gene)) %>% unlist()
print('Metrics for all assembled transcripts annotated with know gene')
eval

## get some more metrics
tptable <- subset %>% dplyr::count(associated_gene, geneid) %>%
  dplyr::filter(n >1) %>% 
  select(geneid, n)

TP <- sum(choose(tptable$n, 2))
# TP=sum(sapply(scl,function(x){ sum(choose(table(x),2)) }))
TPFP=sum(choose(table(subset$geneid),2))
TPFN=sum(choose(table(subset$associated_gene),2))
#show(paste("True Positives:",TP,"TPNP:",TPFP,"TPFN:",TPFN))

FP=TPFP-TP
FN=TPFN-TP
P=TP/TPFP
R=TP/TPFN
all=choose(dim(subset)[1],2)
TN=all-FP-FN-TP
RAND=(TP+TN)/all
Fscore=2*(R*P)/(P+R)
FM=sqrt((TP/TPFP)*(TP/TPFN))

# # set novel tx into separate groups (maybe under-clustering)
# eval_report_clusters(rt_gene %>%
#                        mutate(associated_transcript = ifelse(associated_transcript == 'novel',
#                                                              paste('novel', row_number(), sep = '_'),
#                                                              associated_transcript)) %>%
#                        pull(isoform), 
#                      rt_gene %>%
#                        mutate(associated_transcript = ifelse(associated_transcript == 'novel',
#                                                              paste('novel', row_number(), sep = '_'),
#                                                              associated_transcript)) %>%
#                        pull(associated_transcript)) %>% unlist()
# 
# # on known tx
# eval_report_clusters(rt_gene %>%
#                        dplyr::filter(associated_transcript != 'novel') %>%
#                        pull(isoform), 
#                      rt_gene %>%
#                        dplyr::filter(associated_transcript != 'novel') %>%
#                        pull(associated_transcript)) %>% unlist()


## count how many unique associated_gene in each de novo cluster (divided into different type of associated_gene) 
## multiple count for each de novo cluster
rt_gene %>% 
  group_by(geneid, gene_category) %>% 
  summarise(Unique_Count = n_distinct(associated_gene)) %>% 
  mutate(bin = cut(Unique_Count, breaks = c(0,1:10, Inf))) %>% 
  ggplot(aes(x=bin, fill=gene_category)) +
  geom_bar() +
  facet_grid(.~gene_category) +
  ylab('Number of de novo clusters') +
  xlab('Number of known genes')
ggsave(paste0(dir, '/number_of_clusters_within_genes.pdf'))

## count how many unique known associated_gene in each de novo cluster 
## <=1 count for each de novo cluster
pure_denovo <- subset %>% 
  filter(gene_category=='known_gene') %>%
  group_by(geneid, gene_category) %>% 
  summarise(Unique_Count = n_distinct(associated_gene)) %>%
  ungroup() %>%
  dplyr::count(Unique_Count)

subset %>% 
  filter(gene_category=='known_gene') %>%
  group_by(geneid, gene_category) %>% 
  summarise(Unique_Count = n_distinct(associated_gene)) %>% 
  mutate(bin = cut(Unique_Count, breaks = c(0,1:10, Inf))) %>% 
  ggplot(aes(x=bin, fill=gene_category)) +
  geom_bar() + 
  ylab('Number of de novo clusters') +
  xlab('Number of known genes')
ggsave(paste0(dir, '/number_of_clusters_within_knowngene.pdf'))

## how many associated_gene are split into multiple de novo clusters
rt_gene %>% 
  group_by(associated_gene, gene_category) %>% 
  summarise(Unique_Count = n_distinct(geneid)) %>% 
  mutate(bin = cut(Unique_Count, breaks = c(0,1:10, Inf))) %>% 
  ggplot(aes(x=bin, fill=gene_category)) +
  geom_bar() +
  facet_grid(.~gene_category) +
  ylab('Number of known genes') +
  xlab('Number of de novo clusters')
ggsave(paste0(dir, '/number_of_genes_within_clusters.pdf'))

## how many known associated_gene are split into multiple de novo clusters
pure_known <- subset %>% 
  filter(gene_category=='known_gene') %>%
  group_by(associated_gene, gene_category) %>% 
  summarise(Unique_Count = n_distinct(geneid)) %>%
  ungroup() %>%
  dplyr::count(Unique_Count)

subset %>% 
  filter(gene_category=='known_gene') %>%
  group_by(associated_gene, gene_category) %>% 
  summarise(Unique_Count = n_distinct(geneid)) %>% 
  mutate(bin = cut(Unique_Count, breaks = c(0,1:10, Inf))) %>% 
  ggplot(aes(x=bin, fill=gene_category)) +
  geom_bar() +
  ylab('Number of known genes') +
  xlab('Number of de novo clusters')
ggsave(paste0(dir, '/number_of_knowngenes_within_clusters.pdf'))

## How many clusters are perfect (one to one)

## get de novo clusters which entirely not mapped to any known genes
# geneid_nomap <- rt_gene %>%
#   dplyr::group_by(geneid) %>%
#   dplyr::summarise(NA_Count = sum(is.na(associated_gene)),
#             Total_Count = n()) %>%
#   filter(NA_Count == Total_Count) %>%
#   pull(geneid)

## assign de novo cluster with unique associated_gene
asm_map <- subset %>% 
  # filter(!(geneid %in% geneid_nomap)) %>% # remove nomapped geneid
  add_count(geneid, associated_gene) %>%
  mutate(n = ifelse(is.na(associated_gene), 0, n)) %>% ## give associated_gene == NA 0 count in case too many trasncript from gene cluster not assigned
  group_by(geneid) %>%
  dplyr::filter(n == max(n)) %>% ## majority vote
  dplyr::filter(length == max(length)) %>% ## if tie, only show gene that have longest transcript match (can also filter by RATTLE counts, but RNAbloom2 does not have count)
  slice_head(n = 1) %>% # get the first if there is still ties
  dplyr::select(geneid, associated_gene) %>% 
  ungroup() %>%
  distinct() %>%
  dplyr::rename(assigned_gene = associated_gene)

# ## SKIP THIS when geneid_nomap is empty
# if (!is.null(nrow(geneid_nomap))) {
#   asm_map <- rbind(asm_map,
#                    data.frame(geneid = geneid_nomap,
#                               assigned_gene = paste0('novel', 1:length(geneid_nomap)))) ## add no map geneid back
# }

asm_map <- data.frame(geneid = unique(rt_gene$geneid)) %>% 
  left_join(asm_map)

## get 1-to-1 mapped percentage of only known genes
asm_map_known <- asm_map %>%
  dplyr::filter(!str_detect(assigned_gene, '^novel'))

pct <- (asm_map_known %>% filter(!duplicated(assigned_gene)) %>% nrow()) / nrow(asm_map_known)

summary <- list(eval = eval, Precision = P, Recall = R, Rand = RAND, F1 = Fscore,
                sqanti_anno = rt_gene,
                pct_1to1 = pct, asm_map = asm_map,
                pure_denovo = pure_denovo, pure_known = pure_known,
                total_denovotx = dim(rt_gene)[1], total_denovocluster = length(unique(rt_gene$geneid)),
                known_denovotx = dim(subset)[1], known_denovocluster = length(unique(subset$geneid)))

cat(paste("RAND:", RAND, '\nFM:', FM, '\nF1:', Fscore,
           '\nPrecision:', P, '\nRecall:', R,
           '\n% 1 to 1:', pct,
          '\n'))

saveRDS(summary, paste0(dir, '/summary.rds'))
