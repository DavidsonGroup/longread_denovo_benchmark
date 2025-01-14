
library(tidyverse)
library(plyranges)
library(purrr)
library(scales)
library(fastmatch)
library(UpSetR)
library(data.table)
library(ggpubr)
library(tximport)
library(edgeR)
library(cowplot)
#library(bambu)

source('~/lab_davidson/yan.a/software/scripts_denovo/R/upset_sqanti.R')
source("~/lab_davidson/yan.a/software/scripts_denovo/R/quantile_overlap.R")

# set plotting parameters
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cols <- cbPalette[c(1,2,6,4,7,8)]
names(cols) <- c('limma','bambu','isonform','rattle','rnabloom2','trinity')
shapes <- c(15,15,16,17,17,17)
names(shapes) <- c('sim','bambu','corset','isonclust','rattle','trinity')

col_sqanti <- c(hue_pal()(3), rep(hue_pal()(3)[1], 2), rep(hue_pal()(3)[2], 2), rep(hue_pal()(3)[3], 5))
names(col_sqanti) <- c("Annotated","Novel","Sequin" ,
                       "full-splice_match","incomplete-splice_match",
                       "novel_in_catalog","novel_not_in_catalog",
                       "fusion","genic",
                       "antisense","intergenic","genic_intron" ) 
col_sqanti <- col_sqanti[c(2,1,3, 8:12,6:7,4:5)]

args <- commandArgs(trailingOnly=TRUE)

args <- c('../simulation_1m/bambu/',
          '../simulation_1m/rattle/', 
          '../simulation_1m/rnabloom2/',
          '../simulation_1m/trinity/',
          '../simulation_1m/isonform/'
          #'../bambu/',
          #'../isonform/'
)

dirs <- list.dirs(args, recursive = F, full.names = T)

sqanti_files <- list.files(dirs[grepl('sqanti3$', dirs)], 'unstranded_classification.txt', 
                           recursive = T, full.names = T)

sqanti_summary <- lapply(sqanti_files, function(x) {
  read.table(x, sep = '\t', header = T) %>%
  dplyr::select(c(1:8,15:17,25,29,30,41,46)) %>%
  dplyr::mutate(original = ifelse(str_detect(isoform, '_dup'), str_remove(isoform, '_dup.*'), isoform))
})
names(sqanti_summary) <- c('bambu','isonform', 'rattle','rnabloom2','trinity')

num_per_chr <- lapply(sqanti_summary, function(x) {
  x %>% dplyr::count(chrom) %>% 
    dplyr::filter(str_detect(chrom, 'chr')) ## only canonical chromosome
}) %>% purrr::reduce(full_join, by = 'chrom') 

colnames(num_per_chr) <- c('chrom', names(sqanti_summary) )

# read_gff('/vast/projects/lab_davidson/yan.a/ref/gencode/gencode.v44.annotation.gtf') %>%
#   saveRDS('data/gencodev44.rds')

num_ref <- read_gff('../sqanti_sim/gencode.v44.subset.gtf') %>%
  dplyr::filter(type == 'transcript') %>% 
  data.frame %>% 
  dplyr::count(seqnames)

num_per_chr <- num_per_chr %>%
  left_join(num_ref, by = c('chrom'='seqnames')) %>% 
  dplyr::rename(ref = n) 

num_per_chr %>%
  dplyr::mutate(chrom = factor(chrom, levels = paste0('chr', c(1:22, 'X', 'Y', 'M')))) %>%
  # dplyr::filter(chrom != 'chrY') %>% 
  pivot_longer(-1) %>% 
  #separate(name, c('assembler','depth')) %>%
  ggplot(aes(x = chrom, y = value, color = name)) + 
  geom_point() +
  geom_line(aes(group=name)) +
  # facet_wrap(.~assembler, scales = 'free_y') +
  scale_color_manual(values = cols) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  ylab('Number of transcripts')
ggsave('plot/sqanti_transcript_per_chr.pdf', width = 5, height = 3)

(num_per_chr[,2:(1+length(args))]/num_per_chr$ref) %>%
  mutate(chrom = factor(num_per_chr$chrom, levels = paste0('chr', c(1:22, 'X', 'Y', 'M')))) %>%
  # filter(!(chrom == 'chrY')) %>%
  pivot_longer(1:length(args)) %>% 
  #separate(name, c('assembler','depth')) %>%
  ggplot(aes(x = chrom, y = value, color = name)) + 
  geom_point() +
  geom_line(aes(group=name)) + 
  #facet_wrap(.~assembler, scales = 'free_y') +
  scale_color_manual(values = cols) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  scale_y_log10() +
  ylab('De novo / Reference')
ggsave('plot/sqanti_transcript_odds_per_chr.pdf', width = 5, height = 3)


chimeric_summary <- sapply(sqanti_summary, function(x){
  num =  x$isoform %>% str_detect('_dup') %>% sum
  pct =  x$isoform %>% str_detect('_dup') %>% mean
  c(num, pct)
},
       simplify = T) %>% t() %>% data.frame() %>%
  dplyr::rename(number = X1, pct = X2) %>%
  mutate(method = names(sqanti_summary))
  
chimeric_summary %>%   
  pivot_longer(1:2) %>% 
  #separate(method, c('assembler','depth')) %>%
  #mutate(depth = factor(depth, levels = c('1m','2m','5m','10m'))) %>%
  ggplot(aes(x = method, y = value, color = method)) + 
  geom_point() +
  geom_line(aes(group=method)) +
  scale_color_manual(values = cols) +
  facet_wrap(.~name, scales = 'free_y') +
  theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave('plot/sqanti_number_chimeric.pdf', width = 5, height = 3)



chimeric_summary_chr <- lapply(sqanti_summary, function(x){
  x %>% dplyr::filter(str_detect(isoform, '_dup')) %>% 
    dplyr::count(chrom) %>% 
    dplyr::filter(str_detect(chrom, 'chr'))
}) %>% purrr::reduce(full_join, by = 'chrom') %>%
  replace(is.na(.), 0)

colnames(chimeric_summary_chr)[-1] <-  names(sqanti_summary)
  
chimeric_summary_chr %>% 
  pivot_longer(-1) %>%
  mutate(chrom = factor(chrom, levels = paste0('chr', c(1:22, 'X', 'Y', 'M')))) %>% 
  #separate(name, c('assembler','depth')) %>%
  ggplot(aes(x = chrom, y = value, color = name)) + 
  geom_point() +
  geom_line(aes(group=name)) +
  scale_color_manual(values = cols) +
  #facet_wrap(.~assembler, scales = 'free_y') +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  ylab('Number of chimeric transcripts')
ggsave('plot/sqanti_chimeric_per_chr.pdf', width = 5, height = 3)




chimeric_subcat <- lapply(sqanti_summary, function(x){
  x %>% dplyr::filter(str_detect(isoform, '_dup')) %>% 
    dplyr::count(structural_category, chrom) %>%
    dplyr::filter(str_detect(chrom, 'chr')) %>%
    dplyr::mutate(chrom = factor(chrom, levels = paste0('chr', c(1:22, 'X', 'Y', 'M'))))
}) %>% purrr::reduce(full_join, by = c('chrom', 'structural_category')) 

colnames(chimeric_subcat)[-(1:2)] <- names(sqanti_summary)
  
chimeric_subcat %>% 
  pivot_longer(-(1:2)) %>%
  #separate(name, c('assembler','depth')) %>%
  #mutate(depth = factor(depth, levels = c('1m','2m','5m', '10m'))) %>%
  ggplot(aes(x = chrom, y = value, fill = structural_category)) + 
  geom_bar(stat = 'identity', position = 'stack') +
  facet_wrap(.~name, scales = 'free_y') +
  theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave('plot/sqanti_chimeric_per_chr_category.pdf', width = 10, height = 4)


sum_df <- lapply(sqanti_summary, function(x){
  anno.gene <- x %>% dplyr::count(associated_gene) %>%
    dplyr::mutate(category = ifelse(str_detect(associated_gene, '^ENSG'), 'Annotated', ## potential fusion/read through, very few though
                                     ifelse(str_detect(associated_gene, '^novel'), 'Novel', 'Others'))) %>%
    dplyr::count(category) %>%
    mutate(classification = 'gene') %>%
    dplyr::select(c(3,1,2))
  
  anno.tx <- x %>% dplyr::count(structural_category) %>% 
    mutate(classification = 'transcript') %>%
    dplyr::rename(category = structural_category) %>%
    dplyr::select(c(3,1,2))
  
  df <- rbind(anno.gene, anno.tx)
  
}) %>%
  purrr::reduce(full_join, by = c('classification','category')) 

# sum_df <- list.files(args, pattern = 'sqantisummary.txt', full.names = T)  %>%
#   lapply(., function(x) read.table(x, sep = '\t', header = T)) %>%
#   purrr::reduce(full_join, by = c('classification','category')) 

colnames(sum_df)[-(1:2)] <- names(sqanti_summary)

write.csv(sum_df, 'plot/sum_df.csv')

sum_df %>% 
  pivot_longer(-(1:2)) %>%
  #separate(name, c('assembler','depth'), remove = F) %>%
  mutate(# depth = factor(depth, levels = c('1m','2m','5m','10m')), 
         # name = factor(name, levels = c(
         #                   paste('RNAbloom2', c('1m','2m','5m','10m'), sep = '_'),
         #                   paste('rattle', c( '1m','2m','5m', '10m'), sep = '_'),
         #                   paste('trinityshort', c('1m','2m', '5m','10m'), sep = '_'))),
         category = factor(category, levels = names(col_sqanti))) %>%
  ggplot(aes(x=name, y=value, fill = category)) +
  geom_bar(stat = 'identity') + 
  geom_text(aes(label = value), colour = "white", size = 3.5,
            position=position_stack(vjust=0.5)) +
  facet_wrap(.~classification, scales = 'free_y') +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  scale_fill_manual(values = col_sqanti)
ggsave('plot/sqanti_summary_count.pdf', width = 10, height = 4)

sum_df %>% 
  pivot_longer(-(1:2)) %>%
  mutate(category = factor(category, levels = names(col_sqanti))) %>%
  ggplot(aes(x=name, y=value, fill = category)) +
  geom_bar(stat = 'identity') + 
  facet_wrap(.~classification, scales = 'free_y') +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  scale_fill_manual(values = col_sqanti) +
  ylab('Number')
ggsave('plot/sqanti_summary_count_simplified.pdf', width = 6, height = 2.5)

sum_df %>% 
  pivot_longer(-(1:2)) %>%
  # separate(name, c('assembler','depth'), remove = F) %>%
  mutate(# depth = factor(depth, levels = c('1m','2m','5m','10m')), 
         # name = factor(name, levels = c(
         #                   paste('RNAbloom2', c('1m','2m','5m','10m'), sep = '_'),
         #                   paste('rattle', c( '1m','2m','5m', '10m'), sep = '_'),
         #                   paste('trinityshort', c('1m','2m', '5m','10m'), sep = '_'))),
         category = factor(category, levels = names(col_sqanti))) %>%
  ggplot(aes(x=name, y=value, fill = category)) +
  geom_bar(stat = 'identity', position="fill") + 
  facet_wrap(.~classification, scales = 'free_y') +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  scale_fill_manual(values = col_sqanti)
ggsave('plot/sqanti_summary_pct.pdf', width = 10, height = 4)



pdf('plot/sqanti_upset.pdf', width = 6, height = 4)

## only overlap of know/annotated transcripts and genes
upset_sqanti(sqanti_summary, subcat = 'all', 
             names = names(sqanti_summary) , featuretype = 'gene') %>%
  upset(nsets = 5, order.by = 'freq', set_size.show = T, set_size.scale_max = 20000)
upset_sqanti(sqanti_summary, subcat = 'all', 
             names = names(sqanti_summary) , featuretype = 'transcript') %>%
  upset(nsets = 5, order.by = 'freq', set_size.show = T, set_size.scale_max = 50000)
dev.off()

# get reference expression and transcript length based on bambu
# se1_quant <- readRDS('../bambu/se1_quant.rds')
# txlength <- readRDS('../bambu/txlength.rds')

## load requested true count and DTE, DTU
requested_counts <- read.table('../sqanti_sim/1m/sqanti-sim_index_meta.tsv', header = T, sep = '\t')

true_dtugene <- readRDS("/vast/projects/davidson_longread/yan.a/simulation_20240501/sqanti_sim/design/true_dtugene.rds")
true_dge <- readRDS("/vast/projects/davidson_longread/yan.a/simulation_20240501/sqanti_sim/design/true_dge.rds")
true_dte <- readRDS("/vast/projects/davidson_longread/yan.a/simulation_20240501/sqanti_sim/design/true_dte.rds")
true_dtutx <- readRDS("/vast/projects/davidson_longread/yan.a/simulation_20240501/sqanti_sim/design/true_dtutx.rds")

requested_counts$dtu <- requested_counts$transcript_id %in% true_dtutx
requested_counts$dte <- requested_counts$transcript_id %in% true_dte

requested_counts <- requested_counts %>% filter(sim_counts != 0 ) %>%
  mutate(truecount_quantile = ntile(sim_counts, 10)) %>% # add quantiles for counts, equal number of elements in each quantile
  group_by(truecount_quantile) %>%
  mutate(truecount_quantile_range = paste0("[", round(min(sim_counts), 2), "-", round(max(sim_counts), 2), "]")) %>% # some ties are bread into each bin
  ungroup()

# get salmon quant for each assembly

## for ONT, primary only or include secondary
## note, for ratlle, we manually made quant.sf using rattle counts out of box
## for Ilu, align or map mode

files <- list.files(dirs[grepl('salmon', dirs)], 'quant.sf', recursive = T, full.names = T)  

salmon_quant_list <- lapply(files, function(x){
  tximport(x, 
           type = "salmon", 
           txOut = T)
})

salmon_count <- lapply(salmon_quant_list, function(x) {
  df <- data.frame(counts = x$counts[,1], isoform = rownames(x$counts))
})

names(salmon_count) <- str_split(files, '/', simplify = T)[,6]

saveRDS(salmon_count, 'plot/salmon_count.rds')
saveRDS(requested_counts, 'plot/requested_counts.rds')

## use Salmon expression for all assembly, compared to Bambu expression

sqanti_summary_tmp <- sqanti_summary[c(2,2,3,3,3,4,4,5,5)]

quant_summary <- list() 

# for ONT
quant_summary[1:7] <- lapply(1:7, function(x){
  
  title <- (files[x] %>% str_split('/', simplify = T))[3]
  print(title)
  
  quantile_overlap(sqanti_summary_tmp[[x]], salmon_count[[x]], title, 
                   requested_counts)
  
})
# for Illumina
illumina_counts <- lapply(list.files('../simulation_1m/trinity/newcounts/', pattern = 'counts.txt', full.names = T), 
                          function(x) {
                            read.csv(x, header = F)
                          }) %>% purrr::reduce(full_join, by='V2') %>%
  mutate(across(everything(), .fns = ~replace_na(.,0))) %>%
  dplyr::select(c(2,1,3:7))

illumina_counts <- data.frame(illumina_counts[,1], rowSums(illumina_counts[,-1]))
colnames(illumina_counts) <- c('transcript_id', 'illumina_counts_real')

quant_summary[8:9] <- lapply(8:9, function(x){
  
  title <- (files[x] %>% str_split('/', simplify = T))[3]
  print(title)
  
  quantile_overlap(sqanti_summary_tmp[[x]], salmon_count[[x]], title, 
                   requested_counts )
  
})

names(quant_summary) <- basename(dirname(files)) 

# can add recovery and correlation for bambu if needed
bambucount <- read.table('../simulation_1m/bambu/ONT_merged_dge/default/counts_transcript.txt', header = T)
bambucount$counts <- rowSums(bambucount[,3:8])
bambucount <- data.frame(counts = bambucount$counts, isoform = bambucount$TXNAME)
rownames(bambucount) <- bambucount$isoform

quant_summary$ONT_merged_quant_bambu_count <- quantile_overlap(sqanti_summary$bambu, bambucount, 'bambu', 
                                                                requested_counts)

quant_summary_cpm_recov <- lapply(quant_summary, function(x){
  x$lcpm_tx %>% 
    # group_by(truecpm_bins) %>% 
    dplyr::count(truecpm_bins, isassemble) %>%
    group_by(truecpm_bins) %>%
    summarise(
      Total = sum(n),
      Quantified = sum(n[isassemble == "quantified"]),
      ProportionQuantified = Quantified / Total
    )
}) %>% rbindlist() %>% 
  mutate(method = rep(names(quant_summary), each = 11)) %>%
  separate(method, c(NA,NA, NA, 'assembler', 'quant'), remove = F, sep = '_') %>%
  mutate(assembler = factor(assembler, levels = c('bambu','rattle','rnabloom2','isonform','trinity')))

quant_summary_cpm_recov %>% 
  filter(quant %in% c('onts', 'map') | assembler == 'bambu') %>%
  ggplot(aes(x = truecpm_bins, y = ProportionQuantified, color = assembler)) +
  geom_point() + 
  geom_line(aes(group = method)) +
  scale_color_manual(values = cols) +
  xlab('log2(CPM+1) bins on simulated counts long read data')
ggsave('plot/quant_summary_bin_lineplot.pdf', width = 8, height = 5)

# count quantiles (same as cpm quantile, but show the value in raw counts)
quant_summary_quantile_recov <- lapply(quant_summary, function(x){
  x$lcpm_tx %>% 
    dplyr::count(truecount_quantile, truecount_quantile_range, isassemble) %>%
    group_by(truecount_quantile, truecount_quantile_range) %>%
    summarise(
      Total = sum(n),
      Assembled = sum(n[isassemble != "not assembled"]),
      Quantified = sum(n[isassemble == "quantified"]),
      ProportionQuantified = Quantified / Total,
      ProportionAssembled = Assembled / Total,
      .groups = 'drop'
    ) %>%
    mutate(truecount_quantile_range = factor(truecount_quantile_range,
                                             levels = truecount_quantile_range))
}) %>% rbindlist() %>% 
  mutate(method = rep(names(quant_summary), each = 10)) %>% 
  separate(method, c(NA,NA, NA, 'assembler', 'quant'), remove = F, sep = '_') %>%
  mutate(assembler = factor(assembler, levels = c('bambu','rattle','rnabloom2','isonform','trinity')))
write.csv(quant_summary_quantile_recov, 'plot/quant_summary_quantile_recov.csv')

quant_summary_quantile_recov %>% 
  filter(quant %in% c('onts', 'map') | assembler == 'bambu') %>%
  ggplot(aes(x = truecount_quantile_range, y = ProportionAssembled, color = assembler)) +
  geom_point() + 
  geom_line(aes(group = method)) +
  scale_color_manual(values = cols) +
  ylim(c(0,1)) + 
  xlab('Quantiles on simulated counts long read data')
ggsave('plot/quant_summary_quantile_lineplot.pdf', width = 8, height = 5)

# gene at quantile level
quant_summary_quantile_recov2 <- lapply(quant_summary, function(x){
  x$lcpm_gene %>% 
    mutate(isassemble = ifelse(!isassemble, 'not assembled', ifelse(counts == 0, 'not quantified', 'quantified'))
    ) %>% # change not assembled to 0
    mutate(isassemble = factor(isassemble, levels = c('not assembled','not quantified','quantified'))) %>%
    dplyr::count(truecount_quantile, truecount_quantile_range, isassemble) %>%
    group_by(truecount_quantile, truecount_quantile_range) %>%
    summarise(
      Total = sum(n),
      Assembled = sum(n[isassemble != "not assembled"]),
      Quantified = sum(n[isassemble == "quantified"]),
      ProportionQuantified = Quantified / Total,
      ProportionAssembled = Assembled / Total,
      .groups = 'drop'
    ) %>%
    mutate(truecount_quantile_range = factor(truecount_quantile_range,
                                             levels = truecount_quantile_range))
}) %>% rbindlist() %>% 
  mutate(method = rep(names(quant_summary), each = 10)) %>% 
  separate(method, c(NA,NA, NA, 'assembler', 'quant'), remove = F, sep = '_') %>%
  mutate(assembler = factor(assembler, levels = c('bambu','rattle','rnabloom2','isonform','trinity')))
write.csv(quant_summary_quantile_recov2, 'plot/quant_summary_quantile_recov_gene.csv')

quant_summary_quantile_recov2 %>% 
  filter(quant %in% c('onts', 'map') | assembler == 'bambu') %>%
  ggplot(aes(x = truecount_quantile_range, y = ProportionAssembled, color = assembler)) +
  geom_point() + 
  geom_line(aes(group = method)) +
  scale_color_manual(values = cols) +
  ylim(c(0,1)) + 
  xlab('Quantiles on Bambu counts in 10 million data')
ggsave('plot/quant_summary_quantile_gene_lineplot.pdf', width = 8, height = 5)

# quant_summary_cpm_redun <- lapply(names(quant_summary), function(x){
#   quant_summary[[x]]$lcpm_tx %>% 
#     dplyr::mutate(redundantAsm_bin = cut(redundantAsm, breaks = c(0,1,2,5,10,Inf), right = F)) %>% 
#     dplyr::count(truecpm_bins, redundantAsm_bin) %>%
#     group_by(truecpm_bins) %>%
#     mutate(Total = sum(n)) %>%
#     ungroup() %>%
#     mutate(Proportion = n / Total,
#            method = x)
# }) %>% rbindlist() %>%
#   separate(method, c(NA,NA, NA, 'assembler', 'quant'), remove = F, sep = '_') %>%
#   mutate(assembler = factor(assembler, levels = c('bambu','rattle','rnabloom2','isonform', 'trinity')))
# 
# quant_summary_cpm_redun %>% 
#   filter(quant %in% c('onts', 'map') | assembler == 'bambu') %>%
#   ggplot(aes(x = truecpm_bins, y = Proportion, color = assembler)) +
#   geom_point() + 
#   geom_line(aes(group = method)) +
#   facet_grid(redundantAsm_bin~.) +
#   ylab(c(0,1)) +
#   scale_color_manual(values = cols)
# ggsave('plot/quant_summary_redundancy_lineplot.pdf', width = 8, height = 8)

## simplify redundancy plot

quant_summary_cpm_redun <- lapply(names(quant_summary), function(x){
  quant_summary[[x]]$lcpm_tx %>%
    mutate(method = x)
}) %>% rbindlist() %>%
  separate(method, c(NA,NA, NA, 'assembler', 'quant'), remove = F, sep = '_') %>%
  mutate(assembler = factor(assembler, levels = c('bambu','rattle','rnabloom2','isonform', 'trinity')),
         truecount_quantile_range = factor(truecount_quantile_range,
                                           levels = levels(quant_summary_quantile_recov$truecount_quantile_range)))

quant_summary_cpm_redun %>%
  filter(quant %in% c('onts', 'map'),
         # redundantAsm > 1#,
         #assembler == 'rattle', truecpm_bins == '(10,Inf]'
  ) %>%
  ggplot(aes(x = truecpm_bins, y = redundantAsm, color = assembler)) +
  geom_boxplot() +
  scale_y_log10() +
  scale_color_manual(values = cols) +
  xlab('log2(CPM+1) bins on simulated counts long read data')
ggsave('plot/quant_summary_redundancy_bin_boxplot.pdf', width = 8, height = 8)

quant_summary_cpm_redun %>%
  filter(quant %in% c('onts', 'map'),
         # redundantAsm > 1#,
         #assembler == 'rattle', truecpm_bins == '(10,Inf]'
  ) %>%
  ggplot(aes(x = truecount_quantile_range, y = redundantAsm, color = assembler)) +
  geom_boxplot() +
  scale_y_log10() +
  scale_color_manual(values = cols) +
  xlab('Quantiles on simulated counts long read data')
ggsave('plot/quant_summary_redundancy_quantile_boxplot.pdf', width = 8, height = 8)

lcpm.cor <- data.frame(gene_cor = sapply(quant_summary, function(x) cor(x$lcpm_gene[,1:2] )[1,2], simplify = T), 
                       tx_cor = sapply(quant_summary, function(x) cor(x$lcpm_tx[,1:2] )[1,2], simplify = T),
                       method = names(quant_summary),
                       gene_counts_ref = sapply(quant_summary, function(x) sum(x$lcpm_gene$true_counts), simplify = T),
                       gene_counts_asm = sapply(quant_summary, function(x) sum(x$lcpm_gene$counts), simplify = T),
                       tx_counts_ref = sapply(quant_summary, function(x) sum(x$lcpm_tx$true_counts), simplify = T),
                       tx_counts_asm = sapply(quant_summary, function(x) sum(x$lcpm_tx$counts.max), simplify = T)) %>%
  separate(method, c('col1','col2','salmon','assembler','quant'), remove = F) %>%
  mutate(assembler = factor(assembler, levels = c('bambu','rattle','rnabloom2','isonform', 'trinity')))

write.csv(lcpm.cor, 'plot/cor.csv')

lcpm.cor %>%
  filter(quant %in% c('onts', 'map') | assembler == 'bambu') %>%
  ggplot(aes(x = assembler, y = gene_cor, color = assembler)) +
  geom_point(aes(#shape = quant, 
    size = 1)) + 
  geom_line(aes(group = assembler)) + 
  ylim(c(0,1)) + 
  scale_color_manual(values = cols)
ggsave('plot/quant_summary_gene_cor_lineplot.pdf', width = 4, height = 3)

lcpm.cor %>%
  filter(quant %in% c('onts', 'map') | assembler == 'bambu') %>%
  ggplot(aes(x = assembler, y = tx_cor, color = assembler)) +
  geom_point(aes(#shape = quant, 
                 size = 1)) + 
  geom_line(aes(group = assembler)) + 
  ylim(c(0,1)) + 
  scale_color_manual(values = cols)
ggsave('plot/quant_summary_tx_cor_lineplot.pdf', width = 4, height = 3)


for (i in names(quant_summary)) {
  
  ggsave(paste0('plot/',i,'_quant_summary1.pdf'), plot = quant_summary[[i]]$plot1, width = 8, height = 16)
  ggsave(paste0('plot/',i,'_quant_summary2.pdf'), plot = quant_summary[[i]]$plot2, width = 8, height = 8)
}

saveRDS(quant_summary, 'plot/quant_summary.rds')
saveRDS(sqanti_summary, 'plot/sqanti_summary.rds')
