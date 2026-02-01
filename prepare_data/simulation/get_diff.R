
library(tidyverse)

baselineAbundance <- readRDS("baselineAbundance.rds") * 1e6

table <- read.table('sqanti-sim_index.tsv', header = T, sep = '\t')

table$requested_counts <- baselineAbundance[table$transcript_id]
table$requested_tpm <- table$requested_counts/(sum(table$requested_counts)) * 1e6

set.seed(1123)

geneid <- table %>% 
  filter(requested_counts >= 5) %>%
  count(gene_id) %>% 
  filter(n>=2 & n<=10) %>% 
  pull(gene_id) %>%
  sample(., 3000, replace = F)

## DGE: all expression double
dgeid <- geneid[1:1000]

## DTU: swap 2 of the tx from same gene
dtuid <- geneid[1001:2000]

## DGE&DTU: select 1 tx double
comid <- geneid[2001:3000]
combotxid <- table %>% 
  filter(gene_id %in% comid & requested_counts >= 5) %>%
  group_by(gene_id) %>%
  slice_sample(n=1) %>% pull(transcript_id)

df <- table

# DGE
df <- df %>% mutate(requested_counts2 = ifelse(gene_id %in% dgeid[1:500], 2*requested_counts, 
                                               ifelse(gene_id %in% dgeid[501:1000], 0.5*requested_counts, requested_counts)))
# DGE&DTU
df <- df %>% mutate(requested_counts2 = ifelse(transcript_id %in% combotxid[1:500], 2*requested_counts, 
                                               ifelse(transcript_id %in% combotxid[501:1000], 0.5*requested_counts, requested_counts2)))
# DTU
# Randomly select two rows from each group
selected_rows <- df %>% 
  filter(gene_id %in% dtuid & requested_counts >= 5)  %>% 
  group_by(gene_id) %>%
  slice_sample(n = 2)

# Swap values in a specific column
swapped_rows <- selected_rows %>% 
  group_by(gene_id) %>%
  mutate(temp_id = row_number(),
         requested_counts2 = case_when(
           temp_id == 1 ~ requested_counts2[2],
           temp_id == 2 ~ requested_counts2[1],
           TRUE ~ requested_counts2
         )) %>%
  select(-temp_id)

# Remove the originally selected rows from df
df_remaining <- anti_join(df, selected_rows, by = c("transcript_id", "requested_counts2")) 

# Combine the swapped rows back with the remaining DataFrame
final_df <- bind_rows(df_remaining, swapped_rows) 

total_reads <- sum(final_df$requested_counts2)
final_table <- final_df %>% mutate(requested_counts = requested_counts2,
                                   requested_tpm = requested_counts2/total_reads*1e6) %>%
  select(c(1:16))

## reoder final_table
matches <- match(table$transcript_id, final_table$transcript_id)
final_table <- final_table[matches,]

write.table(final_table, 'sqanti-sim_index_de.tsv',
            quote = F, row.names = F, sep = '\t')
write.table(table, 'sqanti-sim_index_ctrl.tsv',
            quote = F, row.names = F, sep = '\t')
#save.image()

dtutxid <- selected_rows$transcript_id

saveRDS(dtuid, 'dtuid.rds') ## only isoform switching, 1000
saveRDS(dgeid, 'dgeid.rds') ## only gene doubling, 1000
saveRDS(comid, 'comid.rds') ## only 1 tx doubling from the gene, 1000
saveRDS(combotxid, 'combotxid.rds') # selected tx with DTE, 1000
saveRDS(dtutxid, 'dtutxid.rds') # selected tx with switching, 2000

gtf <- table[,1:2]

true_dte <- c(dtutxid, # isoform switching
              combotxid, # selected tx with change
              gtf %>% filter(gene_id %in% dgeid) %>% pull(transcript_id) # all tx in DGE
              ) %>%
  unique() # 6933 transcripts
true_dge <- c(dgeid, comid) # 2000 genes

true_dtutx <- c(dtutxid,
                gtf %>% filter(gene_id %in% comid) %>% pull(transcript_id)) # 5927 transcripts
true_dtugene <- c(dtuid, 
                  comid) # 2000 genes

saveRDS(true_dte, 'true_dte.rds')
saveRDS(true_dge, 'true_dge.rds')
saveRDS(true_dtutx, 'true_dtutx.rds')
saveRDS(true_dtugene, 'true_dtugene.rds')


