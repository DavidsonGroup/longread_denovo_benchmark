get_gene_tx <- function(tx_table, sqanti_class, method = 'rattle', extra = NULL) {
  
  require(data.table)
  
  if (method == 'rattle') {
    
    ## use fx2tab results
    ## command: seqkit fx2tab -nlgGH transcriptome.fasta > fx2tab.txt
    
    table <- read.table(tx_table, header = F, sep = '\t') %>% 
    separate(V1, c('isoform','geneid','origin_txid','counts','label'), sep = ' ') %>% 
      dplyr::mutate(counts = str_remove(counts, 'total_reads=') %>% as.numeric()) 

  } else if (method %in% c('corset','rnabloom2','bambu')) {
    
    ## result after run corset on rnabloom2 assembly

    table <- read.table(tx_table, 
                         sep = '\t', header = F) %>%
      dplyr::rename(isoform = V1, geneid = V2)
    
  } else if (method == 'trinity') {
    
    table <- read.table(tx_table, 
                        sep = '\t', header = F) %>%
      dplyr::rename(isoform = V2, geneid = V1)
    
  } else if (method == 'isonclust') {
    
    ## use fx2tab results
    ## command: seqkit fx2tab -nlgGH transcriptome.fasta > fx2tab.txt
    
    table <- read.table(tx_table, header = F, sep = '\t') %>% 
      separate(V1, c('geneid','batchid','id'), sep = '_', remove = F) %>%  
      dplyr::rename(isoform = V1)
  }
  
  sqanti_summary <- read.table(sqanti_class, sep = '\t', header = T) %>%
    select(c(1:8,15:17,25,29,30,41,46)) %>%
    mutate(original = ifelse(str_detect(isoform, '_dup'), str_remove(isoform, '_dup.*'), isoform))
  
  gene <- left_join(table, sqanti_summary, by = 'isoform') %>%
    dplyr::filter(!str_detect(isoform, '_dup')) %>%
    mutate(geneid = ifelse(is.na(geneid), 
                           paste('rowid', row_number(), sep = '_'),  # for corset results, sometimes not all isforms are spit by corset/minimap2
                           geneid))
  
  if (rlang::is_empty(extra)) {
    
    gene <- gene %>% mutate(mapq=NA, nm=NA)
  
  } else {
    
    meta <- fread(extra) %>%
      mutate(V3 = str_remove(V3, 'NM:i:') %>% as.numeric())
    colnames(meta) <- c('isoform','mapq','nm') 
    
    gene <- left_join(gene, meta, by = 'isoform')
  }
  
  return(gene)
}
  