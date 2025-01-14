quantile_overlap <- function(x, y, title, requested_counts) {
  
  # x is from sqanti3 summary (use knwon gtf as ref)
  # y is salmon quant use assembled fasta
  
  # i=8
  # x <- sqanti_summary_tmp[[i]]
  # y <- salmon_count[[i]]
  
  # requested counts contains 4 columns: requested_counts, sim_counts, transcript_id, length, gene_id, dtu, dte
  
  ## reference transcript centric metrics
  ## join true count and assemble count
  true.counts <- data.frame(#requested_counts = requested_counts$requested_counts,
                            true_counts = requested_counts$sim_counts, 
                            txid = requested_counts$transcript_id, 
                            #txlength = requested_counts$length,
                            dte = requested_counts$dte,
                            dtu = requested_counts$dtu,
                            truecount_quantile = requested_counts$truecount_quantile,
                            truecount_quantile_range = requested_counts$truecount_quantile_range) 
  
  join.counts <- true.counts %>%
    # dplyr::filter(requested_counts != 0) %>% 
    left_join(left_join(x, y, by='isoform') %>%
                dplyr::filter(!str_detect(isoform, '_dup')), ## remove _dup in the sqanti output, they have NA counts,
              by = c('txid' = 'associated_transcript')) %>% 
    add_count(txid) %>%
    group_by(txid) %>%
    summarise(true_counts = mean(true_counts), ## if duplicated txid, bambu_count should be same
              #requested_counts = mean(requested_counts),
              counts.max = if(all(is.na(counts))) NA else max(counts, na.rm = TRUE), 
              counts.sum = if(all(is.na(counts))) NA else sum(counts, na.rm = TRUE), ## add up assembled tx matched to same txid, when some assembled tx is NA, add with na.rm, when all NA, means not captures, so leave as NA ,
              #txlength = mean(txlength), 
              dtu = mean(dtu), dte = mean(dte), 
              truecount_quantile = unique(truecount_quantile), truecount_quantile_range = unique(truecount_quantile_range),
              redundantAsm = mean(n)) %>%  ## if duplicated txid, n should be same, number of assembled tx mapping to same ref
    mutate(isassemble = ifelse(is.na(counts.sum), 'not assembled', ifelse(counts.sum == 0, 'not quantified', 'quantified')),
           redundantAsm = ifelse(is.na(counts.sum), 0, redundantAsm)) %>% # change not assembled to 0
    mutate(isassemble = factor(isassemble, levels = c('not assembled','not quantified','quantified'))) #%>%
    #left_join(requested_counts[,c(1,2,6)], by = c('txid' = 'transcript_id'))
  
  # change NA to 0 to calculate cpm without changing lib size
  join.counts <- join.counts %>% 
    replace_na(list(counts.sum = 0, counts.max = 0)) %>%  
    data.frame()
  
  rownames(join.counts) <- join.counts$txid
  
  # dge <- DGEList(counts = join.counts[,2:4], 
  #                genes = join.counts[, c(1,5:ncol(join.counts))])
  # 
  # dge <- calcNormFactors(dge)
  # 
  # cpm <- cpm(dge, log = T, prior.count = 1)
  
  libsize <- colSums(join.counts[, c('true_counts','counts.max','counts.sum')])
  print(libsize)
  
  cpm <- apply(join.counts[, c('true_counts','counts.max','counts.sum')], 2, 
               function(x) (x/sum(x))*1000000)
  cpm <- log2(cpm + 1)

  print(dim(cpm))
  # library size
  # print(dge$samples)
  
  cpm <- data.frame(cpm) %>%
    mutate(txid = rownames(cpm), 
           truecpm_bins = cut(true_counts, breaks = c(0:10, Inf))) %>% 
    dplyr::rename(true_cpm = true_counts, 
           assemble_cpm.sum = counts.sum,
           assemble_cpm.max = counts.max) %>% 
    left_join(join.counts, by = 'txid')
  rownames(cpm) <- cpm$txid
  
  p1 <- cpm %>% 
    dplyr::count(truecpm_bins, isassemble) %>%
    ggplot(aes(x = truecpm_bins, y = n, fill = isassemble)) +
    geom_bar(stat = 'identity', position = 'fill') +
    theme(legend.position = "bottom") +
    guides(fill=guide_legend(nrow=2, byrow=TRUE)) +
    ggtitle(title)
  
  p2 <- cpm %>% 
    dplyr::count(truecpm_bins, isassemble) %>%
    ggplot(aes(x = truecpm_bins, y = n, fill = isassemble)) +
    geom_bar(stat = 'identity', position = 'stack') +
    theme(legend.position = "bottom") +
    guides(fill=guide_legend(nrow=2, byrow=TRUE)) +
    ggtitle(title)
  
  cor <- cpm %>%
    dplyr::filter(isassemble != 'not assembled') %>%
    group_by(truecpm_bins) %>%
    summarise(correlation = cor(true_counts, counts.max))
  
  # cor <- cpm %>% 
  #   dplyr::filter(isassemble != 'not assembled') %>%
  #   dplyr::select(1:3) %>%
  #   as.matrix() %>%
  #   cor()
  
  p3 <- cor %>%
    mutate(data = title) %>%
    ggplot(aes(x = truecpm_bins, y = correlation, group = data)) +
    geom_point() +
    geom_line() +
    ylim(c(-0.2,1)) +
    theme(legend.position = "bottom") +
    guides(fill=guide_legend(nrow=2, byrow=TRUE)) +
    ggtitle(title)
  
  p4 <- ggscatter(cpm %>% 
                    dplyr::filter(isassemble != 'not assembled'),
                  x = 'true_cpm', y = 'assemble_cpm.max',
                  color = 'gray',
                  conf.int = F,
                  cor.coef = T) + 
    geom_abline(linetype = 2) +
    theme(legend.position = "bottom") +
    guides(fill=guide_legend(nrow=2, byrow=TRUE)) +
    ggtitle(paste(title, 'Pearson correlation'))
  
  p4.1 <- ggscatter(cpm %>% 
                    dplyr::filter(isassemble != 'not assembled'),
                  x = 'true_cpm', y = 'assemble_cpm.max',
                  color = 'dte', facet.by = 'dte',
                  conf.int = F,
                  cor.coef = T) + 
    geom_abline(linetype = 2) +
    theme(legend.position = "bottom") +
    guides(fill=guide_legend(nrow=2, byrow=TRUE)) +
    ggtitle(paste(title, 'Pearson correlation DTE'))
  
  p4.2 <- ggscatter(cpm %>% 
                    dplyr::filter(isassemble != 'not assembled'),
                  x = 'true_cpm', y = 'assemble_cpm.max',
                  color = 'dtu', facet.by = 'dtu',
                  conf.int = F,
                  cor.coef = T) + 
    geom_abline(linetype = 2) +
    theme(legend.position = "bottom") +
    guides(fill=guide_legend(nrow=2, byrow=TRUE)) +
    ggtitle(paste(title, 'Pearson correlation DTU'))
  
  p5 <- cpm %>% 
    dplyr::mutate(redundantAsm_bin = cut(redundantAsm, breaks = c(0,1,2,5,10,Inf), right = F)) %>%
    dplyr::count(truecpm_bins, redundantAsm_bin) %>%
    ggplot(aes(x = truecpm_bins, y = n, fill = redundantAsm_bin)) +
    geom_bar(stat = 'identity', position = 'fill') +
    theme(legend.position = "bottom") +
    guides(fill=guide_legend(nrow=2, byrow=TRUE)) +
    ggtitle(title)
  
  lcpm_tx <- cpm
  
  ## assembly centric metrics, divide assembled tx 
  tmp <- x %>% add_count(original) %>%
    left_join(y, by='isoform') %>%
    mutate(ischimeric = n>1,
           isprimary = !(str_detect(isoform, '_dup'))) %>% 
    left_join(true.counts,
              by = c('associated_transcript'='txid')) %>%
    mutate(isnovel = is.na(true_counts), # same if use associated_transcript=='novel'
           istrueexp = ifelse(is.na(true_counts), 'novel', 
                              ifelse(true_counts == 0, 'not quantified','quantified')), # 0 means not quantified, NA means novel
           isbambunovel = str_detect(associated_transcript, 'BambuTx'),
           ischrM = chrom == 'chrM') %>%
    dplyr::mutate(row_id = row_number()) %>%
    dplyr::mutate(asso_tx = ifelse(associated_transcript=='novel', 
                                   paste(associated_transcript, row_id, sep = '_'), 
                                   associated_transcript)) %>% 
    add_count(asso_tx) 
  
  p6 <- tmp %>% ggplot(aes(x = isprimary)) +
    geom_bar(aes(fill = istrueexp),
             position = 'stack') + 
    ylab('Number of Tx') +
    theme(legend.position = "bottom") +
    guides(fill=guide_legend(nrow=2, byrow=TRUE)) +
    ggtitle(title)
  
  p7 <- tmp %>% ggplot(aes(x = isprimary)) +
    geom_bar(aes(fill = istrueexp),
             position = 'fill') +
    theme(legend.position = "bottom") +
    ylab('Percentage of Tx') +
    guides(fill=guide_legend(nrow=2, byrow=TRUE)) +
    ggtitle(title)
  
  p8 <- tmp %>% ggplot(aes(x = isprimary)) +
    geom_bar(aes(fill = ischrM),
             position = 'fill') +
    theme(legend.position = "bottom") +
    ylab('Percentage of Tx') +
    guides(fill=guide_legend(nrow=2, byrow=TRUE)) +
    ggtitle(title)
  
  p9 <- tmp %>% 
    dplyr::mutate(tx_redundancy = cut(nn, breaks = c(0,1,2,5,10,Inf), right = F)) %>%
    ggplot(aes(isprimary)) + 
    geom_bar(aes(fill = tx_redundancy),
             position = 'stack') +
    theme(legend.position = "bottom") +
    ylab('Number of Tx') +
    guides(fill=guide_legend(nrow=2, byrow=TRUE)) +
    ggtitle(title)
  
  p_2 <- plot_grid(p6, p7, p8, p9, nrow = 2,
                 labels = c('A', 'B', 'C', 'D'), label_size = 8)
  
  ## reference gene centric metrics
  ## join true count and assemble count
  true.counts <- data.frame(#requested_counts = requested_counts$requested_counts,
                            true_counts = requested_counts$sim_counts, 
                            gene_id = requested_counts$gene_id) %>%
    filter(!str_detect(gene_id, '^BambuGene')) %>%
    group_by(gene_id) %>%
    summarise_all(sum) %>%
    ungroup() %>%
    mutate(truecount_quantile = ntile(true_counts, 10)) %>% # add quantiles for counts, equal number of elements in each quantile
    group_by(truecount_quantile) %>%
    mutate(truecount_quantile_range = paste0("[", round(min(true_counts), 2), "-", round(max(true_counts), 2), "]")) %>% # some ties are bread into each bin
    ungroup() 
  
  join.counts <- true.counts %>%
    # dplyr::filter(requested_counts != 0) %>% 
    left_join(left_join(x, y, by='isoform') %>%
                dplyr::filter(!str_detect(isoform, '_dup')) %>%
                group_by(associated_gene) %>%
                summarise(counts = sum(counts, na.rm = T)), ## remove _dup in the sqanti output, they have NA counts,
              by = c('gene_id' = 'associated_gene')) %>%
    mutate(isassemble = !(is.na(counts)))
  
  # change NA to 0 to calculate cpm without changeing lib size
  join.counts <- join.counts %>% 
    replace_na(list(counts=0)) %>%  
    data.frame()
  
  rownames(join.counts) <- join.counts$gene_id
  
  # dge <- DGEList(counts = join.counts[,2:4], 
  #                genes = join.counts[, 1])
  # 
  # dge <- calcNormFactors(dge)
  # 
  # cpm <- cpm(dge, log = T, prior.count = 1)
  
  libsize <- colSums(join.counts[, c('true_counts', 'counts')])
  print(libsize)
  
  cpm <- apply(join.counts[,c(c('true_counts', 'counts'))], 2, 
               function(x) (x/sum(x))*1000000)
  cpm <- log2(cpm + 1)
  
  print(dim(cpm))
  
  # library size
  # print(dge$samples)
  
  cpm <- data.frame(cpm) %>%
    mutate(gene_id = rownames(cpm), 
           truecpm_bins = cut(true_counts, breaks = c(0:10, Inf))) %>% 
    dplyr::rename(true_cpm = true_counts, 
                  assemble_cpm = counts) %>% 
    left_join(join.counts, by = 'gene_id')
  rownames(cpm) <- cpm$gene_id
  
  lcpm_gene <- cpm
  # cor <- cpm %>% 
  #   dplyr::select(1:3) %>%
  #   as.matrix() %>%
  #   cor()
  
  p.gene <- ggscatter(cpm,
                 x = 'true_cpm', y = 'assemble_cpm',
                 color = 'gray',
                 conf.int = F,
                 cor.coef = T) + 
    geom_abline(linetype = 2) +
    theme(legend.position = "bottom") +
    guides(fill=guide_legend(nrow=2, byrow=TRUE)) +
    ggtitle(title)
  
  p_1 <- plot_grid(p1, p2, 
                   p4, p.gene,
                   p4.1, p4.2, 
                   p3, p5, 
                   nrow = 4,
                   labels = c('A', 'B', 'C', 'D','E','F','G', 'H'), label_size = 8)
  
  return(list(lcpm_gene = lcpm_gene, lcpm_tx = lcpm_tx, plot1 = p_1, plot2 = p_2))
  
}
