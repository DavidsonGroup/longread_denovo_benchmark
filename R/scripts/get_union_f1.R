get_union_f1 <- function(sqanti_anno, truth = tx_union, keep.all.truth = T, remove.nomatch = T) {
  
  if (keep.all.truth) {
    
    if (remove.nomatch) {
      
      subset <- truth %>%
        left_join(sqanti_anno, by = c('associated_transcript', 'associated_gene') ) %>%
        filter(!is.na(isoform))
    
      } else if (!remove.nomatch) {
      
      subset <- truth %>%
        left_join(sqanti_anno, by = c('associated_transcript', 'associated_gene') ) 
    
      }
    
  } else if (!keep.all.truth) {
    # this is to match the previous analysis from eval_cluster.R which use all de novo transcritps that matches to know gene (FSM/ISM, NNC/NIC, genic)
    # keep de novo transcripts that match to the known set of transcripts
    # the number of de novo can still be different
    # mutiple de novo can be assigned to same known tx
    
    subset <- sqanti_anno %>% 
      filter(associated_transcript %in% truth$associated_transcript & gene_category == 'known_gene')
    
  }
  
  # eval_report_clusters(subset$geneid, 
  #                      subset$associated_gene) %>% 
  #   unlist()  
  
  ## get some more metrics
  tptable <- subset %>% dplyr::count(associated_gene, geneid) %>%
    dplyr::filter(n >1) %>% 
    dplyr::select(geneid, n)
  
  TP <- sum(choose(tptable$n, 2))
  # TP=sum(sapply(scl,function(x){ sum(choose(table(x),2)) }))
  TPFP=sum(choose(table(subset$geneid),2))
  TPFN=sum(choose(table(subset$associated_gene),2))
  #show(paste("True Positives:",TP,"TPNP:",TPFP,"TPFN:",TPFN))
  
  FP <- TPFP-TP
  FN <- TPFN-TP
  all <- choose(dim(subset)[1],2)
  TN <- all-FP-FN-TP
  P <- TP/TPFP
  R <- TP/TPFN
  
  res <- data.frame(number = dim(subset)[1],
                    FP=FP,
                    FN=FN,
                    TN=TN,
                    TP=TP,
                    Precision=P,
                    Recall=R,
                    #RAND=(TP+TN)/all,
                    Fscore=2*(R*P)/(P+R)
                    #FM=sqrt((TP/TPFP)*(TP/TPFN))
                    )
  
  return(res)
  # cat(paste("RAND:", RAND, '\nFM:', FM, '\nF1:', Fscore,
  #           '\nPrecision:', P, '\nRecall:', R,
  #           '\n'))
}
