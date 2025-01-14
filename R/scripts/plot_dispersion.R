plot_dispersion <- function(x, minCount = 3) {
  require(data.table)
  require(tidyverse)
  require(magrittr)
  
  keep <- filterByExpr(x, group=x$samples$group, min.count = minCount)
  x <- x[keep, , keep.lib.sizes=FALSE]
  
  dt.mao.plot <- as.data.table(x$genes %>%
                                 add_count(GENEID, name = 'NTranscriptPerGene'))
  
  dt.mao.plot$TranscriptID <- rownames(x)
  
  dt.mao.plot$Overdispersion <- 
    salmon_tx$annotation$Overdispersion[match(dt.mao.plot$TranscriptID, rownames(salmon_tx$annotation))]
  
  dt.mao.plot[,NTranscriptPerGeneTrunc := 
                ifelse(NTranscriptPerGene<10,NTranscriptPerGene,paste0('>=10'))]
  
  dt.mao.plot$NTranscriptPerGeneTrunc %<>% factor(levels = paste0(c(1:10,'>=10')))
  
  # Number of transcripts from single-transcript genes
  # dt.mao.plot[NTranscriptPerGene == 1,
  #             .(.N,sum(NTranscriptPerGene == 1 & Overdispersion>(1/0.9))/sum(NTranscriptPerGene == 1))]
  # 
  # dt.mao.plot[NTranscriptPerGene > 1,
  #             .(.N,sum(NTranscriptPerGene > 1 & Overdispersion>(1/0.9))/sum(NTranscriptPerGene > 1))]
  # 
  # dt.mao.plot[NTranscriptPerGene >= 10,
  #             .(NGeneID = length(unique(GENEID)),
  #               NTranscriptID = .N,mean(Overdispersion))]
  
  dt.mao.plot.long <- 
    melt(dt.mao.plot,
         id.vars = c('TranscriptID','NTranscriptPerGeneTrunc'),
         measure.vars = c('Overdispersion'),
         variable.name = 'Type',value.name = 'Overdispersion')
  
  plot.mao <-
    ggplot(data = dt.mao.plot.long,
           aes(x = NTranscriptPerGeneTrunc,y = log10(Overdispersion))) +
    geom_smooth(aes(group = Type,color = Type),se = FALSE,span = 0.8,method = 'loess',show.legend = FALSE) +
    geom_boxplot(aes(fill = Type),outlier.alpha = 0.2,col = 'black',alpha = 0.5) +
    labs(x = 'Number of transcripts per gene', y = 'Overdispersion (log10 scale)') +
    scale_y_continuous(limits = c(0,2.5),breaks = seq(0,3,0.5)) +
    # scale_fill_manual(values = c('Illumina short reads' = 'red','ONT long reads' = 'blue')) +
    # scale_color_manual(values = c('Illumina short reads' = 'red','ONT long reads' = 'blue')) +
    theme_bw(base_size = 8,base_family = 'sans') +
    theme(panel.grid = element_blank(),
          axis.text = element_text(colour = 'black', size = 8),
          legend.title = element_blank(),
          legend.position = c(0.175, 0.925),
          legend.text = element_text(colour = 'black', size = 8),
          legend.background = element_rect(fill = alpha("white", 0)))
  
  return(plot.mao)
  
}
