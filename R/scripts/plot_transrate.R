plot_transrate <- function(df, depth = TRUE) {
  
  ## df is a data frame containing columns: length, prop_gc, reference_coverage, assembler, depth
  ## if only 1 dataset with no different depth, set depth = FALSE to suppress alpha

  alpha <- c('2' = 0.1, '5' = 0.5, '10' = 1)
    
  df %>% 
    ggplot(aes(x=assembler, y=length, fill = assembler, alpha = factor(depth))) +
    geom_boxplot(position = position_dodge2(preserve = "single")) +
    theme(axis.text.x=element_text(angle=45, hjust=1)) +
    stat_summary(fun.y = function(x) log10(mean(10^x)), geom = "point", 
                 position = position_dodge(width = 0.75)) +  # Add points to plot ,  col = "red"
    stat_summary(fun.y = function(x) log10(mean(10^x)), geom = "text", # col = "red",     # Add text to plot
                 vjust = 1.5, aes(label = round(..y.., digits = 0)),
                 position = position_dodge(width = 0.75)) +
    scale_y_log10() +
    scale_fill_manual(values = cols) + 
    scale_alpha_manual(values = alpha) +
    ylab('Transcript length (bp)')
  ggsave('plot/transcript_length.pdf', width = 6, height = 4)
  
  df %>% 
    ggplot(aes(x=assembler, y=prop_gc, fill = assembler, alpha = factor(depth))) +
    geom_boxplot(position = position_dodge2(preserve = "single")) + 
    theme(axis.text.x=element_text(angle=45, hjust=1)) +
    stat_summary(fun = mean, geom = "point",
                 position = position_dodge(width = 0.75)) +  # Add points to plot ,  col = "red"
    stat_summary(fun = mean, geom = "text", # col = "red",     # Add text to plot
                 vjust = 1.5, aes(label = round(..y.., digits = 2)),
                 position = position_dodge(width = 0.75)) +
    scale_fill_manual(values = cols) + 
    scale_alpha_manual(values = alpha) +
    ylab('GC (%)')
  ggsave('plot/transcript_gc.pdf', width = 6, height = 4)
  
  df %>% 
    ggplot(aes(x=assembler, y=reference_coverage, fill = assembler, alpha = factor(depth))) +
    geom_boxplot(position = position_dodge2(preserve = "single")) + 
    theme(axis.text.x=element_text(angle=45, hjust=1)) +
    stat_summary(fun = mean, geom = "point",
                 position = position_dodge(width = 0.75)) +  # Add points to plot ,  col = "red"
    stat_summary(fun = mean, geom = "text", # col = "red",     # Add text to plot
                 vjust = 1.5, aes(label = round(..y.., digits = 2)),
                 position = position_dodge(width = 0.75)) +
    scale_fill_manual(values = cols) + 
    scale_alpha_manual(values = alpha) +
    ylab('Reference coverage (%))')
  ggsave('plot/transcript_refcov.pdf', width = 6, height = 4)
  
  if (!(depth)) {
    df$depth <- 1 # for plotting bar
  }
  
  df %>% 
    mutate(refcov_cat = cut(reference_coverage, breaks = c(0,0.25,0.5,0.75,1))) %>%
    ggplot(aes(x=depth, fill = refcov_cat)) +
    geom_bar() + 
    theme(axis.text.x=element_text(angle=45, hjust=1)) +
    stat_count(geom = "text", colour = "white", size = 3.5,
               aes(label = ..count..), position=position_stack(vjust=0.5)) +
    facet_grid(.~assembler, scales = 'free_x', space = "free") +
    ylab('Number of transcript')
  ggsave('plot/transcript_count_refcov.pdf', width = 8, height = 4)
  
  df %>% 
    mutate(refcov_cat = cut(reference_coverage, breaks = c(0,0.25,0.5,0.75,1))) %>%
    ggplot(aes(x=depth, fill = refcov_cat)) +
    geom_bar(position="fill") +
    theme(axis.text.x=element_text(angle=45, hjust=1)) +
    facet_grid(.~assembler, scales = 'free_x', space = "free") +
    ylab('Percentage of transcript')
  ggsave('plot/transcript_pct_refcov.pdf', width = 8, height = 4)
  
  df %>%
    mutate(length_cat = cut(length, breaks = c(0,500,1000, 3000, Inf))) %>%
    ggplot(aes(x=depth, fill = length_cat)) +
    geom_bar() + 
    theme(axis.text.x=element_text(angle=45, hjust=1)) +
    stat_count(geom = "text", colour = "white", size = 3.5,
               aes(label = ..count..), position=position_stack(vjust=0.5)) +
    facet_grid(.~assembler, scales = 'free_x', space = "free") +
    ylab('Number of transcript')
  ggsave('plot/transcript_count_length.pdf', width = 8, height = 4)
  
  df %>%
    mutate(length_cat = cut(length, breaks = c(0,500,1000, 3000, Inf))) %>%
    ggplot(aes(x=depth, fill = length_cat)) +
    geom_bar(position="fill") +
    theme(axis.text.x=element_text(angle=45, hjust=1)) +
    facet_grid(.~assembler, scales = 'free_x', space = "free") +
    ylab('Percentage of transcript')
  ggsave('plot/transcript_pct_length.pdf', width = 8, height = 4)
}
