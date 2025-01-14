get_stats <- function(x) {
  # Author: Feng Yan yan.a@wehi.edu.a
  # Date: 30/08/2023
  
  ## get stats from transrate output for summary and sequence level
  
  ## number of tx
  transrate <-
    read.table(paste0(x, '/assemblies.csv'),
               header = T,
               sep = ',')
  
  ## length and gc of tx distribution
  contigpath <-
    list.files(x,
               pattern = 'contigs.csv',
               recursive = T,
               full.names = T)
  
  seqstats <-
    read.csv(contigpath)
  
  return(list(summary = transrate, seqstats = seqstats))
  
}
