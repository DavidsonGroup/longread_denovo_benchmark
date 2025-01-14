
upset_sqanti <- function(list, featuretype = 'transcript', 
                         subcat = 'all', 
                         names = c('10m','1m','2m','5m')) {
  require(fastmatch)
  require(UpSetR)
  
  ## subcat can be 
  ## "antisense" "full-splice_match" "fusion" "genic" "genic_intron" "incomplete-splice_match" "intergenic" "novel_in_catalog""novel_not_in_catalog"   
  
  if ('all' %in% subcat) {
    listsub <- list
  } else {
    listsub <- lapply(list, function(x) x %>% filter(structural_category %in% subcat))
  }
  
  if (featuretype == 'transcript') {
    fid <- rbindlist(listsub)$associated_transcript %>% unique()
    fid <- fid[!(fid == 'novel')] # remove 'novel'
    
    # Get a list of associated_transcript/gene columns
    feature_lists <- lapply(listsub, function(df) df$associated_transcript)
    
  } else if (featuretype == 'gene') {
    fid <- rbindlist(listsub)$associated_gene %>% unique()
    fid <- fid[!grepl('novel', fid)] %>% unique()
    
    # Get a list of associated_transcript/gene columns
    feature_lists <- lapply(listsub, function(df) df$associated_gene)
  }
  
  # Apply the function to each element in fid
  feature_count <- t(vapply(fid, function(x) {
    check_in_list(x, feature_lists)
  }, FUN.VALUE = logical(length(feature_lists))))
  
  # convert to numeric matrix
  feature_count_df <- feature_count %>% apply(., 2, as.numeric) %>% data.frame() 
  colnames(feature_count_df) <- names
  rownames(feature_count_df) <- rownames(feature_count)
  upset(feature_count_df, nsets = 4)
  return(feature_count_df)
}


# Function to check if an element is in a list of vectors
check_in_list <- function(x, list_of_vectors) {
  vapply(list_of_vectors, function(vector) {
    any(fmatch(x, vector, nomatch = 0))
  }, FUN.VALUE = logical(1))
}





