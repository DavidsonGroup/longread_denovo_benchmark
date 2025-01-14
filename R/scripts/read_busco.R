read_busco <- function (result_dir = NULL) 
{
  ## from https://github.com/almeidasilvaf/cogeqc
  
  file <- dir(result_dir)
  file <- file[grepl("txt", file)]
  if (length(file) != 1) {
    message("More than 1 BUSCO summary file found. Using only the first.")
    file <- file[1]
  }
  full_path <- paste0(result_dir, "/", file)
  cl <- c("Complete_SC", "Complete_duplicate", "Fragmented", 
          "Missing")
  if (startsWith(file, "short")) {
    lines <- readLines(paste0(result_dir, "/", file))
    lineage <- lines[grepl("dataset", lines)]
    lineage <- gsub(".*dataset is: | \\(Creation.*", "", 
                    lineage)
    complete_sc <- lines[grepl("single-copy BUSCOs", lines)]
    complete_sc <- gsub("\\tComplete.*|\\t", "", complete_sc)
    complete_dup <- lines[grepl("duplicated BUSCOs", lines)]
    complete_dup <- gsub("\\tComplete.*|\\t", "", complete_dup)
    fragmented <- lines[grepl("Fragmented BUSCOs", lines)]
    fragmented <- gsub("\\tFragmented.*|\\t", "", fragmented)
    missing <- lines[grepl("Missing BUSCOs", lines)]
    missing <- gsub("\\tMissing.*|\\t", "", missing)
    final_df <- data.frame(Class = factor(cl, levels = cl), 
                           Frequency = as.numeric(c(complete_sc, complete_dup, 
                                                    fragmented, missing)), Lineage = lineage)
  }
  else if (startsWith(file, "batch")) {
    df <- read.csv(full_path, sep = "\t", header = FALSE, 
                   skip = 1)
    df <- df[, c(1, 2, 4:7)]
    colnames(df) <- c("File", "Lineage", "Complete_SC", 
                      "Complete_duplicate", "Fragmented", "Missing")
    final_df <- reshape2::melt(df, id = c("File", "Lineage"))
    colnames(final_df) <- c("File", "Lineage", "Class", 
                            "Frequency")
    final_df <- final_df[, c("Class", "Frequency", "Lineage", 
                             "File")]
    final_df$Class <- factor(final_df$Class, levels = cl)
  }
  else {
    stop("Wrong results directory. Could not find BUSCO summary file.")
  }
  return(final_df)
}
