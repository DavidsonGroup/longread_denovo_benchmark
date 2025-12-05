# all quant
# tiled sample level
setwd('~/lab_davidson/yan.a/Dong_2023/R_revise/')
all_quant <- readRDS("plot/all_quant.rds")

tx_counts_bambu <- readRDS('../bambu/merged_10m_dge/R/y.rds')
tx_counts_bambu <- tx_counts_bambu[rowSums(tx_counts_bambu$counts) != 0, ]
dim(tx_counts_bambu)

# for non sequin
txorder <- all_quant[[1]]$tx_exp_tiled %>%
  filter(!str_detect(tx, '^R|^Bambu')) %>% pull(tx) %>% unique()

truecpm_tx <- log2(tx_counts_bambu$counts + 1)[txorder, ] %>%
  as.matrix() %>%
  c()

names(truecpm_tx) <- rep(txorder, 6)

mat <- lapply(all_quant, function(x) {
  x$tx_exp_tiled %>%
    filter(!str_detect(tx, '^R|^Bambu')) %>%
    pull(assemble_cpm)
}) %>% Reduce(cbind, .) 
colnames(mat) <- names(all_quant)
df <- data.frame(mat) 

libsize <- lapply(all_quant, function(x) {
  x$tx_exp_tiled$counts.max
}) %>% Reduce(cbind, .) 

colnames(libsize) <- names(all_quant)

# rattle or rnabloom2 vs truth
# tiled expression level
plots <- lapply(names(df), function(colname) {
  ggplot(df, aes(x = .data[[colname]], y = truecpm_tx)) +
    geom_hex(bins = 50) +
    scale_fill_viridis_c(trans = "log10") +
    labs(
      x = colname,
      y = 'true_lcpm',
      title = paste("Correlation =", round(cor(df[[colname]], truecpm_tx), 3))
    ) +
    theme_bw()
})

library(patchwork)
wrap_plots(plots[16:17], ncol = 2)

library(ggplot2)
library(gridExtra) # Required for arranging multiple plots

chunk_size   <- 89305
plot_list    <- list()

## Index of the transcript you want to highlight
highlight_idx <- which(names(truecpm_tx) == "ENST00000199764.7")

for (i in 1:6) {
  # Define indices for the current chunk
  start_idx <- (i - 1) * chunk_size + 1
  end_idx   <- i * chunk_size
  # Be safe in case you run past the end
  end_idx   <- min(end_idx, length(truecpm_tx))
  
  idx_range <- start_idx:end_idx
  
  # Subset df and vector, and align them
  chunk_df        <- df[idx_range, ]
  chunk_df$true_lcpm <- truecpm_tx[idx_range]
  
  # Correlation for this chunk
  chunk_cor <- round(cor(chunk_df$rattle.10m.sec, chunk_df$true_lcpm), 3)
  
  # Base plot
  p <- ggplot(chunk_df, aes(x = rattle.10m.sec, y = true_lcpm)) +
    #geom_bin_2d() +
    geom_hex(bins = 50) +
    scale_fill_viridis_c(trans = "log10") +
    labs(
      x = "rattle",
      y = "true_lcpm",
      title = paste("RATTLE Sample", i, "| Cor =", chunk_cor)
    ) +
    xlim(c(0, 15)) +
    ylim(c(0, 15)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    theme_bw()
  
  # idx_in_chunk <- highlight_idx[highlight_idx >= start_idx & highlight_idx <= end_idx]
  # 
  # ## 2) If any, convert global index -> row number in this chunk
  # if (length(idx_in_chunk) > 0) {
  #   rows_in_chunk  <- idx_in_chunk - start_idx + 1
  #   highlight_rows <- chunk_df[rows_in_chunk, , drop = FALSE]
  #   
  #   ## 3) Add highlighted point(s)
  #   p <- p +
  #     geom_point(
  #       data   = highlight_rows,
  #       aes(x = rattle.10m.sec, y = true_lcpm),
  #       colour = "red",
  #       size   = 3
  #     )
  # }
  # 
  plot_list[[i]] <- p
}

# 3. Arrange and display all 6 plots in a grid (e.g., 2 rows of 3)
wrap_plots(grobs = plot_list, nrow = 2, ncol = 3)


chunk_size   <- 89305
plot_list    <- list()

## Index of the transcript you want to highlight
highlight_idx <- which(names(truecpm_tx) == "ENST00000199764.7")

for (i in 1:6) {
  # Define indices for the current chunk
  start_idx <- (i - 1) * chunk_size + 1
  end_idx   <- i * chunk_size
  # Be safe in case you run past the end
  end_idx   <- min(end_idx, length(truecpm_tx))
  
  idx_range <- start_idx:end_idx
  
  # Subset df and vector, and align them
  chunk_df        <- df[idx_range, ]
  chunk_df$true_lcpm <- truecpm_tx[idx_range]
  
  # Correlation for this chunk
  chunk_cor <- round(cor(chunk_df$rnabloom2.10m.sec, chunk_df$true_lcpm), 3)
  
  # Base plot
  p <- ggplot(chunk_df, aes(x = rnabloom2.10m.sec, y = true_lcpm)) +
    #geom_bin_2d() +
    geom_hex(bins = 50) +
    scale_fill_viridis_c(trans = "log10") +
    labs(
      x = "rnabloom2",
      y = "true_lcpm",
      title = paste("RNABloom2 Sample", i, "| Cor =", chunk_cor)
    ) +
    xlim(c(0, 15)) +
    ylim(c(0, 15)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    theme_bw()
  
  # idx_in_chunk <- highlight_idx[highlight_idx >= start_idx & highlight_idx <= end_idx]
  # 
  # ## 2) If any, convert global index -> row number in this chunk
  # if (length(idx_in_chunk) > 0) {
  #   rows_in_chunk  <- idx_in_chunk - start_idx + 1
  #   highlight_rows <- chunk_df[rows_in_chunk, , drop = FALSE]
  #   
  #   ## 3) Add highlighted point(s)
  #   p <- p +
  #     geom_point(
  #       data   = highlight_rows,
  #       aes(x = rattle.10m.sec, y = true_lcpm),
  #       colour = "red",
  #       size   = 3
  #     )
  # }
  
  plot_list[[i]] <- p
}

# 3. Arrange and display all 6 plots in a grid (e.g., 2 rows of 3)
wrap_plots(grobs = plot_list, nrow = 2, ncol = 3)
# 
# # only look at assembled
# # all quant 2
# all_quant2 <- readRDS("/vast/projects/lab_davidson/yan.a/dRNA_gencode/R/plot/all_quant2.rds")
# 
# mat2 <- lapply(all_quant2, function(x){
#   
#   txid <- unique(x$tx_exp_tiled %>%
#                    filter(!str_detect(tx, '^R')) %>%
#                    pull(tx))
#   
#   df_ass <- x$tx_exp_tiled %>%
#     filter(tx %in% txid) 
#   
#   stopifnot(df_ass$tx == rep(txid, 6))
#   
#   tx <- data.frame(denovo_tx = df_ass %>%
#                      pull(assemble_cpm), 
#                    txid = paste(df_ass$tx, rep(1:6, each = length(txid))), 
#                    denovo_cts = df_ass$counts.max,
#                    truth = log2(tx_counts_bambu$counts + 1)[txid, ] %>%
#                      as.matrix() %>%
#                      c(), 
#                    truth_counts = tx_counts_bambu$counts[txid, ] %>%
#                      as.matrix() %>%
#                      c())
#   
#   res <- list(tx = tx)
#   
# }) 
# 
# plots2 <- lapply(names(mat2), function(colname) {
#   df <- mat2[[colname]]$tx
#   ggplot(df, aes(x = denovo_tx, y = truth)) +
#     geom_hex(bins = 50) +
#     scale_fill_viridis_c(trans = "log10") +
#     labs(
#       x = colname,
#       y = 'true_lcpm',
#       title = paste("Correlation =", round(cor(df$denovo_tx, df$truth), 3))
#     ) +
#     theme_bw()
# })
# 
# library(patchwork)
# wrap_plots(plots2[6:7], ncol = 2)
# 
# mat2$rattle.sec$tx %>% 
#   inner_join(mat2$rnabloom2.sec$tx, by = 'txid', suffix = c('.rattle','.rnabloom2')) %>% 
#   select(denovo_tx.rattle, denovo_tx.rnabloom2) %>% 
#   cor()
# 
# # rnabloom2 and rattle shared assembled transcript only
# ggplot(mat2$rattle.sec$tx %>% 
#          inner_join(mat2$rnabloom2.sec$tx, by = 'txid', suffix = c('.rattle','.rnabloom2')) %>% 
#          select(denovo_tx.rattle, denovo_tx.rnabloom2) , 
#        aes(x = denovo_tx.rnabloom2, y = denovo_tx.rattle)) +
#   geom_hex(bins = 50) +
#   scale_fill_viridis_c(trans = "log10") +
#   labs(
#     x = 'rnabloom2',
#     y = 'rattle',
#     title = paste("Correlation =", 0.687)
#   ) +
#   theme_bw()

