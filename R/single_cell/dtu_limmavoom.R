
#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(GGally))

args = commandArgs(trailingOnly=TRUE)

# outputdir <- "rattle_corset_oarfish"
# group_mon <- c(0,3,5,8,9,10,11,12)
# group_lym <- c(1,2,4,6,7)

outputdir <- args[1]
group_mon <- args[2]
group_lym <- args[3]

# print(outputdir)
# print(group_mon)
# print(group_lym)
group_mon <- as.numeric(strsplit(group_mon, ",")[[1]])
group_lym <- as.numeric(strsplit(group_lym, ",")[[1]])

dge_tx <- readRDS(paste0(outputdir, "/dgelist_tx.rds"))
dge_gene <- readRDS(paste0(outputdir, "/dgelist_gene.rds")) # to get gene level clustering

# dge_tx <- readRDS("rattle_corset_oarfish/dgelist_tx.rds")
# dge_gene <- readRDS("rattle_corset_oarfish/dgelist_gene.rds") # to get gene level clustering

counts <- dge_tx$counts
samples <- dge_tx$samples

# Assuming your count matrix is called 'counts' and sample info is 'samples'

# Create a grouping column
samples$cluster_group <- ifelse(samples$cluster %in% group_mon, "Monocytic",
                                ifelse(samples$cluster %in% group_lym, "Lymphoid", NA))

sum(is.na(samples$cluster_group))

# Create combined sample-cluster group identifier
samples$new_group <- paste(samples$sample, samples$cluster_group, sep="_")

# Aggregate counts by summing columns belonging to the same group
aggregated_counts <- sapply(unique(samples$new_group), function(group) {
  cols_to_sum <- which(samples$new_group == group)
  if(length(cols_to_sum) == 1) {
    return(counts[, cols_to_sum])
  } else {
    return(rowSums(counts[, cols_to_sum]))
  }
})

dim(aggregated_counts)

y <- DGEList(counts = aggregated_counts, 
             genes = dge_tx$genes,
             group = rep(c('Lymphoid','Monocytic'), 3))

saveRDS(y, paste0(outputdir, '/y.rds'))

keep <- filterByExpr(y, group=y$samples$group, 
                    min.count = 20, min.total.count = 50)
sum(keep)
y <- y[keep,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)

# cpm <- cpm(y, log = T)
# # cor(cpm)
# 
# ggpairs(cpm %>% data.frame())
# ggsave(paste0(outputdir, '/cor_tx.pdf'), width = 5, height = 5)

pdf(paste0(outputdir, '/mds_tx.pdf'), width = 5, height = 5)
plotMDS(y,col=c(1:2)[y$samples$group])
dev.off()

design <- model.matrix(~0 + y$samples$group)
colnames(design) <- c('Lymphoid','Monocytic')

contr <- makeContrasts(Lymphoid - Monocytic, levels = design)

v <- voom(y, design, plot=TRUE)

fit <- lmFit(v, design)
fit <- contrasts.fit(fit, contrasts = contr)
fit.de <- eBayes(fit, robust=TRUE)

is.de <- decideTests(fit.de, p.value=0.05)
summary(is.de)

plotMD(fit.de, status = is.de)
ggsave(paste0(outputdir, '/dte_mdplot.pdf'))

# dte <- topTable(fit.de, n = Inf #, p.value = 0.05
# )
fit.de <- treat(fit.de, fc = 1.5)
dte <- topTreat(fit.de, n = Inf)

sp <- diffSplice(fit, 
                 geneid="geneid",exonid="isoform")

dtu_tx <- topSplice(sp, test="t", number = Inf #, FDR = 0.05
)

dtu_gene <- topSplice(sp, test="simes", number = Inf#, FDR = 0.05
)

saveRDS(dte, paste0(outputdir, '/dte_tx.rds'))
saveRDS(dtu_gene, paste0(outputdir, '/dtu_gene.rds'))
saveRDS(dtu_tx, paste0(outputdir, '/dtu_tx.rds'))

