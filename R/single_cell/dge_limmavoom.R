
#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(GGally))

args = commandArgs(trailingOnly=TRUE)

# outputdir <- "bambu_bambu_oarfish/"
# group_mon <- "0,1,6,7,9,10,11"
# group_lym <- "2,3,4,5,8"

outputdir <- args[1]
group_mon <- args[2]
group_lym <- args[3]

# print(outputdir)
# print(group_mon)
# print(group_lym)
group_mon <- as.numeric(strsplit(group_mon, ",")[[1]])
group_lym <- as.numeric(strsplit(group_lym, ",")[[1]])

dge_gene <- readRDS(paste0(outputdir, "/dgelist_gene.rds"))
counts <- dge_gene$counts
samples <- dge_gene$samples

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

x <- DGEList(counts = aggregated_counts,
             group = rep(c('Lymphoid','Monocytic'), 3),
             genes = dge_gene$genes)

saveRDS(x, paste0(outputdir, '/x.rds'))

keep <- filterByExpr(x, group=x$samples$group,
                     min.count = 50, min.total.count = 100)
sum(keep)

x <- x[keep,,keep.lib.sizes=FALSE]
x <- calcNormFactors(x)

# cpm <- cpm(x, log = T)
# cor(cpm)

# ggpairs(cpm %>% data.frame())
# ggsave(paste0(outputdir, '/cor_gene.pdf'), width = 5, height = 5)

pdf(paste0(outputdir, '/mds_gene.pdf'), width = 5, height = 5)
plotMDS(x,col=c(1:2)[x$samples$group])
dev.off()

design <- model.matrix(~0 + x$samples$group)
colnames(design) <- c('Lymphoid','Monocytic')
x <- estimateDisp(x, design, robust=TRUE)

pdf(paste0(outputdir, '/dge_bcvplot.pdf'), width = 5, height = 5)
plotBCV(x)
dev.off()

contr2 <- makeContrasts(Lymphoid - Monocytic, levels = design)

pdf(paste0(outputdir, '/dge_mvplot.pdf'), width = 5, height = 5)
v <- voom(x, design, plot=TRUE)
dev.off()

fit <- lmFit(v, design)
fit <- contrasts.fit(fit, contrasts = contr2)
fit.de <- eBayes(fit, robust=TRUE)
# dge <- topTable(fit.de, n = Inf #, p.value = 0.05
# )
fit.de <- treat(fit.de, fc = 1.5)
dge <- topTreat(fit.de, n = Inf)

is.de <- decideTests(fit.de,p.value=0.05)
summary(is.de)

pdf(paste0(outputdir, '/dge_mdplot.pdf'), width = 5, height = 5)
plotMD(fit.de, status = is.de)
dev.off()

saveRDS(dge, paste0(outputdir, '/dge_gene.rds'))

