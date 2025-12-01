#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(GGally))

source('~/lab_davidson/yan.a/software/scripts_denovo/R/get_gene_tx.R')

args = commandArgs(trailingOnly=TRUE)

# args[1] <- c('../simulation_1m/rattle/') # parent folder to get method and sqanti3 files
# args[2] <- '../simulation_1m/rattle/fx2tab.txt' # gene tx mapping
# args[3] <- 'rattle_rattle' # output prefix
# args[4] <- 'onts' # quant file to use

# args[1] <- c('../simulation_1m/rattle/') # parent folder to get method and sqanti3 files
# args[2] <- '../simulation_1m/rattle/ONT_merged_corset/ONT_merged-clusters_mod.txt' # gene tx mapping
# args[3] <- 'rattle_corset' # output prefix
# args[4] <- 'onts' # quant file to use

# args[1] <- '../simulation_1m/trinity/'
# args[2] <- "../simulation_1m/trinity/Ilu_merged/Ilu_merged_trinity.Trinity.fasta.gene_trans_map"
# args[3] <- 'trinity_trinity'
# args[1] <- '../simulation_1m/rnabloom2/'
# args[2] <- "../simulation_1m/rnabloom2/ONT_merged_corset/ONT_merged-clusters_mod.txt"
# args[3] <- 'rnabloom2_corset'

dirs <- list.dirs(args[1], recursive = F, full.names = T)

method <- str_split(args[3], pattern = '_', simplify = T)[2]

rt_gene <- get_gene_tx(args[2], 
                       list.files(dirs[grepl('sqanti3$', dirs)], pattern = 'classification.txt', full.names = T), 
                       method = method)
rownames(rt_gene) <- rt_gene$isoform

# if args[5] is not NA, means hybrid we need to add up count matrix
if (is.na(args[5])) {
  outputdir <- paste(args[3], args[4], sep = '_')
  
  dir.create(outputdir)
  
  files <- list.files(dirs[grepl('dge', dirs)], 
                      pattern = 'quant.sf', full.names = T, recursive = T)
  
  files <- files[grepl(args[4], files)]
  
  salmon_tx <- catchSalmon(dirname(files))
  colnames(salmon_tx$counts) <- basename(dirname(files))
  pdf(paste0(outputdir, '/hist_disp.pdf'))
  hist(salmon_tx$annotation$Overdispersion)
  dev.off()
  
  salmon_counts <- salmon_tx$counts/salmon_tx$annotation$Overdispersion
  
} else {
  
  ext <- "sum"
  outputdir <- paste(args[3], ext, sep = '_')
  
  dir.create(outputdir)
  
  files <- list.files(dirs[grepl('dge', dirs)], 
                      pattern = 'quant.sf', full.names = T, recursive = T)
  
  files1 <- files[grepl(args[4], files)]
  files2 <- files[grepl(args[5], files)]
  
  salmon_tx1 <- catchSalmon(dirname(files1))
  colnames(salmon_tx1$counts) <- basename(dirname(files1))
  pdf(paste0(outputdir, '/hist_disp1.pdf'))
  hist(salmon_tx1$annotation$Overdispersion)
  dev.off()
  
  salmon_tx2 <- catchSalmon(dirname(files2))
  colnames(salmon_tx2$counts) <- basename(dirname(files2))
  pdf(paste0(outputdir, '/hist_disp2.pdf'))
  hist(salmon_tx2$annotation$Overdispersion)
  dev.off()
  
  stopifnot(all.equal(rownames(salmon_tx1$counts), rownames(salmon_tx2$counts)))
  # salmon_tx <- salmon_tx1
  salmon_counts <- salmon_tx1$counts/salmon_tx1$annotation$Overdispersion + salmon_tx2$counts/salmon_tx2$annotation$Overdispersion
}

y <- DGEList(counts = salmon_counts, 
             genes = rt_gene[rownames(salmon_counts),],
             group = rep(c('ctrl','de'), each = 3))
colnames(y) <- c('ctrl1','ctrl2', 'ctrl3','de1','de2','de3')

saveRDS(y, paste0(outputdir, '/y.rds'))

keep <-filterByExpr(y, group=y$samples$group, min.count = 2, min.total.count = 5)
y <- y[keep,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)

cpm <- cpm(y, log = T)
# cor(cpm)

ggpairs(cpm %>% data.frame())
ggsave(paste0(outputdir, '/cor_tx.pdf'), width = 5, height = 5)

pdf(paste0(outputdir, '/mds_tx.pdf'), width = 5, height = 5)
plotMDS(y,col=c(1:2)[y$samples$group])
dev.off()

design <- model.matrix(~0 + y$samples$group)
colnames(design) <- c('ctrl','de')

contr <- makeContrasts(ctrl - de, levels = design)

v <- voom(y, design, plot=TRUE)

fit <- lmFit(v, design)
fit <- contrasts.fit(fit, contrasts = contr)
fit.de <- eBayes(fit, robust=TRUE)

is.de <- decideTests(fit.de, p.value=0.05)
summary(is.de)

plotMD(fit.de, status = is.de)
ggsave(paste0(outputdir, '/dte_mdplot.pdf'))

dte <- topTable(fit.de, n = Inf #, p.value = 0.05
                )

sp <- diffSplice(fit, 
                 geneid="geneid",exonid="isoform")

dtu_tx <- topSplice(sp, test="t", number = Inf #, FDR = 0.05
                    )

dtu_gene <- topSplice(sp, test="simes", number = Inf#, FDR = 0.05
                      )

saveRDS(dte, paste0(outputdir, '/dte_tx.rds'))
saveRDS(dtu_gene, paste0(outputdir, '/dtu_gene.rds'))
saveRDS(dtu_tx, paste0(outputdir, '/dtu_tx.rds'))

