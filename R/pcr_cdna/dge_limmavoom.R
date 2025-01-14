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

# args[1] <- '../rattle/'
# args[2] <- "../rattle/A549_MCF7_directRNA_corset/A549_MCF7_directRNA-clusters_mod.txt" 
# args[3] <- 'rattle_corset'

# args <- c('../rattle/', '../rattle/merged_2m/fx2tab.txt',  'rattle_rattle_2m', 'onts')
outputdir <- paste(args[3], args[4], sep = '_')

dir.create(outputdir)

dirs <- list.dirs(args[1], recursive = F, full.names = T)
data <- str_split(args[3],'_', simplify = T)
data <- ifelse(data[1] == 'trinity' & data[3] == '10m', 'ilu', data[3]) # for short read name
dirs <- dirs[grepl(data, dirs)]

files <- list.files(dirs[grepl('dge', dirs)], 
                    pattern = 'quant.sf', full.names = T, recursive = T)

files <- files[grepl(args[4], files)]

method <- str_split(args[3], pattern = '_', simplify = T)[2]

rt_gene <- get_gene_tx(args[2], 
                       list.files(dirs[grepl('sqanti3$', dirs)], 
                                  pattern = 'classification.txt', 
                                  full.names = T,recursive = T), 
                       method = method)
rownames(rt_gene) <- rt_gene$isoform

salmon_gene <- tximport(files, 
                        type = "salmon", 
                        tx2gene = rt_gene[,c('isoform', 'geneid')])

x <- DGEList(counts = salmon_gene$counts, 
             #genes = rt_gene,
             group = rep(c('H1975','HCC827'), each = 3))
colnames(x) <- c('ctrl1','ctrl2', 'ctrl3','de1','de2','de3')

saveRDS(x, paste0(outputdir, '/x.rds'))

keep <- filterByExpr(x, group=x$samples$group, min.count = 5)
x <- x[keep,,keep.lib.sizes=FALSE]
x <- calcNormFactors(x)

cpm <- cpm(x, log = T)
# cor(cpm)
  
ggpairs(cpm %>% data.frame())
ggsave(paste0(outputdir, '/cor_gene.pdf'), width = 5, height = 5)

pdf(paste0(outputdir, '/mds_gene.pdf'), width = 5, height = 5)
plotMDS(x,col=c(1:2)[x$samples$group])
dev.off()

design <- model.matrix(~0 + x$samples$group)
colnames(design) <- c('H1975','HCC827')
x <- estimateDisp(x, design, robust=TRUE)

pdf(paste0(outputdir, '/dge_bcvplot.pdf'), width = 5, height = 5)
plotBCV(x)
dev.off()

contr2 <- makeContrasts(H1975 - HCC827, levels = design)

pdf(paste0(outputdir, '/dge_mvplot.pdf'), width = 5, height = 5)
v <- voom(x, design, plot=TRUE)
dev.off()

fit <- lmFit(v, design)
fit <- contrasts.fit(fit, contrasts = contr2)
fit.de <- eBayes(fit, robust=TRUE)
dge <- topTable(fit.de, n = Inf #, p.value = 0.05
                )

is.de <- decideTests(fit.de,p.value=0.05)
summary(is.de)

pdf(paste0(outputdir, '/dge_mdplot.pdf'),width = 5, height = 5)
plotMD(fit.de, status = is.de)
dev.off()

saveRDS(dge, paste0(outputdir, '/dge_gene.rds'))




