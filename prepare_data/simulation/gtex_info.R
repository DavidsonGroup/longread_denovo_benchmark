# get baseline expression from GTex

# setwd('~/davidson_longread/yan.a/gtex_quant/')

library(tidyverse)
library(edgeR)
library(tximport)

files.quant <- list.files('.', pattern = 'quant.sf', recursive = T, full.names = T)
names <- basename(dirname(files.quant)) %>% str_remove('\\.FAK.*onts')

df <- tximport::tximport(files.quant,
                         type = 'salmon', txIn = T, txOut = T,
                         ignoreAfterBar = T)

colnames(df$counts) <- names
saveRDS(df,'df.rds')

cts <- df$counts
rownames(cts) <- rownames(df$counts)
is.expr <- rowSums(cts > 5) >= 15
sum(is.expr)
# [1] 40581

cts <- cts[is.expr, ]
prop <- goodTuringProportions(round(cts))
baselineAbundance <- rowMeans(prop)

# saveRDS(baselineAbundance, 'baselineAbundance.rds')
# write.table(names(baselineAbundance), 'txid.txt', row.names = F, quote = F, col.names = F)

## use catchSalmon to get transcript length 
df2 <- catchSalmon(dirname(files.quant))
colnames(df2$counts) <- names
rownames(df2$counts) <- rownames(df$counts)

identical(df$counts, df2$counts)

cts <- df2$counts[df2$annotation$Length > 200, ] ## remove small transcript < 200bp 
# rownames(cts) <- rownames(df$counts)
is.expr <- rowSums(cts > 5) >= 15 ## remove lowly expressed transcript
sum(is.expr)
# [1] 40509

cts <- cts[is.expr, ]
prop <- goodTuringProportions(round(cts))
baselineAbundance <- rowMeans(prop)

saveRDS(baselineAbundance, 'baselineAbundance.rds')
write.table(names(baselineAbundance), 'txid.txt', row.names = F, quote = F, col.names = F)



