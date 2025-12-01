library(tidyverse)

set.seed(100)

args <- c('sqanti-sim_index_ctrl.tsv', 'sqanti-sim_index_de.tsv')

n.feat <- 40509
n.libs <- c(3,3)
lib.sizes <- rep(1000000, sum(n.libs))
bcv.true <- 0.2
df.bcv <- 40

n.groups <- length(n.libs)
n.samples <- sum(n.libs)

baselineAbundance <- lapply(args, function(x) {
  baselineAbundance <- read.table(x, header = T)$requested_counts
  baselineAbundance <- baselineAbundance/sum(baselineAbundance)
})

# baselineAbundance <- read.table(args[1], header = T)$requested_counts
# baselineAbundance <- baselineAbundance/sum(baselineAbundance)
# 
# baselineAbundancePerGroup <- matrix(baselineAbundance,
#                                     ncol = n.groups,nrow = n.feat)
  
baselineAbundancePerGroup <- cbind(baselineAbundance[[1]], 
                                   baselineAbundance[[2]])
  
# Generating expected counts
mu0 <- lapply(seq_len(n.groups),function(x){
  size <- lib.sizes[seq(1 + n.libs[x] * (x - 1), n.libs[x] * x)]
  matrix(baselineAbundancePerGroup[,x],n.feat,1) %*% matrix(size,1,length(size))
})
mu0 <- do.call(cbind,mu0)

# Generating random noise around dispersion trend. I am generating 1 RV per
# transcript per group. In the voom simulation, we had 1 RV per gene per sample.
# Here, I am arguing that the expression level and dispersion should be exactly
# the same among libraries of the same group. Since we already generated DE
# status in the steps above, it makes sense to have 1 random shift around the
# dispersion trend per group and a fixed resulting dispersion per group.
chisq <- lapply(seq_len(n.groups),function(x){
  rv <- df.bcv / rchisq(n.feat, df = df.bcv)
  matrix(rv,ncol = n.libs[x],nrow = n.feat)
})
chisq <- do.call(cbind,chisq)

# Biological variation and Dispersion trend
bcv0 <- bcv.true + 1/sqrt(mu0)
disp <- bcv0 ^ 2 * chisq

# Biological variation
shape <- 1/disp
scale <- mu0/shape

expr <- rgamma(n.feat * n.samples, shape = shape, scale = scale)

expr <- matrix(expr,nrow =  n.feat, ncol =  n.samples)

# sanity check for BCV plot
library(edgeR)
y <- DGEList(counts = cbind(rpois(n.feat, expr[,1]),
                            rpois(n.feat, expr[,2]),
                            rpois(n.feat, expr[,3]),
                            rpois(n.feat, expr[,4]),
                            rpois(n.feat, expr[,5]),
                            rpois(n.feat, expr[,6])),
             group = gl(2,3))

y <- normLibSizes(y)
y <- estimateDisp(y, robust = T)

plotBCV(y)

# write output for sqanti-sim
expr.l <- list(expr[,1:3], expr[,4:6])

for (i in 1:2) {
  for (j in 1:3) {
    
    table <- read.table(args[i], header = T)
    df1 <- table %>% mutate(requested_counts = rpois(n.feat, expr.l[[i]][,j]), 
                            requested_tpm = requested_counts/sum(requested_counts)*1e6)
    filename <- tools::file_path_sans_ext(args[i])
    write.table(df1, paste0(filename, j, '.tsv'),
                quote = F, row.names = F, sep = '\t')
  }
}
