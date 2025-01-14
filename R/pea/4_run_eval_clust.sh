#!/bin/bash
module load R/4.2.3

Rscript --vanilla eval_cluster.R ../rattle/pea_sqanti3 ../rattle/fx2tab.txt rattle_rattle
Rscript --vanilla eval_cluster.R ../rattle/pea_sqanti3 ../rattle/pea_corset/pea-clusters_mod.txt rattle_corset

Rscript --vanilla eval_cluster.R ../rnabloom2/pea_sqanti3 ../rnabloom2/pea_corset/pea-clusters_mod.txt rnabloom2_corset

Rscript --vanilla eval_cluster.R ../trinitystranded/pea_sqanti3 ../trinitystranded/pea/pea_trinity.Trinity.fasta.gene_trans_map trinity_trinity
Rscript --vanilla eval_cluster.R ../trinitystranded/pea_sqanti3 ../trinitystranded/pea_corset/pea-clusters_mod.txt trinity_corset

Rscript --vanilla eval_cluster.R ../bambu/pea_sqanti3 ../bambu/pea_dge/filter/tx2gene.txt bambu_bambu
