#!/bin/bash
module load R/4.2.3

## 2m
Rscript --vanilla eval_cluster.R ../isonform/merged_2m_sqanti3/ ../isonform/merged_2m_corset/merged_2m-clusters_mod.txt isonform_corset_2m
Rscript --vanilla eval_cluster.R ../isonform/merged_2m_sqanti3/ ../isonform/merged_2m/fx2tab.txt isonform_isonclust_2m

Rscript --vanilla eval_cluster.R ../rattle/merged_2m_sqanti3/ ../rattle/merged_2m_corset/merged_2m-clusters_mod.txt rattle_corset_2m
Rscript --vanilla eval_cluster.R ../rattle/merged_2m_sqanti3/ ../rattle/merged_2m/fx2tab.txt  rattle_rattle_2m

Rscript --vanilla eval_cluster.R ../rnabloom2/merged_2m_sqanti3/ ../rnabloom2/merged_2m_corset/merged_2m-clusters_mod.txt rnabloom2_corset_2m

Rscript --vanilla eval_cluster.R ../bambu/merged_2m_sqanti3/ ../bambu/merged_2m_dge/filter/tx2gene.txt bambu_bambu_2m
# 
# ## 5m
Rscript --vanilla eval_cluster.R ../isonform/merged_5m_sqanti3/ ../isonform/merged_5m_corset/merged_5m-clusters_mod.txt isonform_corset_5m
Rscript --vanilla eval_cluster.R ../isonform/merged_5m_sqanti3/ ../isonform/merged_5m/fx2tab.txt isonform_isonclust_5m

Rscript --vanilla eval_cluster.R ../rattle/merged_5m_sqanti3/ ../rattle/merged_5m_corset/merged_5m-clusters_mod.txt rattle_corset_5m
Rscript --vanilla eval_cluster.R ../rattle/merged_5m_sqanti3/ ../rattle/merged_5m/fx2tab.txt  rattle_rattle_5m

Rscript --vanilla eval_cluster.R ../rnabloom2/merged_5m_sqanti3/ ../rnabloom2/merged_5m_corset/merged_5m-clusters_mod.txt rnabloom2_corset_5m

Rscript --vanilla eval_cluster.R ../bambu/merged_5m_sqanti3/ ../bambu/merged_5m_dge/filter/tx2gene.txt bambu_bambu_5m
# 
# ## 10m
Rscript --vanilla eval_cluster.R ../rattle/merged_10m_sqanti3/ ../rattle/merged_10m/fx2tab.txt  rattle_rattle_10m
Rscript --vanilla eval_cluster.R ../rattle/merged_10m_sqanti3/ ../rattle/merged_10m_corset/merged_10m-clusters_mod.txt rattle_corset_10m

Rscript --vanilla eval_cluster.R ../rnabloom2/merged_10m_sqanti3/ ../rnabloom2/merged_10m_corset/merged_10m-clusters_mod.txt rnabloom2_corset_10m

Rscript --vanilla eval_cluster.R ../bambu/merged_10m_sqanti3/ ../bambu/merged_10m_dge/filter/tx2gene.txt bambu_bambu_10m

Rscript --vanilla eval_cluster.R ../trinitystranded/merged_ilu_sqanti3/ ../trinitystranded/merged_ilu/merged_ilu_trinity.Trinity.fasta.gene_trans_map trinity_trinity_10m
Rscript --vanilla eval_cluster.R ../trinitystranded/merged_ilu_sqanti3/ ../trinitystranded/merged_ilu_corset/merged_ilu-clusters_mod.txt trinity_corset_10m


