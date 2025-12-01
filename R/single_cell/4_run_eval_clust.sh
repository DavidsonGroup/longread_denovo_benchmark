#!/bin/bash
module load R/4.2.3

Rscript --vanilla eval_cluster.R ../rattle/PB_merged_sqanti3 ../rattle/fx2tab.txt rattle_rattle
Rscript --vanilla eval_cluster.R ../rattle/PB_merged_sqanti3 ../rattle/PB_merged_corset/PB_merged-clusters_mod.txt rattle_corset

Rscript --vanilla eval_cluster.R ../rnabloom2/PB_merged_sqanti3 ../rnabloom2/PB_merged_corset/PB_merged-clusters_mod.txt rnabloom2_corset

Rscript --vanilla eval_cluster.R ../isonform/PB_merged_sqanti3 ../isonform/PB_merged_corset/PB_merged-clusters_mod.txt isonform_corset
Rscript --vanilla eval_cluster.R ../isonform/PB_merged_sqanti3 ../isonform/fx2tab.txt isonform_isonclust

Rscript --vanilla eval_cluster.R ../bambu/PB_merged_sqanti3 ../bambu/PB_merged_dge/filter/tx2gene.txt bambu_bambu
