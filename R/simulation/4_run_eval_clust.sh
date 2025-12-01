#!/bin/bash
module load R/4.2.3

Rscript --vanilla eval_cluster.R ../simulation_1m/rattle/ONT_merged_sqanti3 ../simulation_1m/rattle/fx2tab.txt rattle_rattle
Rscript --vanilla eval_cluster.R ../simulation_1m/rattle/ONT_merged_sqanti3 ../simulation_1m/rattle/ONT_merged_corset/ONT_merged-clusters_mod.txt rattle_corset

Rscript --vanilla eval_cluster.R ../simulation_1m/rnabloom2/ONT_merged_sqanti3 ../simulation_1m/rnabloom2/ONT_merged_corset/ONT_merged-clusters_mod.txt rnabloom2_corset

Rscript --vanilla eval_cluster.R ../simulation_1m/trinity/Ilu_merged_sqanti3 ../simulation_1m/trinity/Ilu_merged/Ilu_merged_trinity.Trinity.fasta.gene_trans_map trinity_trinity
Rscript --vanilla eval_cluster.R ../simulation_1m/trinity/Ilu_merged_sqanti3 ../simulation_1m/trinity/Ilu_merged_corset/Ilu_merged-clusters_mod.txt trinity_corset

Rscript --vanilla eval_cluster.R ../simulation_1m/isonform/ONT_merged_sqanti3 ../simulation_1m/isonform/ONT_merged_corset/ONT_merged-clusters_mod.txt isonform_corset
Rscript --vanilla eval_cluster.R ../simulation_1m/isonform/ONT_merged_sqanti3 ../simulation_1m/isonform/ONT_merged/fx2tab.txt isonform_isonclust

Rscript --vanilla eval_cluster.R ../simulation_1m/bambu/ONT_merged_sqanti3 ../simulation_1m/bambu/ONT_merged_dge/filter/tx2gene.txt bambu_bambu

# hybrid
Rscript --vanilla eval_cluster.R ../simulation_1m/rnaspades/hybrid_merged_sqanti3/ ../simulation_1m/rnaspades/hybrid_merged_corset/hybrid_merged-clusters_mod.txt rnaspades_corset
Rscript --vanilla eval_cluster.R ../simulation_1m/rnaspades/hybrid_merged_sqanti3/ ../simulation_1m/rnaspades/tx2gene.txt rnaspades_rnaspades

Rscript --vanilla eval_cluster.R ../simulation_1m/rnabloom2hybrid/hybrid_merged_sqanti3/ ../simulation_1m/rnabloom2hybrid/hybrid_merged_corset/hybrid_merged-clusters_mod.txt rnabloom2hybrid_corset

# bambu denovo 
Rscript --vanilla eval_cluster.R ../simulation_1m/bambudenovo/ONT_merged_sqanti3 ../simulation_1m/bambudenovo/ONT_merged_dge/filter/tx2gene.txt bambudenovo_bambu
Rscript --vanilla eval_cluster.R ../simulation_1m/bambudenovo/ONT_merged_sqanti3 ../simulation_1m/bambudenovo/ONT_merged_corset/ONT_merged-clusters_mod.txt bambudenovo_corset
