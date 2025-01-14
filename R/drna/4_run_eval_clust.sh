#!/bin/bash
module load R/4.2.3

Rscript --vanilla eval_cluster.R ../rattle/A549_MCF7_directRNA_sqanti3 ../rattle/A549_MCF7_directRNA/fx2tab.txt rattle_rattle
Rscript --vanilla eval_cluster.R ../rattle/A549_MCF7_directRNA_sqanti3 ../rattle/A549_MCF7_directRNA_corset/A549_MCF7_directRNA-clusters_mod.txt rattle_corset

Rscript --vanilla eval_cluster.R ../rnabloom2/A549_MCF7_directRNA_sqanti3 ../rnabloom2/A549_MCF7_directRNA_corset/A549_MCF7_directRNA-clusters_mod.txt rnabloom2_corset

# Rscript --vanilla eval_cluster.R ../trinity/A549_MCF7_Illimina_sqanti3 ../trinity/A549_MCF7_Illimina/A549_MCF7_Illimina_merged_trinity.Trinity.fasta.gene_trans_map trinity_trinity
# Rscript --vanilla eval_cluster.R ../trinity/A549_MCF7_Illimina_sqanti3 ../trinity/A549_MCF7_Illimina_corset/A549_MCF7_Illimina-clusters_mod.txt trinity_corset

Rscript --vanilla eval_cluster.R ../trinitystranded/A549_MCF7_Illimina_sqanti3 ../trinitystranded/A549_MCF7_Illimina/A549_MCF7_Illimina_merged_trinity.Trinity.fasta.gene_trans_map trinity_trinity
Rscript --vanilla eval_cluster.R ../trinitystranded/A549_MCF7_Illimina_sqanti3 ../trinitystranded/A549_MCF7_Illimina_corset/A549_MCF7_Illimina-clusters_mod.txt trinity_corset

Rscript --vanilla eval_cluster.R ../isonform/A549_MCF7_directRNA_sqanti3 ../isonform/A549_MCF7_directRNA_corset/A549_MCF7_directRNA-clusters_mod.txt isonform_corset
Rscript --vanilla eval_cluster.R ../isonform/A549_MCF7_directRNA_sqanti3 ../isonform/A549_MCF7_directRNA/fx2tab.txt isonform_isonclust

Rscript --vanilla eval_cluster.R ../bambu/A549_MCF7_directRNA_sqanti3 ../bambu/A549_MCF7_directRNA_dge/filter/tx2gene.txt bambu_bambu
