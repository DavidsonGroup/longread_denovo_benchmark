#!/bin/bash
module load R/4.2.3

# for isonform: 2 counting method, 2 clustering method, 4 combo in total
for i in {2,5}
do
Rscript --vanilla dtu_limmavoom.R ../isonform/ ../isonform/merged_${i}m/fx2tab.txt  isonform_isonclust_${i}m onts
Rscript --vanilla dge_limmavoom.R ../isonform/ ../isonform/merged_${i}m/fx2tab.txt  isonform_isonclust_${i}m onts

Rscript --vanilla dtu_limmavoom.R ../isonform/ ../isonform/merged_${i}m_corset/merged_${i}m-clusters_mod.txt isonform_corset_${i}m onts
Rscript --vanilla dge_limmavoom.R ../isonform/ ../isonform/merged_${i}m_corset/merged_${i}m-clusters_mod.txt isonform_corset_${i}m onts

Rscript --vanilla dtu_limmavoom.R ../isonform/ ../isonform/merged_${i}m/fx2tab.txt  isonform_isonclust_${i}m ontp
Rscript --vanilla dge_limmavoom.R ../isonform/ ../isonform/merged_${i}m/fx2tab.txt  isonform_isonclust_${i}m ontp

Rscript --vanilla dtu_limmavoom.R ../isonform/ ../isonform/merged_${i}m_corset/merged_${i}m-clusters_mod.txt isonform_corset_${i}m ontp
Rscript --vanilla dge_limmavoom.R ../isonform/ ../isonform/merged_${i}m_corset/merged_${i}m-clusters_mod.txt isonform_corset_${i}m ontp
done

# for rattle: 3 counting method, 2 clustering method, 6 combo in total
for i in {2,5,10}
do
# Rscript --vanilla dtu_limmavoom.R ../simulation_1m/rattle/ ../simulation_1m/rattle/fx2tab.txt rattle_rattle count
# Rscript --vanilla dge_limmavoom.R ../simulation_1m/rattle/ ../simulation_1m/rattle/fx2tab.txt rattle_rattle count
# 
# Rscript --vanilla dtu_limmavoom.R ../simulation_1m/rattle/ ../simulation_1m/rattle/ONT_merged_corset/ONT_merged-clusters_mod.txt rattle_corset count
# Rscript --vanilla dge_limmavoom.R ../simulation_1m/rattle/ ../simulation_1m/rattle/ONT_merged_corset/ONT_merged-clusters_mod.txt rattle_corset count

Rscript --vanilla dtu_limmavoom.R ../rattle/ ../rattle/merged_${i}m/fx2tab.txt  rattle_rattle_${i}m onts
Rscript --vanilla dge_limmavoom.R ../rattle/ ../rattle/merged_${i}m/fx2tab.txt  rattle_rattle_${i}m onts

Rscript --vanilla dtu_limmavoom.R ../rattle/ ../rattle/merged_${i}m_corset/merged_${i}m-clusters_mod.txt rattle_corset_${i}m onts
Rscript --vanilla dge_limmavoom.R ../rattle/ ../rattle/merged_${i}m_corset/merged_${i}m-clusters_mod.txt rattle_corset_${i}m onts

Rscript --vanilla dtu_limmavoom.R ../rattle/ ../rattle/merged_${i}m/fx2tab.txt  rattle_rattle_${i}m ontp
Rscript --vanilla dge_limmavoom.R ../rattle/ ../rattle/merged_${i}m/fx2tab.txt  rattle_rattle_${i}m ontp

Rscript --vanilla dtu_limmavoom.R ../rattle/ ../rattle/merged_${i}m_corset/merged_${i}m-clusters_mod.txt rattle_corset_${i}m ontp
Rscript --vanilla dge_limmavoom.R ../rattle/ ../rattle/merged_${i}m_corset/merged_${i}m-clusters_mod.txt rattle_corset_${i}m ontp
done

# for rnabloom2: 2 counting methods, 1 clustering method, 2 combo intotal
for i in {2,5,10}
do
Rscript --vanilla dtu_limmavoom.R ../rnabloom2/ ../rnabloom2/merged_${i}m_corset/merged_${i}m-clusters_mod.txt rnabloom2_corset_${i}m onts
Rscript --vanilla dge_limmavoom.R ../rnabloom2/ ../rnabloom2/merged_${i}m_corset/merged_${i}m-clusters_mod.txt rnabloom2_corset_${i}m onts

Rscript --vanilla dtu_limmavoom.R ../rnabloom2/ ../rnabloom2/merged_${i}m_corset/merged_${i}m-clusters_mod.txt rnabloom2_corset_${i}m ontp
Rscript --vanilla dge_limmavoom.R ../rnabloom2/ ../rnabloom2/merged_${i}m_corset/merged_${i}m-clusters_mod.txt rnabloom2_corset_${i}m ontp
done

# for trinity: 2 counting methods, 2 clustering method, 4 combo intotal

Rscript --vanilla dtu_limmavoom.R ../trinitystranded/ ../trinitystranded/merged_ilu/merged_ilu_trinity.Trinity.fasta.gene_trans_map trinity_trinity_10m align
Rscript --vanilla dge_limmavoom.R ../trinitystranded/ ../trinitystranded/merged_ilu/merged_ilu_trinity.Trinity.fasta.gene_trans_map trinity_trinity_10m align

Rscript --vanilla dtu_limmavoom.R ../trinitystranded/ ../trinitystranded/merged_ilu_corset/merged_ilu-clusters_mod.txt trinity_corset_10m align
Rscript --vanilla dge_limmavoom.R ../trinitystranded/ ../trinitystranded/merged_ilu_corset/merged_ilu-clusters_mod.txt trinity_corset_10m align

Rscript --vanilla dtu_limmavoom.R ../trinitystranded/ ../trinitystranded/merged_ilu/merged_ilu_trinity.Trinity.fasta.gene_trans_map trinity_trinity_10m map
Rscript --vanilla dge_limmavoom.R ../trinitystranded/ ../trinitystranded/merged_ilu/merged_ilu_trinity.Trinity.fasta.gene_trans_map trinity_trinity_10m map

Rscript --vanilla dtu_limmavoom.R ../trinitystranded/ ../trinitystranded/merged_ilu_corset/merged_ilu-clusters_mod.txt trinity_corset_10m map
Rscript --vanilla dge_limmavoom.R ../trinitystranded/ ../trinitystranded/merged_ilu_corset/merged_ilu-clusters_mod.txt trinity_corset_10m map

