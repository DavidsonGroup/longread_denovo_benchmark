#!/bin/bash
module load R/4.2.3

# for rattle: 3 counting method, 2 clustering method, 6 combo in total

# Rscript --vanilla dtu_limmavoom.R ../rattle/ ../rattle/fx2tab.txt rattle_rattle count
# Rscript --vanilla dge_limmavoom.R ../rattle/ ../rattle/fx2tab.txt rattle_rattle count
# 
# Rscript --vanilla dtu_limmavoom.R ../rattle/ ../rattle/ONT_merged_corset/ONT_merged-clusters_mod.txt rattle_corset count
# Rscript --vanilla dge_limmavoom.R ../rattle/ ../rattle/ONT_merged_corset/ONT_merged-clusters_mod.txt rattle_corset count

Rscript --vanilla dtu_limmavoom.R ../rattle/ ../rattle/fx2tab.txt rattle_rattle onts
Rscript --vanilla dge_limmavoom.R ../rattle/ ../rattle/fx2tab.txt rattle_rattle onts

Rscript --vanilla dtu_limmavoom.R ../rattle/ ../rattle/pea_corset/pea-clusters_mod.txt rattle_corset onts
Rscript --vanilla dge_limmavoom.R ../rattle/ ../rattle/pea_corset/pea-clusters_mod.txt rattle_corset onts

Rscript --vanilla dtu_limmavoom.R ../rattle/ ../rattle/fx2tab.txt rattle_rattle ontp
Rscript --vanilla dge_limmavoom.R ../rattle/ ../rattle/fx2tab.txt rattle_rattle ontp

Rscript --vanilla dtu_limmavoom.R ../rattle/ ../rattle/pea_corset/pea-clusters_mod.txt rattle_corset ontp
Rscript --vanilla dge_limmavoom.R ../rattle/ ../rattle/pea_corset/pea-clusters_mod.txt rattle_corset ontp

# for rnabloom2: 2 counting methods, 1 clustering method, 2 combo intotal

Rscript --vanilla dtu_limmavoom.R ../rnabloom2/ ../rnabloom2/pea_corset/pea-clusters_mod.txt rnabloom2_corset onts
Rscript --vanilla dge_limmavoom.R ../rnabloom2/ ../rnabloom2/pea_corset/pea-clusters_mod.txt rnabloom2_corset onts

Rscript --vanilla dtu_limmavoom.R ../rnabloom2/ ../rnabloom2/pea_corset/pea-clusters_mod.txt rnabloom2_corset ontp
Rscript --vanilla dge_limmavoom.R ../rnabloom2/ ../rnabloom2/pea_corset/pea-clusters_mod.txt rnabloom2_corset ontp

# for trinity: 2 counting methods, 2 clustering method, 4 combo intotal

Rscript --vanilla dtu_limmavoom.R ../trinitystranded/ ../trinitystranded/pea/pea_trinity.Trinity.fasta.gene_trans_map trinity_trinity align
Rscript --vanilla dge_limmavoom.R ../trinitystranded/ ../trinitystranded/pea/pea_trinity.Trinity.fasta.gene_trans_map trinity_trinity align

Rscript --vanilla dtu_limmavoom.R ../trinitystranded/ ../trinitystranded/pea_corset/pea-clusters_mod.txt trinity_corset align
Rscript --vanilla dge_limmavoom.R ../trinitystranded/ ../trinitystranded/pea_corset/pea-clusters_mod.txt trinity_corset align

Rscript --vanilla dtu_limmavoom.R ../trinitystranded/ ../trinitystranded/pea/pea_trinity.Trinity.fasta.gene_trans_map trinity_trinity map
Rscript --vanilla dge_limmavoom.R ../trinitystranded/ ../trinitystranded/pea/pea_trinity.Trinity.fasta.gene_trans_map trinity_trinity map

Rscript --vanilla dtu_limmavoom.R ../trinitystranded/ ../trinitystranded/pea_corset/pea-clusters_mod.txt trinity_corset map
Rscript --vanilla dge_limmavoom.R ../trinitystranded/ ../trinitystranded/pea_corset/pea-clusters_mod.txt trinity_corset map


