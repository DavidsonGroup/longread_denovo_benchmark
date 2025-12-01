#!/bin/bash
module load R/4.4.2

# Rscript --vanilla 0_preparecounts.R ../rattle/ ../rattle/PB_merged_corset/PB_merged-clusters_mod.txt rattle_corset oarfish
Rscript --vanilla 0_preparecounts.R ../rattle/ ../rattle/fx2tab.txt rattle_rattle oarfish

# Rscript --vanilla 0_preparecounts.R ../rnabloom2/ ../rnabloom2/PB_merged_corset/PB_merged-clusters_mod.txt rnabloom2_corset oarfish

# Rscript --vanilla 0_preparecounts.R ../isonform/ ../isonform/PB_merged_corset/PB_merged-clusters_mod.txt isonform_corset oarfish
# Rscript --vanilla 0_preparecounts.R ../isonform/ ../isonform/fx2tab.txt isonform_isonclust oarfish

# for gene level, get clustering and pseudobulk gene DGE object
# Rscript --vanilla dge.R rattle_corset_oarfish
Rscript --vanilla dge.R rattle_rattle_oarfish

# Rscript --vanilla dge.R rnabloom2_corset_oarfish
# 
# Rscript --vanilla dge.R isonform_corset_oarfish
# Rscript --vanilla dge.R isonform_isonclust_oarfish

# for transcript level, get clustering and pseudobulk transcript DGE object
# Rscript --vanilla dte.R rattle_corset_oarfish
Rscript --vanilla dte.R rattle_rattle_oarfish

# Rscript --vanilla dte.R rnabloom2_corset_oarfish
# 
# Rscript --vanilla dte.R isonform_corset_oarfish
# Rscript --vanilla dte.R isonform_isonclust_oarfish

# add bambu
# Rscript --vanilla 0_preparecounts.R ../bambu/ ../bambu/PB_merged_dge/filter/tx2gene.txt bambu_bambu oarfish
# Rscript --vanilla dge.R bambu_bambu_oarfish
# Rscript --vanilla dte.R bambu_bambu_oarfish