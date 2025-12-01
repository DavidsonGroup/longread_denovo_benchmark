#!/bin/bash

#SBATCH --job-name=prepare_counts
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=10:00:00
#SBATCH --array=1-6
#SBATCH --output=logs/prepare_counts_%A_%a.out
#SBATCH --error=logs/prepare_counts_%A_%a.err

# Create logs directory if it doesn't exist
mkdir -p logs

module load R/4.4.2

# Define the parameters for each array job
case $SLURM_ARRAY_TASK_ID in
    1)  
        PREFIX="rattle"
        METHOD="corset"
        ;;
    2)  
        PREFIX="rattle"
        METHOD="rattle"
        ;;
    3)  
        PREFIX="rnabloom2"
        METHOD="corset"
        ;;
    4)  
        PREFIX="isonform"
        METHOD="corset"
        ;;
    5)  
        PREFIX="isonform"
        METHOD="isonclust"
        ;;
    6)  
        PREFIX="bambu"
        METHOD="bambu"
        ;;
    7)  
        PREFIX="gencode"
        METHOD="gencode"
        ;;    
    *)
        echo "Invalid array task ID"
        exit 1
        ;;
esac

# Set the cluster file path based on method
if [ "$METHOD" = "rattle" ] || [ "$METHOD" = "isonclust" ]; then
    CLUSTER_FILE="../${PREFIX}/fx2tab.txt"
elif [ "$METHOD" = "bambu" ]; then
    CLUSTER_FILE="../${PREFIX}/PB_merged_dge/filter/tx2gene.txt"
elif [ "$METHOD" = "gencode" ]; then
    CLUSTER_FILE="../${PREFIX}/tx2gene.txt"
else
    CLUSTER_FILE="../${PREFIX}/PB_merged_corset/PB_merged-clusters_mod.txt"
fi

# Run the analysis pipeline
NAME="${PREFIX}_${METHOD}_oarfish"

# Prepare counts
Rscript --vanilla 0_preparecounts.R "../${PREFIX}/" "${CLUSTER_FILE}" "${PREFIX}_${METHOD}" "oarfish"

# Run DGE and DTE analysis
Rscript --vanilla dge.R "${NAME}"
Rscript --vanilla dte.R "${NAME}"