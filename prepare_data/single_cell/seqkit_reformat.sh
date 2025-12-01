#!/bin/bash
# Usage: sbatch slurm-serial-job-script
# Prepared By: Alex Yan
#              yan.a@wehi.edu.au

# NOTE: To activate a SLURM option, remove the whitespace between the '#' and 'SBATCH'

# To give your job a name, replace "MyJob" with an appropriate name
#SBATCH --job-name=seqkit

# To set a project account for credit charging,
# SBATCH --account=ls25

# Request CPU resource for a serial job
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8

# Memory usage (MB)
#SBATCH --mem-per-cpu=8000

# Set your minimum acceptable walltime, format: day-hours:minutes:seconds
#SBATCH --time=48:00:00

# To receive an email when job completes or fails
#SBATCH --mail-user=yan.a@wehi.edu.au
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=BEGIN

#SBATCH --array=1-3

# Set the file for output (stdout)
#SBATCH --output=reformat/%x_%A_%a.out

# Set the file for error log (stderr)
#SBATCH --error=reformat/%x_%A_%a.err

# Use reserved node to run job when a node reservation is made for you already
# SBATCH --reservation=reservation_name

# SBATCH --partition=genomics
# SBATCH --qos=genomics

# SBATCH --dependency=afterok:10567760
# SBATCH --reservation=highmem

# Command to run a serial job

module purge
module load micromamba

eval "$(micromamba shell hook --shell bash)"

micromamba activate /stornext/Bioinf/data/lab_davidson/yan.a/conda_env/isoncorrect

mkdir -p reformat

J=$SLURM_ARRAY_TASK_ID

reads=$(sed "${J}q;d" full_fasta.txt)
# name=$(basename "$(dirname "$raw_reads")")
name="$(basename "$reads")"
name="${name%%.*}"

seqkit replace -j 8 \
  -p '^(\S+).*XM=([^;]+);CB=([^;]+).*' \
  -r $'$1\tUB:Z:$2\tCB:Z:$3' "$reads" > "reformat/${name}.fa.gz"

seqkit stats -j 8 "$reads" "reformat/${name}.fa.gz" > "reformat/${name}_stats.txt"