#!/bin/bash
# Usage: sbatch slurm-serial-job-script
# Prepared By: Alex Yan
#              yan.a@wehi.edu.au

# NOTE: To activate a SLURM option, remove the whitespace between the '#' and 'SBATCH'

# To give your job a name, replace "MyJob" with an appropriate name
#SBATCH --job-name=rnabloom2

# To set a project account for credit charging,
# SBATCH --account=ls25

# Request CPU resource for a serial job
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48

# Memory usage (MB)
#SBATCH --mem-per-cpu=8000

# Set your minimum acceptable walltime, format: day-hours:minutes:seconds
#SBATCH --time=48:00:00

# To receive an email when job completes or fails
#SBATCH --mail-user=yan.a@wehi.edu.au
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=BEGIN

# SBATCH --array=1-4

# Set the file for output (stdout)
# SBATCH --output=stdout

# Set the file for error log (stderr)
# SBATCH --error=stderr

# Use reserved node to run job when a node reservation is made for you already
# SBATCH --reservation=reservation_name

# SBATCH --partition=genomics
# SBATCH --qos=genomics

# SBATCH --dependency=afterok:10567760
# SBATCH --reservation=highmem

# Command to run a serial job

module purge
module load anaconda3/2019.03
conda activate /stornext/Bioinf/data/lab_davidson/yan.a/conda_env/rnabloom

# raw_reads='/vast/scratch/users/yan.a/vast_scratch/pea_fastq_merged/pea_fulllength_rescued_merged.fastq.gz'
raw_reads='/home/users/allstaff/yan.a/lab_davidson/yan.a/pea_fastq_merged/pea_fulllength_rescued_merged.fastq' # because rattle used uncompressed 
num_cores=48
base=$(basename $raw_reads _fulllength_rescued_merged.fastq)

export JAVA_TOOL_OPTIONS="-Xmx300g"
rnabloom -long ${raw_reads} -stranded -t 48 -outdir ${base}

