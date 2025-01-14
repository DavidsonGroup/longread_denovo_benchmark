#!/bin/bash
# Usage: sbatch slurm-serial-job-script
# Prepared By: Alex Yan
#              yan.a@wehi.edu.au

# NOTE: To activate a SLURM option, remove the whitespace between the '#' and 'SBATCH'

# To give your job a name, replace "MyJob" with an appropriate name
#SBATCH --job-name=trinity

# To set a project account for credit charging,
# SBATCH --account=ls25

# Request CPU resource for a serial job
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48

# Memory usage (MB)
#SBATCH --mem-per-cpu=10000

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
module load singularity

read1='/home/users/allstaff/yan.a/vast_scratch/drna/trinitystranded/A549_MCF7_Illimina_merged_R1.fastq.gz'
read2=$(echo $read1 | sed "s/R1/R2/" )

num_cores=48
base=$(basename $read1 _R1.fastq.gz)

binddir='/vast/scratch/users/yan.a/vast_scratch/'

mkdir $base
cd $base

# remove --monitoring
singularity exec --bind $binddir -e /stornext/Bioinf/data/lab_davidson/yan.a/software/trinityrnaseq.v2.15.1.simg Trinity \
          --seqType fq --SS_lib_type RF --max_memory 480G --CPU 48 --trimmomatic --full_cleanup \
          --left ${read1} --right ${read2} \
          --output ${base}_trinity
