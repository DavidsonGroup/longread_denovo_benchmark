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
#SBATCH --mem-per-cpu=10000

# Set your minimum acceptable walltime, format: day-hours:minutes:seconds
#SBATCH --time=48:00:00

# To receive an email when job completes or fails
#SBATCH --mail-user=yan.a@wehi.edu.au
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=BEGIN

# SBATCH --array=3

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
module load micromamba

eval "$(micromamba shell hook --shell bash)"
micromamba activate /stornext/Bioinf/data/lab_davidson/yan.a/conda_env/rnabloom

#J=1
read1='/vast/projects/lab_davidson/yan.a/Dong_2023/hybrid_reads/hybrid_sub_1.fastq'
read2='/vast/projects/lab_davidson/yan.a/Dong_2023/hybrid_reads/hybrid_sub_2.fastq'
long_reads="/vast/projects/lab_davidson/yan.a/Dong_2023/hybrid_reads/hybrid_LR.fq"

num_cores=48
base='hybrid_merged'

export JAVA_TOOL_OPTIONS="-Xmx480g"

# remove stranded setting
rnabloom -long ${long_reads} -left ${read1} -right ${read2} -t ${num_cores} -outdir ${base} -stranded

