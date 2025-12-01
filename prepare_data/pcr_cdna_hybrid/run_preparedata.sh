#!/bin/bash
# Usage: sbatch slurm-serial-job-script
# Prepared By: Alex Yan
#              yan.a@wehi.edu.au

# NOTE: To activate a SLURM option, remove the whitespace between the '#' and 'SBATCH'

# To give your job a name, replace "MyJob" with an appropriate name
#SBATCH --job-name=hybrid

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

# SBATCH --partition=bigmem
# SBATCH --qos=genomics

# SBATCH --dependency=afterok:10567760
# SBATCH --reservation=highmem

# Command to run a serial job


module purge

## individual reads for DE
for i in /vast/projects/lab_davidson/yan.a/Dong_2023/subsampleshort/trim/*1P.fastq.gz; 
do
    /home/users/allstaff/yan.a/yan.a/software/seqtk/seqtk sample -s100 ${i} 0.5 > $(basename $i _1P.fastq.gz)_hybrid_sub_1.fastq.gz
done

for i in /vast/projects/lab_davidson/yan.a/Dong_2023/subsampleshort/trim/*2P.fastq.gz; 
do
    /home/users/allstaff/yan.a/yan.a/software/seqtk/seqtk sample -s100 ${i} 0.5 > $(basename $i _2P.fastq.gz)_hybrid_sub_2.fastq.gz
done

for i in /vast/projects/lab_davidson/yan.a/Dong_2023/subsample/barcode*_10m.fq.gz; 
do
    base=$(basename $i _10m.fq.gz)
    /home/users/allstaff/yan.a/yan.a/software/seqtk/seqtk sample -s100 ${i} 0.5 > ${base}_hybrid_LR.fq.gz
done

module purge
module load micromamba
eval "$(micromamba shell hook --shell bash)"
micromamba activate isoncorrect

cat *_hybrid_sub_1.fastq.gz | seqkit rename -j 12 -n -1 > hybrid_sub_1.fastq.gz
cat *_hybrid_sub_2.fastq.gz | seqkit rename -j 12 -n -1 > hybrid_sub_2.fastq.gz
cat *_hybrid_LR.fq.gz | seqkit rename -j 12 -n -1 > hybrid_LR.fq.gz 

seqkit stats hybrid_sub_1.fastq.gz hybrid_sub_2.fastq.gz hybrid_LR.fq.gz