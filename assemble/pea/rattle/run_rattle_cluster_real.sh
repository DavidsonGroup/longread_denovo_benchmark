#!/bin/bash
# Usage: sbatch slurm-serial-job-script
# Prepared By: Alex Yan
#              yan.a@wehi.edu.au

# NOTE: To activate a SLURM option, remove the whitespace between the '#' and 'SBATCH'

# To give your job a name, replace "MyJob" with an appropriate name
#SBATCH --job-name=rattle

# To set a project account for credit charging,
# SBATCH --account=ls25

# Request CPU resource for a serial job
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48

# Memory usage (MB)
#SBATCH --mem-per-cpu=10000

# Set your minimum acceptable walltime, format: day-hours:minutes:seconds
#SBATCH --time=14-00:00:00

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

#SBATCH --partition=long
# SBATCH --qos=genomics

# SBATCH --dependency=afterok:10567760
# SBATCH --reservation=highmem

# Command to run a serial job

module load gcc/9.1.0

raw_reads='/vast/scratch/users/yan.a/vast_scratch/pea_fastq_merged/pea_fulllength_rescued_merged.fastq.gz'
num_cores=48
base=$(basename $raw_reads _fulllength_rescued_merged.fastq.gz)

mkdir ${base}

files=$(awk '{printf (NR==1 ? "" : ",") $0}' ../raw_reads.txt)
labels="Twt1,Twt2,Twt3,Trog1,Trog2,Trog3,Dwt1,Dwt2,Dwt3,Drog1,Drog2,Drog3"

## use pre-processed fastq (adapter, length, and strand oriented)

## cluster
time /home/users/allstaff/yan.a/yan.a/software/RATTLE/rattle cluster -i ${files} -l ${labels} -t ${num_cores} -o ${base} --iso --rna

## cluster summary
time /home/users/allstaff/yan.a/yan.a/software/RATTLE/rattle cluster_summary -i ${files} -l ${labels} -c ${base}/clusters.out > ${base}/cluster_summary.tsv

## correct (min reads 5 default)
time /home/users/allstaff/yan.a/yan.a/software/RATTLE/rattle correct -i ${files} -l ${labels} -c ${base}/clusters.out -o ${base}/ -t ${num_cores} -r 5 

## polish (stranded cDNA data enable --rna)
time /home/users/allstaff/yan.a/yan.a/software/RATTLE/rattle polish -i ${base}/consensi.fq -l ${labels} -o ${base} -t ${num_cores} --rna --summary 

## some cleaning and formating
time ~/yan.a/software/seqtk/seqtk seq -a ${base}/transcriptome.fq | awk '/^[^>]/{ gsub(/U/,"T") }1' > ${base}/transcriptome.fasta

module purge
module load anaconda3/2019.03

conda activate /stornext/Bioinf/data/lab_davidson/yan.a/conda_env/isoncorrect
seqkit fx2tab -nlgGH ${base}/transcriptome.fasta > fx2tab.txt
