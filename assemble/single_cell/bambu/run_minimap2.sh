#!/bin/bash
# Usage: sbatch slurm-serial-job-script
# Prepared By: Alex Yan
#              yan.a@wehi.edu.au

# NOTE: To activate a SLURM option, remove the whitespace between the '#' and 'SBATCH'

# To give your job a name, replace "MyJob" with an appropriate name
#SBATCH --job-name=bambu

# To set a project account for credit charging,
# SBATCH --account=ls25

# Request CPU resource for a serial job
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48

# Memory usage (MB)
#SBATCH --mem-per-cpu=4000

# Set your minimum acceptable walltime, format: day-hours:minutes:seconds
#SBATCH --time=7:00:00

# To receive an email when job completes or fails
#SBATCH --mail-user=yan.a@wehi.edu.au
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=BEGIN

# BATCH --array=1-3

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

micromamba activate /stornext/Bioinf/data/lab_davidson/yan.a/conda_env/isoncorrect
minimap2 --version
# 2.26-r1175
module load salmon/1.10.2
salmon --version

module load samtools/1.19.2
samtools --version 

## use raw fasta from output, not SQANTI corrected

# raw_reads=$(sed "${J}q;d" ../raw_reads.txt)
num_cores=48
base='PB_merged'

method=$(basename ${PWD})
genome="/vast/projects/lab_davidson/yan.a/ref/gencode/GRCh38.primary_assembly.genome.fa"
gtf="/vast/projects/lab_davidson/yan.a/ref/gencode/gencode.v44.annotation.gtf"

paftools.js gff2bed $gtf > anno.bed
junc=$(realpath anno.bed)

mkdir -p ${base}_dge
cd ${base}_dge

# do indexing on the fly
#minimap2 -x map-ont -I 1000G -t $num_cores -d index.mmi $genome

for i in {1..3}
do
raw_reads=$(sed "${i}q;d" ../../raw_reads.txt)
# name=$(basename "$(dirname "$raw_reads")")
name=$(basename $raw_reads | awk -F'.' '{print $1}')

# -uf for full-length cDNA, use PacBio hq mode
minimap2 -ax splice:hq -I 1000G -uf -t 48 --junc-bed $junc $genome $raw_reads | samtools sort -@ 48 -O BAM -o ${name}_sorted.bam
samtools index -@ 48 ${name}_sorted.bam

done

cd ..

# rm index.mmi

module purge
module load R/4.2.3

Rscript --vanilla bambu_analysis.R ${base}_dge