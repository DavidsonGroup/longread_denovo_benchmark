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

# SBATCH --array=1-3

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
conda activate /stornext/Bioinf/data/lab_davidson/yan.a/conda_env/isoncorrect
minimap2 --version
# 2.26-r1175
module load salmon/1.10.2
salmon --version

module load samtools/1.19.2
samtools --version 

## use raw fasta from output, not SQANTI corrected
raw_reads='/vast/scratch/users/yan.a/vast_scratch/pea_fastq_merged/pea_fulllength_rescued_merged.fastq.gz'
num_cores=48
base=$(basename $raw_reads _fulllength_rescued_merged.fastq.gz)

gtf="/home/users/allstaff/yan.a/lab_davidson/yan.a/ref/ncbi/pea_zw6/GCF_024323335.1_CAAS_Psat_ZW6_1.0_genomic.gtf"
genome="/home/users/allstaff/yan.a/lab_davidson/yan.a/ref/ncbi/pea_zw6/GCF_024323335.1_CAAS_Psat_ZW6_1.0_genomic.fna"

# paftools.js gff2bed $gtf > /home/users/allstaff/yan.a/lab_davidson/yan.a/ref/ncbi/pea_zw6/anno.bed
junc="/home/users/allstaff/yan.a/lab_davidson/yan.a/ref/ncbi/pea_zw6/anno.bed"

mkdir -p ${base}_dge
cd ${base}_dge

# do indexing on the fly
#minimap2 -x map-ont -I 1000G -t $num_cores -d index.mmi $genome

for i in {1..12}
do
raw_reads=$(sed "${i}q;d" ../../raw_reads.txt)
name=$(basename $raw_reads | awk -F'_' '{print $1}')

## -uf for full-length cDNA
minimap2 -ax splice -I 1000G -t $num_cores --junc-bed $junc $genome $raw_reads | samtools sort -@ $num_cores -O BAM -o ${name}_sorted.bam
samtools index -@ $num_cores ${name}_sorted.bam

done

cd ..

# rm index.mmi

module purge
module load R/4.2.0

Rscript --vanilla bambu_analysis.R ${base}_dge