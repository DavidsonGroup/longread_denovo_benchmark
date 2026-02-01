#!/bin/bash
#SBATCH --job-name=ref_prep
#SBATCH --account=ls25
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --mem-per-cpu=8000
#SBATCH --time=48:00:00
#SBATCH --mail-user=yan.a@wehi.edu.au
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=BEGIN
#SBATCH --output=stdout_prep
#SBATCH --error=stderr_prep
#SBATCH --partition=genomics
#SBATCH --qos=genomics

module purge
module load anaconda3/2019.03
conda activate /stornext/Bioinf/data/lab_davidson/yan.a/conda_env/isoncorrect
module load salmon/1.10.2
module load samtools/1.19.2

num_cores=48
# Using the vast path for reference as seen in script 3
ref='/vast/projects/lab_davidson/yan.a/ref/gencode/gencode.v44.transcripts.fa' 
annotation_gtf='/vast/projects/lab_davidson/yan.a/ref/gencode/gencode.v44.annotation.gtf'

minimap2 -x map-ont -I 1000G -t $num_cores -d index.mmi $ref

for raw_reads in ../../dataset/GTeX/long-read/sequence_data/*gz
do
    name=$(basename $raw_reads .fastq.gz)
    
    minimap2 -ax map-ont -Y -p 1.0 -N 100 -t $num_cores index.mmi $raw_reads | samtools view -@ $num_cores -Sb > ${name}_${method}_unsorted.bam
    salmon quant --ont -p $num_cores -t $ref -l A --numBootstraps 100 -a ${name}_${method}_unsorted.bam -o ${name}_quant_${method}_onts
    
    rm ${name}_${method}_unsorted.bam 
done

Rscript gtex_info.R

# extract gtf
grep -Ff txid.txt $annotation_gtf > gencode.v44.subset.gtf 

# extract fasta and modify the header
seqkit grep -r -f txid.txt $ref | awk 'BEGIN{FS=OFS="|"} /^>/{print $1; next} {print}' > gencode.v44.subset.fasta

# convert to bed
paftools.js gff2bed gencode.v44.subset.gtf > gencode.v44.subset.bed 
