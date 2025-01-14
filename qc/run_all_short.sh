#!/bin/bash
# Usage: sbatch slurm-serial-job-script
# Prepared By: Alex Yan
#              yan.a@wehi.edu.au

# NOTE: To activate a SLURM option, remove the whitespace between the '#' and 'SBATCH'

# To give your job a name, replace "MyJob" with an appropriate name
#SBATCH --job-name=transrate

# To set a project account for credit charging,
# SBATCH --account=ls25

# Request CPU resource for a serial job
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48

# Memory usage (MB)
#SBATCH --mem-per-cpu=4000

# Set your minimum acceptable walltime, format: day-hours:minutes:seconds
#SBATCH --time=48:00:00

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


## set up input and output directory prefix
# raw_reads='/vast/scratch/users/yan.a/vast_scratch/dRNA/A549_MCF7_directRNA_merged.fastq.gz'
read1='/vast/projects/davidson_longread/yan.a/Dong_2023/subsampleshort/Ilu_merged_R1.fq.gz'
read2='/vast/projects/davidson_longread/yan.a/Dong_2023/subsampleshort/Ilu_merged_R2.fq.gz'

num_cores=48
base="merged_ilu"

assembly="${base}/${base}_trinity.Trinity.fasta" # can be fq or fa file
stranded="stranded"

refgtf="/vast/projects/lab_davidson/yan.a/Dong_2023/reference/annotation_gencodev44_sequin.gtf"
reffa="/vast/projects/lab_davidson/yan.a/Dong_2023/reference/transcriptome_gencodev44_sequin.fa"
genome="/vast/projects/lab_davidson/yan.a/Dong_2023/reference/genome_GRCh38_sequin.fa"
busco="/home/users/allstaff/yan.a/lab_davidson/yan.a/dRNA/rattle/busco_downloads/lineages/primates_odb10"
buscodb=$(basename $busco)

# setting SQANTI3 enviroment
SQANTI="/home/users/allstaff/yan.a/lab_davidson/yan.a/Dong_2023/SQANTI3-5.1.2/" ## use modified version
cupcake="$HOME/yan.a/software/cDNA_Cupcake/"

echo Step1: Transrate

module purge
module load miniconda3
conda activate /stornext/Bioinf/data/lab_davidson/yan.a/conda_env/transrate

# # in case it is a dRNA assembly, do cleaning converting U to T
~/yan.a/software/seqtk/seqtk seq -a $assembly | awk '/^[^>]/{ gsub(/U/,"T") }1' > ${base}/transcriptome.reformat.fasta

assembly="${PWD}/${base}/transcriptome.reformat.fasta"
transrate --assembly $assembly --threads $num_cores --output ${base}_transrate --reference $reffa

echo Step2: BUSCO

module purge
module load miniconda3
conda activate /stornext/Bioinf/data/lab_davidson/yan.a/conda_env/busco

busco -i $assembly -l $busco -o ${base}_busco -m transcriptome -c $num_cores -f --offline # force overwrite

tar -czf ${base}_busco/run_${buscodb}.tar.gz ${base}_busco/run_${buscodb}
rm -rf ${base}_busco/run_${buscodb} 

echo Step3: SQANTI3

module purge
module load miniconda3
conda activate /stornext/Bioinf/data/lab_davidson/yan.a/conda_env/squanti3

export PYTHONPATH=$PYTHONPATH:${cupcake}
export PYTHONPATH=$PYTHONPATH:${cupcake}/sequence

rm -rf ${base}_sqanti3 # incase already generated corrected.fasta been used

mkdir -p ${base}_sqanti3
cd ${base}_sqanti3

echo starting first pass

python ${SQANTI}/sqanti3_qc.py $assembly $refgtf $genome \
                     --CAGE_peak ${SQANTI}/data/ref_TSS_annotation/human.refTSS_v3.1.hg38.bed    \
                     --polyA_motif_list ${SQANTI}/data/polyA_motifs/mouse_and_human.polyA_motif.txt    \
                     -o firstpass -d firstpass \
                     --cpus $num_cores --report skip --fasta --force_id_ignore --skipORF

# python ~/yan.a/software/rnaseq_utils/scripts/extract_sqanti_summary.py ${base}_sqanti3/${base}_SQANTI3_report.html > sqantisummary.txt

if [ $stranded != "stranded" ]
then
# if sequence data is unstranded (or unknown strandness)
# get all transcripts annotated as antisense, and RT, then do a sqanti3 again

    awk '$6=="antisense" {print $1}' firstpass/firstpass_classification.txt | sed 's/^/"/;s/$/"/' > asid.txt

    grep -F -f asid.txt firstpass/firstpass_corrected.gtf | awk 'BEGIN{FS=OFS="\t"} {if($7 == "+") $7="-"; else if($7 == "-") $7="+"; print}' > antisense.gtf

    echo starting second pass for antisense transcripts

    python ${SQANTI}/sqanti3_qc.py antisense.gtf $gtf $genome \
                        --CAGE_peak ${SQANTI}/data/ref_TSS_annotation/human.refTSS_v3.1.hg38.bed    \
                        --polyA_motif_list ${SQANTI}/data/polyA_motifs/mouse_and_human.polyA_motif.txt    \
                        -o secondpass -d secondpass \
                        --cpus $num_cores --report skip --force_id_ignore --skipORF

# add a code to avoid overlap novelGene in sense and antisense
    awk 'BEGIN{FS=OFS="\t"} {if($7 ~ /novelGene/) $7 = $7 "_as"; print}' secondpass/secondpass_classification.txt > secondpass/secondpass_classification_tmp.txt
# replace tx in sense with results from antisense run 
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$0; next} $1 in a{$0=a[$1]} {print $0}' secondpass/secondpass_classification_tmp.txt firstpass/firstpass_classification.txt > unstranded_classification.txt

fi
cd ..

echo Step4: BUSCO corrected

module purge
module load miniconda3
conda activate /stornext/Bioinf/data/lab_davidson/yan.a/conda_env/busco

# use first pass fasta, second pass fasta will just be a subset and reverse complement
# umapped tx will be ignored, and chimeric tx will be splitted with same header, duplcated header will be add a suffix
awk '/^>/ {seen[$0]++; if(seen[$0] > 1) {print $0 "_" seen[$0]-1} else {print}; next} {print}' ${base}_sqanti3/firstpass/firstpass_corrected.fasta > ${base}_sqanti3/${base}_corrected_dedup.fasta

busco -i ${base}_sqanti3/${base}_corrected_dedup.fasta -l $busco -o ${base}_busco_corrected -m transcriptome -c $num_cores -f

tar -czf ${base}_busco_corrected/run_${buscodb}.tar.gz ${base}_busco_corrected/run_${buscodb}
rm -rf ${base}_busco_corrected/run_${buscodb} 

echo Step5: Minimap2 and Salmon quantification merged 

module purge
module load miniconda3
conda activate /stornext/Bioinf/data/lab_davidson/yan.a/conda_env/isoncorrect
minimap2 --version
# 2.26-r1175
module load salmon/1.10.2
salmon --version
module load samtools/1.19.2 
module load bwa/0.7.17

## use raw fasta from output, not SQANTI corrected

tx="${PWD}/${base}/transcriptome.reformat.renamed.fasta"
method=$(basename ${PWD})

mkdir ${base}_salmon
cd ${base}_salmon

# bwa-mem and quant
bwa index $tx
bwa mem -t $num_cores $tx $read1 $read2 | samtools view -@ $num_cores -Sb > ${base}_${method}_unsorted.bam

samtools sort -@ $num_cores ${base}_${method}_unsorted.bam > ${base}_${method}_sorted.bam
samtools index -@ $num_cores ${base}_${method}_sorted.bam
samtools flagstat ${base}_${method}_sorted.bam > ${base}_${method}.flagstat

samtools view -@ $num_cores -f 2 -F3840 -Sb ${base}_${method}_unsorted.bam > ${base}_${method}_filtered.bam

salmon quant -p $num_cores -t $tx -l A -a ${base}_${method}_filtered.bam -o ${base}_quant_${method}_align
rm ${base}_${method}_sorted.bam ${base}_${method}_sorted.bam.bai ${base}_${method}_unsorted.bam ${base}_${method}_filtered.bam 

# salmon quant directly

salmon index -t $tx -i salmon_index -p $num_cores
salmon quant -i salmon_index -l A -1 $read1 -2 $read2 -o ${base}_quant_${method}_map -p $num_cores --gcBias --seqBias --posBias --thinningFactor 64
# rm read1.fastq read2.fastq 
rm -r salmon_index
cd ..

echo Step6: Corset clustering

module purge
module load miniconda3
conda activate /stornext/Bioinf/data/lab_davidson/yan.a/conda_env/isoncorrect
module load samtools/1.19.2 

minimap2 --version

## first pass
mkdir ${base}_corset
cd ${base}_corset

# echo $file

time minimap2 -t $num_cores -x asm5 -a -L -P --no-long-join -r2k -s40 -I 1000G $tx $tx | samtools sort -@ $num_cores -O BAM -o ${base}.bam

samtools index ${base}.bam
samtools idxstats ${base}.bam | cut -f1 | grep -v '^*' > refid.txt

# time samtools flagstat ${base}.bam > ${base}.flagstat
# samtools view -f 4 -F 256 ${base}.bam | cut -f 1 | awk -v OFS='\t' '{print $1, "nomap_" NR, 0}' > ${base}_unmapped_table.txt

time ~/yan.a/software/corset-1.09-linux64/corset -f true -m 1 -r true -p $base ${base}.bam

comm -13 <(cut -f1 ${base}-clusters.txt  | sort) <(sort refid.txt) | awk -v OFS='\t' '{print $1, "nomap_" NR, 0}' > ${base}_unmapped_table.txt

cat ${base}-clusters.txt <(cut -f1,2 ${base}_unmapped_table.txt) > ${base}-clusters_mod.txt
cd ..


echo Step7: Quantification of individual samples

module purge
module load miniconda3
conda activate /stornext/Bioinf/data/lab_davidson/yan.a/conda_env/isoncorrect
minimap2 --version
# 2.26-r1175
module load salmon/1.10.2
salmon --version
module load samtools/1.19.2
samtools --version 
module load bwa/0.7.17

mkdir -p ${base}_dge
cd ${base}_dge

bwa index $tx 
salmon index -t $tx -i salmon_index -p $num_cores

for i in {1..6}
do
read1=$(sed "${i}q;d" ../../raw_reads_short.txt)
name=$(basename $read1 _R1_sub.fastq.gz)
read2=$(echo $read1 | sed 's/R1/R2/g')

bwa mem -t $num_cores $tx $read1 $read2 | samtools view -@ $num_cores -Sb > ${name}_${method}_unsorted.bam

# samtools sort -@ 48 ${name}_${method}_unsorted.bam > ${name}_${method}_sorted.bam
samtools view -@ $num_cores -f 2 -F3840 -Sb ${name}_${method}_unsorted.bam > ${name}_${method}_filtered.bam

salmon quant -p $num_cores -t $tx -l A --numBootstraps 100 -a ${name}_${method}_filtered.bam -o ${name}_quant_${method}_align
rm  ${name}_${method}_unsorted.bam ${name}_${method}_filtered.bam 

# salmon quant directly

salmon quant -i salmon_index -l A --numBootstraps 100 -1 $read1 -2 $read2 -o ${name}_quant_${method}_map -p $num_cores --gcBias --seqBias --posBias --thinningFactor 64

done

rm -r salmon_index

