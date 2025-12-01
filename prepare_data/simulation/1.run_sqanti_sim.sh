#!/bin/bash
module load anaconda3/2019.03

conda activate /stornext/Bioinf/data/lab_davidson/yan.a/conda_env/squanti3

SQANTI="$HOME/davidson_longread/yan.a/simulation_20240501/SQANTI3-5.1.2/" ## use modified version
cupcake="$HOME/yan.a/software/cDNA_Cupcake/"

export PYTHONPATH=$PYTHONPATH:${cupcake}
export PYTHONPATH=$PYTHONPATH:${cupcake}/sequence

export PATH="$PATH:/vast/projects/lab_davidson/yan.a/software/SQANTI-SIM/"

num_cores=48

gtf="${PWD}/gencode.v44.subset.gtf"
genome="/vast/projects/lab_davidson/yan.a/ref/gencode/GRCh38.primary_assembly.genome.fa"

## 1, design sqanti-index

mkdir -p design
cd design

# classify
sqanti-sim.py classif \
    --gtf $gtf \
    --cores ${num_cores}

# design
sqanti-sim.py design sample \
    -i sqanti-sim_index.tsv \
    --gtf $gtf \
    --genome $genome\
    --ont \
    --trans_number 45000 \
    --ISM 2000 --NIC 2000 --NNC 2000 \
    --seed 1234 \
    --cores ${num_cores} \
    --read_type cDNA

## 2, select DE 

## select diff gene and transcripts, change expression, save baseline in 2 conditions, and true DE/DT
Rscript --vanilla ../get_diff.R

## simulate sqanti-sim-index files for replicates using Gamma-Poission
Rscript --vanilla ../simulate_dispersion.R

## 3, sqanti-sim simulation

for J in {1..3}
do

# simulate ctrl 
mkdir -p ctrl${J} 

cp ../design/sqanti-sim_index_ctrl${J}.tsv ctrl${J}/sqanti-sim_index.tsv

cd ctrl${J} 
# simulation
sqanti-sim.py sim \
    -i sqanti-sim_index.tsv \
    --gtf $gtf \
    --genome $genome \
    --read_type cDNA \
    --ont \
    --long_count 1000000 \
    --illumina \
    --short_count 15000000 \
    --seed ${J} \
    --cores ${num_cores}

cd ..

# simulate de 
mkdir -p de${J} 

cp ../design/sqanti-sim_index_de${J}.tsv de${J}/sqanti-sim_index_de.tsv

cd de${J} 
# simulation
sqanti-sim.py sim \
    -i sqanti-sim_index_de.tsv \
    --gtf $gtf \
    --genome $genome \
    --read_type cDNA \
    --ont \
    --long_count 1000000 \
    --illumina \
    --short_count 15000000 \
    --seed ${J} \
    --cores ${num_cores}

cd ..

done
