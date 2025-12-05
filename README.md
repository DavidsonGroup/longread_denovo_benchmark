- [Longread\_denovo\_benchmark](#longread_denovo_benchmark)
  - [Generate simulation data](#generate-simulation-data)
  - [Assemble](#assemble)
    - [ONT simulation data (unstraned)](#ont-simulation-data-unstraned)
    - [ONT PCR-cDNA data from cancer cell lines (stranded)](#ont-pcr-cdna-data-from-cancer-cell-lines-stranded)
    - [ONT Direct RNA data from cancer cell lines (stranded)](#ont-direct-rna-data-from-cancer-cell-lines-stranded)
    - [PacBio kinnex Human PBMC single-cell data 10x 3' kit (stranded)](#pacbio-kinnex-human-pbmc-single-cell-data-10x-3-kit-stranded)
    - [ONT Pea PCR-cDNA data (stranded)](#ont-pea-pcr-cdna-data-stranded)
  - [Generate summary and quality metrics](#generate-summary-and-quality-metrics)
  - [Summarise metric and differential analysis in R](#summarise-metric-and-differential-analysis-in-r)
  - [Generate plot](#generate-plot)
  - [Reference](#reference)
  - [License](#license)


# Longread_denovo_benchmark
Code to generate, process and quality check long read *de novo* transcriptome assembly

## Generate simulation data
[Code for simulation](prepare_data/simulation/)

We first obtained a subset of transcripts that are widely expressed in the GTEx v9 long read dataset (92 samples) using Gencode comprehensive annotation (v44). We kept transcripts with more than 5 reads in at least 15 samples after Salmon quantification (18145 genes, 40509 transcripts), and stored their mean count per million (CPM) values as the control group’s baseline expression. 

We then generated a perturbed set of CPM values where transcript expression was changed by: (1) randomly selecting 1000 genes and changing all transcripts belonging to that gene concordantly (500 genes 2 fold up and 500 genes 2 fold down), (2) selected another 1000 genes randomly, and then select 2 random transcripts from the gene and swap their expression, (3) selected another 1000 genes randomly, and then select 1 random transcript to change its expression (500 transcripts 2 fold up and 500 transcripts 2 fold down). The updated CPM were stored as the perturbed group baseline expression. 

We then generated a count matrix and CPM matrix for 3 control replicates and 3 perturbed replicates with gamma distribution, followed by a Poisson distribution (Baldoni et al., [2024](https://doi.org/10.1093/nar/gkad1167)). Both long-read and short-read FASTQ files were simulated using SQANTI-SIM with default settings and ONT R9.4 cDNA error profile (v 0.2.1) (Mestre-Tomás et al., [2023](https://doi.org/10.1186/s13059-023-03127-0)). The long read data contained 6 million reads in total, and an average read length of 1085 bp, and short read data was 100 bp paired-end. We then subsampled the short-read data to match the total number of base pairs in the long read data (6.5 billion bases). 

The simulated data was non-stranded, and contains 2000 DE genes, 2000 genes with DTU, 5927 transcripts with DTU and 6933 DE transcripts. It is available for download at https://doi.org/10.5281/zenodo.14263456.  


## Assemble 
All assemblies are uploaded to https://doi.org/10.5281/zenodo.17538009.- [Longread\_denovo\_benchmark](#longread_denovo_benchmark)

### ONT simulation data (unstraned)
[Code for assembling](assemble/simulation/)

### ONT PCR-cDNA data from cancer cell lines (stranded)
[Code for assembling](assemble/pcr_cdna/)

### ONT Direct RNA data from cancer cell lines (stranded)
[Code for assembling](assemble/drna/)

### PacBio kinnex Human PBMC single-cell data 10x 3' kit (stranded)
Data downloaded and processed use [code](prepare_data/single_cell/)
[Code for assembling](assemble/single_cell/)

### ONT Pea PCR-cDNA data (stranded)
[Code for assembling](assemble/pea/)

## Generate summary and quality metrics
We now provide code to run all quanlity checks, this can also be done using the [nextflow pipeline denovo_qc_nextflow](https://github.com/alexyfyf/denovo_qc_nextflow). 

[Code for generating quality metrics](qc/)

We have generated the following measures to compare the quality of assemblies.
1. Transrate analysis of *de novo* transcriptome
2. BUSCO analysis using raw *de novo* transcriptome
3. SQANTI3 analysis of *de novo* transcriptome
4. BUSCO analysis using genome corrected *de novo* transcriptome
5. Salmon quantification of pooled samples using *de novo* transcriptome
6. (Optional) Oarfish quantification of pooled samples or single-cell using *de novo* transcriptome
7. Corset clustering of *de novo* transcriptome
8. Salmon quantification of individual samples using *de novo* transcriptome

## Summarise metric and differential analysis in R
[Code for summarising and DE analysis](R/)

We then performed DE analysis (DGE, DTE, DTU) using the count matrix from step 8 for bulk RNAseq, and from step 6 for single-cell analysis (clustering and pseudobulking first).

## Generate plot

## Reference
Towards accurate, reference-free differential expression: A comprehensive evaluation of long-read *de novo* transcriptome assembly.
Feng Yan, Pedro L. Baldoni, James Lancaster, Matthew E. Ritchie, Mathew G. Lewsey, Quentin Gouil, Nadia M. Davidson.
**bioRxiv** 2025.02.02.635999; doi: https://doi.org/10.1101/2025.02.02.635999

## License
This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
