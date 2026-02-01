# Long-Read Simulation Pipeline

This pipeline prepares reference data and performs long-read simulation using `SQANTI-SIM`. This directory contains the scripts required to simulate Long-Read transcriptomic data with controlled Differential Gene Expression (DGE), Differential Transcript Expression (DTE), and Differential Transcript Usage (DTU). It is organized into two sequential steps.

## Pipeline Steps

### Step 1: Reference Preparation
**Script:** `1.run_ref_prep.sh`

This steps performs the following:
1.  **Alignment & Quantification:** Aligns raw GTEx long-read data to the reference transcriptom using `minimap2` and quantifies abundance using `salmon`.
2.  **Baseline Estimation:** Runs `gtex_info.R` to generate `baselineAbundance.rds` and a list of expressed transcripts (`txid.txt`).
3.  **Subsetting:** specific GTF, FASTA, and BED files for the identified transcripts.

**Usage:**
```bash
sbatch 1.run_ref_prep.sh
# Ensure this completes successfully before running Step 2
```

### Step 2: Simulation
**Script:** `2.run_simulation.sh`

This step performs the actual simulation:
1.  **Design:** Runs `sqanti-sim.py design` to create the simulation index.
2.  **Differential Expression:** Calls `get_diff.R` to establish DGE/DTU/DTE ground truth lists.
3.  **Dispersion:** Calls `simulate_dispersion.R` to simulate biological variation across replicates.
4.  **Run Simulation:** Executes `sqanti-sim.py sim` to generate synthetic FASTQ reads for Control and DE conditions.

**Usage:**
```bash
sbatch 2.run_simulation.sh
```

## Dependencies
*   **Input Data:** Expects raw data at `../../dataset/GTeX/long-read/sequence_data/`
*   **Reference Files:**
    *   **Transcript Reference FASTA:** `gencode.v44.transcripts.fa` (Required for Step 1)
    *   **Genome Reference FASTA:** `GRCh38.primary_assembly.genome.fa` (Required for Step 2)
    *   **Annotation GTF:** `gencode.v44.annotation.gtf`
*   **Software & Versions:**
    *   `minimap2` (2.26-r1175)
    *   `salmon` (1.10.2)
    *   `samtools` (1.19.2)
    *   `seqkit` (2.5.1)
    *   `SQANTI3` (5.1.2)
    *   `SQANTI-SIM` (0.2.1)
