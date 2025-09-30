#!/bin/bash

#SBATCH --job-name=plantrnafm_preprocessing
#SBATCH --output=nextflow_%j.log
#SBATCH --error=nextflow_%j.err
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G

# Load required modules
module load nextflow

# Run the Nextflow pipeline
nextflow run main.nf 
