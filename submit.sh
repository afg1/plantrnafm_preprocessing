#!/bin/bash

#SBATCH --job-name=plantrnafm_preprocessing
#SBATCH --output=out_plrnafm
#SBATCH --error=err_plrnafm
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G

# Load required modules
module load nextflow

# Run the Nextflow pipeline
nextflow run -resume main.nf 
