## Overall Goal
We are reproducing the preprocessing pipeline from a recent publication which uses a large transcriptome dataset to build an RNA language model. Part of this preprocessing requires us to run tblastn between RNA transcripts and translated proteins to find the CDS, 5' UTR and 3' UTR regions in the transcribed RNA. 

## Specific approach
We have transcriptome data for just over 1400 plants, comprising the assembled transcripts, and the translated proteins. These are in FASTA files, and are stored in a directory structure like this:

├── assemblies
│   ├── AALA-Meliosma_cuneifolia
│   │   ├── AALA-SOAPdenovo-Trans-assembly.fa
│   │   ├── AALA-translated-protein.fa
│   ├── AAXJ-Atriplex_prostrata-2_samples_combined
│   │   ├── AAXJ-SOAPdenovo-Trans-assembly.fa
│   │   ├── AAXJ-translated-protein.fa
│   ├── ABCD-Racomitrium_elongatum
│   │   ├── ABCD-SOAPdenovo-Trans-assembly.fa
│   │   ├── ABCD-translated-protein.fa
...
(Note, filenames may not foillow this exact pattern, but the transcripts will aways have *Trans-assembly.fa, and the proteins will always have *-translated-protein.fa)

We aim to do a mapping per species, so for each of these directories we need to:

1. Run `makeblastdb -in {species transcript file} -dbtype nucl -out {species}_db`
2. Run `tblastn -query {species protein file} -db {species}_db - out {species}_tbn_results.txt -outfmt 6 -evalue 1e-10 -qcov_hsp_perc 70.0`

This will result in ~1400 results files, which contain the auery and subject IDs linking the nucleotide and protein sequences, and the start and end coordinates of the alignment, relative to the start of the transcript.

We will then need to process the alignment output to label each nucleotide in a nucleotide sequence with matching protein, such that all nucleotides prior to the start coordinate are labelled '5' (for 5' UTR), the nucleotides within the start and stop coordinate are labelled 'C' (for CDS) and the nucleotides after the stop coordinate are labelled '3' (for 3' UTR'). The output of this will therefore be a sequence something like 55555555CCCCCCCCCCCC33333. This should be written into a polars dataframe with three columns: the transcript ID (taken from column 2 of the tblastn output), the nucleotide sequence (taken from the transcript FASTA file), and the nucleotide label created from the output processing code. The dataframe should be written to disk in parquet format, one for each species, meaning we will have ~1400 parquet files with the 5'UTR-CDS-3'UTR data.

These parquet files will be the output of our pipeline.

## Orchestration
We will use nextflow for this project. We will have a process for making the blastdb, and another for running tblastn. These processes must run for each species individually. The pipeline will need to keep track of which output goes with which species. There are a few ways to do this, one of which may include some kind of pre-flight script that produces a CSV file linking species prefix (the four letters, e.g. AALA) to the full path of transcript and protein FASTA files, which is then taken as input to the pipeline. It is absolutely critical that the species files do not get mixed up, so be sure to use some kind of metadata passing in the processes to ensure that.

The output of the final process should copy the parquet files to a results directory.

## Additional Requirements
- The pipeline will be running on a SLURM cluster, so we will need to have a singularity container with the installed python environment. It will not be possible to build this container locally, so I suggest you write a github workflow to handle building it and pushing it to the ghcr.

- You will have to write both nextflow and python code to accomplish this task

- For the python post-processing, you should use the polars python package to handle the data processing. Use what you know about the output format of tblastn to specify the schema and read directly from the .txt file (using read_csv with delimiter='\t'). Write output in parquet format.

- This project has been initialised with uv, and we are using that for dependency management. There is a virtual environment available should you need to test anything. However, the data is not available for testing due to its size.

- Your (python) code should be checked by black