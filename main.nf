#!/usr/bin/env nextflow
nextflow.enable.dsl=2

log.info """
         P L A N T  R N A F M  P R E P R O C E S S I N G
         ===================================================
         Input file: ${params.input}
         Output dir: ${params.outdir}
         """

// --- Main Workflow ---
workflow {
    // Channel structure: [species_id, transcriptome_path, proteome_path]
    Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row -> tuple(row.species_id, file(row.transcriptome), file(row.proteome)) }
        .set { species_ch }

    mmseqs_search(species_ch)
    post_process(mmseqs_search.out)
}

// --- Processes ---

process mmseqs_search {
    label 'cpu_heavy'
    errorStrategy 'ignore'
    tag "Running mmseqs search for ${species_id}"
    publishDir "${params.outdir}/mmseqs_results", mode: 'copy', pattern: "${species_id}_results.tsv"

    input:
    tuple val(species_id), path(transcriptome), path(proteome)

    output:
    tuple val(species_id), path("${species_id}_results.tsv"), path(transcriptome)

    script:
    """
    mmseqs easy-search ${proteome} ${transcriptome} ${species_id}_results.tsv tmp \
            --search-type 3 \
            --format-output "qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore" \
            --threads ${task.cpus}
    """
}

process post_process {
    label 'cpu_medium'
    tag "Processing alignments for ${species_id}"
    publishDir "${params.outdir}/labeled_sequences", mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(species_id), path(mmseqs_file), path(transcriptome_fasta)

    output:
    path("${species_id}.parquet")

    script:
    """
    process_alignments.py \
        --tblastn ${mmseqs_file} \
        --fasta ${transcriptome_fasta} \
        --output ${species_id}.parquet
    """
}