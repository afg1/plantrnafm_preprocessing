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

    makeblastdb(species_ch)
    tblastn(makeblastdb.out.db_ch)
    post_process(tblastn.out.tblastn_ch.join(makeblastdb.out.fasta_ch))
}

// --- Processes ---

process makeblastdb {
    label 'cpu_light'
    tag "Make BLAST DB for ${species_id}"
    publishDir "${params.outdir}/blast_db", mode: 'copy', pattern: "${species_id}_db.*"

    input:
    tuple val(species_id), path(transcriptome), path(proteome)

    output:
    // Pass database files and proteome to the next step
    tuple val(species_id), path("${species_id}_db.*"), path(proteome), emit: 'db_ch'
    // Pass transcriptome fasta for post-processing step
    tuple val(species_id), path(transcriptome), emit: 'fasta_ch'

    script:
    """
    makeblastdb -in ${transcriptome} -dbtype nucl -out ${species_id}_db
    """
}

process tblastn {
    label 'cpu_heavy'
    tag "Running tblastn for ${species_id}"
    publishDir "${params.outdir}/tblastn_results", mode: 'copy'

    input:
    tuple val(species_id), path(db), path(proteome)

    output:
    tuple val(species_id), path("${species_id}_tbn_results.txt"), emit: 'tblastn_ch'

    script:
    """
    tblastn -query ${proteome} \
            -db ${db[0].simpleName} \
            -out ${species_id}_tbn_results.txt \
            -outfmt 6 \
            -evalue 1e-10 \
            -qcov_hsp_perc 70.0
    """
}

process post_process {
    label 'cpu_medium'
    tag "Processing alignments for ${species_id}"
    publishDir "${params.outdir}/labeled_sequences", mode: 'copy'

    input:
    tuple val(species_id), path(tblastn_file), path(transcriptome_fasta)

    output:
    path("${species_id}.parquet")

    script:
    """
    process_alignments.py \
        --tblastn ${tblastn_file} \
        --fasta ${transcriptome_fasta} \
        --output ${species_id}.parquet
    """
}