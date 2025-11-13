#!/usr/bin/env python3

import argparse
from pathlib import Path

import polars as pl
import polars.selectors as cs
import pyfastx


def process_alignments(tblastn_output: Path, transcriptome_fasta: Path, proteome_fasta: Path, output_parquet: Path):
    """Processes mmseqs output to create a labeled nucleotide sequence dataframe."""

    # Define the schema for the mmseqs output format
    mmseqs_schema = {
        "query": pl.Utf8,
        "target": pl.Utf8,
        "pident": pl.Float64,
        "alnlen": pl.Int64,
        "mismatch": pl.Int64,
        "gapopen": pl.Int64,
        "qstart": pl.Int64,
        "qend": pl.Int64,
        "tstart": pl.Int64,
        "tend": pl.Int64,
        "evalue": pl.Float64,
        "bits": pl.Float64,
    }

    # Read the mmseqs output
    df = pl.read_csv(
        tblastn_output,
        has_header=False,
        separator="\t",
        new_columns=list(mmseqs_schema.keys()),
        schema=mmseqs_schema,
    )
    # Filter to keep only self-hits (where protein matches its own source transcript)
    df = df.filter(pl.col("query") == pl.col("target"))
    # Get the best hit for each transcript (target)
    best_hits = df.group_by("target").agg(pl.col("bits").arg_max()).select(
        pl.col("target"), pl.col("bits").alias("best_hit_idx")
    )
    grouped = df.group_by("target").agg(pl.col("*"))
    good_hits = grouped.join(best_hits, on='target').with_columns(cs.list().list.get(pl.col("best_hit_idx")))

    # Load the transcriptome and proteome
    transcripts = {record.name : record.seq for record in pyfastx.Fasta(str(transcriptome_fasta))}
    proteins    = {record.name: record.seq  for record in pyfastx.Fasta(str(proteome_fasta))}

    records = []
    for row in good_hits.iter_rows(named=True):
        transcript_id = row["target"]
        if transcript_id not in transcripts:
            continue

        transcript_seq = transcripts[transcript_id]
        seq_len = len(transcript_seq)

        # mmseqs gives 1-based coordinates
        start = row["tstart"]
        end = row["tend"]

        # Ensure start is always less than end
        if start > end:
            start, end = end, start

        start -= 1  # Convert to 0-based

        # Generate the label sequence
        label_seq = (
            "5" * start
            + "C" * (end - start)
            + "3" * (seq_len - end)
        )

        # Check if protein starts with Methionine
        starts_with_met = False
        if transcript_id in proteins:
            protein_seq = proteins[transcript_id]
            starts_with_met = protein_seq[0] == 'M' if len(protein_seq) > 0 else False

        records.append(
            {
                "transcript_id": transcript_id,
                "nucleotide_sequence": transcript_seq,
                "nucleotide_label": label_seq,
                "starts_with_methionine": starts_with_met,
            }
        )

    # Create the final dataframe and write to parquet
    if records:
        output_df = pl.DataFrame(records)
        output_df.write_parquet(output_parquet)


def main():
    parser = argparse.ArgumentParser(description="Process mmseqs output and generate labeled sequences.")
    parser.add_argument("--tblastn", type=Path, required=True, help="Path to the mmseqs output file.")
    parser.add_argument("--fasta", type=Path, required=True, help="Path to the transcriptome FASTA file.")
    parser.add_argument("--proteome", type=Path, required=True, help="Path to the proteome FASTA file.")
    parser.add_argument("--output", type=Path, required=True, help="Path to the output Parquet file.")
    args = parser.parse_args()

    process_alignments(args.tblastn, args.fasta, args.proteome, args.output)


if __name__ == "__main__":
    main()
