#!/usr/bin/env python3

import argparse
from pathlib import Path

import polars as pl
import pyfastx


def process_alignments(tblastn_output: Path, transcriptome_fasta: Path, output_parquet: Path):
    """Processes tblastn output to create a labeled nucleotide sequence dataframe."""

    # Define the schema for the tblastn output format 6
    tblastn_schema = {
        "qseqid": pl.Utf8,
        "sseqid": pl.Utf8,
        "pident": pl.Float64,
        "length": pl.Int64,
        "mismatch": pl.Int64,
        "gapopen": pl.Int64,
        "qstart": pl.Int64,
        "qend": pl.Int64,
        "sstart": pl.Int64,
        "send": pl.Int64,
        "evalue": pl.Float64,
        "bitscore": pl.Float64,
    }

    # Read the tblastn output
    df = pl.read_csv(
        tblastn_output,
        has_header=False,
        separator="\t",
        new_columns=list(tblastn_schema.keys()),
        schema=tblastn_schema,
    )

    # Get the best hit for each transcript (sseqid)
    best_hits = df.group_by("sseqid").agg(pl.col("bitscore").arg_max()).select(
        pl.col("sseqid"), pl.col("bitscore").alias("best_hit_idx")
    )
    df = df.join(best_hits, on="sseqid").filter(pl.col("bitscore").arg_max() == pl.col("best_hit_idx"))

    # Load the transcriptome
    transcripts = pyfastx.Fasta(str(transcriptome_fasta))

    records = []
    for row in df.iter_rows(named=True):
        transcript_id = row["sseqid"]
        if transcript_id not in transcripts:
            continue

        transcript_seq = transcripts[transcript_id].seq
        seq_len = len(transcript_seq)

        # tblastn gives 1-based coordinates
        start = row["sstart"]
        end = row["send"]

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

        records.append(
            {
                "transcript_id": transcript_id,
                "nucleotide_sequence": transcript_seq,
                "nucleotide_label": label_seq,
            }
        )

    # Create the final dataframe and write to parquet
    if records:
        output_df = pl.DataFrame(records)
        output_df.write_parquet(output_parquet)


def main():
    parser = argparse.ArgumentParser(description="Process tblastn output and generate labeled sequences.")
    parser.add_argument("--tblastn", type=Path, required=True, help="Path to the tblastn output file.")
    parser.add_argument("--fasta", type=Path, required=True, help="Path to the transcriptome FASTA file.")
    parser.add_argument("--output", type=Path, required=True, help="Path to the output Parquet file.")
    args = parser.parse_args()

    process_alignments(args.tblastn, args.fasta, args.output)


if __name__ == "__main__":
    main()
