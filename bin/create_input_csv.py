import csv
import argparse
from pathlib import Path
import sys

def create_input_sheet(assemblies_dir: Path, output_csv: Path):
    """
    Scans a directory of assemblies and creates a CSV file mapping species
    to their transcript and protein FASTA files.
    """
    if not assemblies_dir.is_dir():
        print(f"Error: Directory not found: {assemblies_dir}", file=sys.stderr)
        sys.exit(1)

    header = ["species_id", "transcriptome", "proteome"]
    with open(output_csv, "w", newline="") as f_out:
        writer = csv.writer(f_out)
        writer.writerow(header)

        for species_dir in sorted(assemblies_dir.iterdir()):
            if not species_dir.is_dir():
                continue

            species_id = species_dir.name.split("-")[0]
            transcript_file = next(species_dir.glob("*-Trans-assembly.fa"), None)
            protein_file = next(species_dir.glob("*-translated-protein.fa"), None)

            if transcript_file and protein_file:
                writer.writerow([species_id, str(transcript_file.resolve()), str(protein_file.resolve())])
            else:
                print(f"Warning: Missing files for {species_dir.name}", file=sys.stderr)

def main():
    parser = argparse.ArgumentParser(description="Create an input CSV for the Nextflow pipeline.")
    parser.add_argument(
        "-i", "--input-dir",
        type=Path,
        default=Path("assemblies"),
        help="Path to the directory containing species subdirectories."
    )
    parser.add_argument(
        "-o", "--output-csv",
        type=Path,
        default=Path("input.csv"),
        help="Path to write the output CSV file."
    )
    args = parser.parse_args()

    create_input_sheet(args.input_dir, args.output_csv)
    print(f"Successfully created input sheet: {args.output_csv}")

if __name__ == "__main__":
    main()