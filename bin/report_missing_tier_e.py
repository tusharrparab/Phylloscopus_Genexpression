#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--species-plan", required=True)
    parser.add_argument("--ortholog-targets", required=True)
    parser.add_argument("--outdir", required=True)
    return parser.parse_args()


def read_tsv(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows, fieldnames):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    species_rows = read_tsv(Path(args.species_plan))
    target_rows = read_tsv(Path(args.ortholog_targets))

    missing_species = []
    status_rows = []

    for species in species_rows:
        if species["evidence_tier"].strip() != "E":
            continue

        missing_species.append(
            {
                "species_id": species["species_id"].strip(),
                "scientific_name": species["scientific_name"].strip(),
                "reason": "no_usable_public_sequence_input",
            }
        )

        for target in target_rows:
            status_rows.append(
                {
                    "species_id": species["species_id"].strip(),
                    "scientific_name": species["scientific_name"].strip(),
                    "gene_id": target["gene_id"].strip(),
                    "gene_symbol": target["gene_symbol"].strip(),
                    "evidence_tier": "E",
                    "reconstruction_status": "missing_data",
                    "sequence_length": "0",
                    "confidence_tier": "none",
                    "method": "unavailable",
                    "notes": "no assembly, transcriptome, or usable raw read source in manifest",
                }
            )

    write_tsv(outdir / "missing_species.tsv", missing_species, ["species_id", "scientific_name", "reason"])
    write_tsv(
        outdir / "ortholog_status.tsv",
        status_rows,
        [
            "species_id",
            "scientific_name",
            "gene_id",
            "gene_symbol",
            "evidence_tier",
            "reconstruction_status",
            "sequence_length",
            "confidence_tier",
            "method",
            "notes",
        ],
    )


if __name__ == "__main__":
    main()

