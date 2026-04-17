#!/usr/bin/env python3

import argparse
import csv
import shutil
from pathlib import Path


REQUIRED_COLUMNS = {
    "species_manifest": ["species_id", "scientific_name"],
    "ortholog_targets": ["gene_id", "gene_symbol", "ref_species", "cds_length"],
    "reference_manifest": ["reference_id", "scientific_name", "reference_role"],
}

UNIQUE_COLUMNS = {
    "species_manifest": "species_id",
    "ortholog_targets": "gene_id",
    "reference_manifest": "reference_id",
}


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--species-manifest", required=True)
    parser.add_argument("--ortholog-targets", required=True)
    parser.add_argument("--reference-manifest", required=True)
    parser.add_argument("--outdir", required=True)
    return parser.parse_args()


def read_rows(path: Path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        rows = list(reader)
    return fieldnames, rows


def validate_table(label: str, path: Path):
    fieldnames, rows = read_rows(path)
    missing = [column for column in REQUIRED_COLUMNS[label] if column not in fieldnames]
    if missing:
        raise SystemExit(f"{label} is missing required columns: {', '.join(missing)}")

    for index, row in enumerate(rows, start=2):
        if None in row:
            raise SystemExit(
                f"{label} has too many columns on line {index}; check delimiter alignment in {path.name}"
            )

    unique_column = UNIQUE_COLUMNS[label]
    seen = set()
    for row in rows:
        key = (row.get(unique_column) or "").strip()
        if not key:
            raise SystemExit(f"{label} contains an empty value for {unique_column}")
        if key in seen:
            raise SystemExit(f"{label} contains a duplicate value for {unique_column}: {key}")
        seen.add(key)

    return rows


def write_summary(outdir: Path, counts):
    summary_path = outdir / "validation_summary.tsv"
    with summary_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["file", "records", "status"])
        for name, count in counts:
            writer.writerow([name, count, "ok"])


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    sources = {
        "species_manifest": Path(args.species_manifest),
        "ortholog_targets": Path(args.ortholog_targets),
        "reference_manifest": Path(args.reference_manifest),
    }

    counts = []
    for label, path in sources.items():
        rows = validate_table(label, path)
        counts.append((path.name, len(rows)))

    shutil.copyfile(sources["species_manifest"], outdir / "species_manifest.tsv")
    shutil.copyfile(sources["ortholog_targets"], outdir / "ortholog_targets.tsv")
    shutil.copyfile(sources["reference_manifest"], outdir / "reference_manifest.tsv")
    write_summary(outdir, counts)


if __name__ == "__main__":
    main()
