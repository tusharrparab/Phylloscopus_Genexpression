#!/usr/bin/env python3

import argparse
import csv
import shutil
from pathlib import Path


REQUIRED_COLUMNS = {
    "species_manifest": ["species_id", "scientific_name"],
    "ortholog_targets": [
        "gene_id",
        "gene_symbol",
        "category",
        "orthology_basis",
        "copy_number_expectation",
        "ref_species",
        "cds_length",
        "rationale",
    ],
    "reference_manifest": ["reference_id", "scientific_name", "reference_role"],
}

UNIQUE_COLUMNS = {
    "species_manifest": "species_id",
    "ortholog_targets": "gene_id",
    "reference_manifest": "reference_id",
}

TARGET_CATEGORIES = {
    "phylogenetic_backbone",
    "migration_candidate",
    "vocalization_neural_candidate",
    "hypoxia_elevation_candidate",
    "housekeeping_control",
}

ORTHOLOGY_BASIS = {
    "historical_single_locus_marker",
    "single_copy_preferred",
    "candidate_gene_requires_validation",
    "expression_control_not_for_phylogeny",
}

COPY_NUMBER_EXPECTATIONS = {
    "single_copy_preferred",
    "screen_for_paralogs",
    "control_not_interpreted_as_phylogeny",
}

REFERENCE_ROLES = {"primary", "secondary"}
EVIDENCE_CONFIDENCE = {"high", "medium", "low", "none", "unknown"}


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

    if label == "ortholog_targets":
        for row in rows:
            gene_id = row.get("gene_id", "").strip()
            category = (row.get("category") or "").strip()
            if category not in TARGET_CATEGORIES:
                raise SystemExit(f"ortholog_targets has invalid category for {gene_id}: {category}")

            orthology_basis = (row.get("orthology_basis") or "").strip()
            if orthology_basis not in ORTHOLOGY_BASIS:
                raise SystemExit(f"ortholog_targets has invalid orthology_basis for {gene_id}: {orthology_basis}")

            copy_number = (row.get("copy_number_expectation") or "").strip()
            if copy_number not in COPY_NUMBER_EXPECTATIONS:
                raise SystemExit(
                    f"ortholog_targets has invalid copy_number_expectation for {gene_id}: {copy_number}"
                )

            try:
                cds_length = int((row.get("cds_length") or "").strip())
            except ValueError as exc:
                raise SystemExit(f"ortholog_targets has non-integer cds_length for {gene_id}") from exc
            if cds_length <= 0:
                raise SystemExit(f"ortholog_targets has non-positive cds_length for {gene_id}")

            if not (row.get("rationale") or "").strip():
                raise SystemExit(f"ortholog_targets is missing a rationale for {gene_id}")

    if label == "reference_manifest":
        primary_count = sum(1 for row in rows if (row.get("reference_role") or "").strip() == "primary")
        if primary_count != 1:
            raise SystemExit(f"reference_manifest must contain exactly one primary reference; found {primary_count}")

        for row in rows:
            role = (row.get("reference_role") or "").strip()
            if role not in REFERENCE_ROLES:
                raise SystemExit(
                    f"reference_manifest has invalid reference_role for {row.get('reference_id', '').strip()}: {role}"
                )

    if label == "species_manifest":
        for row in rows:
            value = (row.get("evidence_confidence") or "").strip()
            if value and value not in EVIDENCE_CONFIDENCE:
                raise SystemExit(
                    f"species_manifest has invalid evidence_confidence for {row.get('species_id', '').strip()}: {value}"
                )

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
