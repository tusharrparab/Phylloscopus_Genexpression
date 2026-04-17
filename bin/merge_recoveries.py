#!/usr/bin/env python3

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--species-plan", required=True)
    parser.add_argument("--ortholog-targets", required=True)
    parser.add_argument("--tier-dir", action="append", default=[])
    parser.add_argument("--tier-status", action="append", default=[])
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


def read_fasta(path: Path):
    header = None
    chunks = []
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:]
                chunks = []
            else:
                chunks.append(line)
    if header is not None:
        yield header, "".join(chunks)


def chunk_string(text, width):
    return [text[index:index + width] for index in range(0, len(text), width)]


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    merged_sequences_dir = outdir / "ortholog_sequences"
    merged_sequences_dir.mkdir(parents=True, exist_ok=True)

    species_rows = read_tsv(Path(args.species_plan))
    species_lookup = {row["species_id"].strip(): row for row in species_rows}
    target_rows = read_tsv(Path(args.ortholog_targets))
    gene_order = [row["gene_id"].strip() for row in target_rows]

    all_status_rows = []
    for status_path in args.tier_status:
        all_status_rows.extend(read_tsv(Path(status_path)))

    all_status_rows.sort(key=lambda row: (row["species_id"], row["gene_id"]))

    write_tsv(
        outdir / "ortholog_status_long.tsv",
        all_status_rows,
        [
            "species_id",
            "scientific_name",
            "gene_id",
            "gene_symbol",
            "target_category",
            "orthology_basis",
            "copy_number_expectation",
            "target_rationale",
            "evidence_tier",
            "evidence_confidence",
            "reconstruction_status",
            "sequence_length",
            "confidence_tier",
            "method",
            "notes",
        ],
    )

    matrix = {}
    for species_id, species in species_lookup.items():
        matrix[species_id] = {
            "species_id": species_id,
            "scientific_name": species["scientific_name"].strip(),
            "evidence_tier": species["evidence_tier"].strip(),
        }
        for gene_id in gene_order:
            matrix[species_id][gene_id] = "not_evaluated"

    for row in all_status_rows:
        matrix[row["species_id"]][row["gene_id"]] = row["reconstruction_status"]

    matrix_rows = [matrix[species_id] for species_id in sorted(matrix)]
    write_tsv(
        outdir / "ortholog_status_matrix.tsv",
        matrix_rows,
        ["species_id", "scientific_name", "evidence_tier", *gene_order],
    )

    gene_to_entries = defaultdict(list)
    for tier_dir in args.tier_dir:
        fasta_dir = Path(tier_dir) / "ortholog_sequences"
        if not fasta_dir.exists():
            continue
        for fasta_path in fasta_dir.glob("*/*.cds.fasta"):
            for header, sequence in read_fasta(fasta_path):
                parts = header.split("|")
                if len(parts) < 2:
                    continue
                species_id = parts[0]
                gene_id = parts[1]
                gene_to_entries[gene_id].append((species_id, sequence))

    for gene_id in gene_order:
        entries = sorted(gene_to_entries.get(gene_id, []))
        gene_fasta_path = merged_sequences_dir / f"{gene_id}.fna"
        with gene_fasta_path.open("w") as handle:
            for species_id, sequence in entries:
                handle.write(f">{species_id}|{gene_id}\n")
                handle.write("\n".join(chunk_string(sequence, 60)))
                handle.write("\n")

    status_counts = Counter(row["reconstruction_status"] for row in all_status_rows)
    category_counts = Counter(row["category"].strip() for row in target_rows)
    strict_orthology_targets = sum(
        1
        for row in target_rows
        if (row.get("copy_number_expectation") or "").strip() == "single_copy_preferred"
    )
    summary_rows = [
        {"metric": "species_count", "value": str(len(species_rows))},
        {"metric": "gene_count", "value": str(len(gene_order))},
        {"metric": "targets_single_copy_preferred", "value": str(strict_orthology_targets)},
        {
            "metric": "targets_paralogy_screen_required",
            "value": str(
                sum(
                    1
                    for row in target_rows
                    if (row.get("copy_number_expectation") or "").strip() == "screen_for_paralogs"
                )
            ),
        },
        {"metric": "category_phylogenetic_backbone", "value": str(category_counts.get("phylogenetic_backbone", 0))},
        {"metric": "category_migration_candidate", "value": str(category_counts.get("migration_candidate", 0))},
        {
            "metric": "category_vocalization_neural_candidate",
            "value": str(category_counts.get("vocalization_neural_candidate", 0)),
        },
        {
            "metric": "category_hypoxia_elevation_candidate",
            "value": str(category_counts.get("hypoxia_elevation_candidate", 0)),
        },
        {"metric": "category_housekeeping_control", "value": str(category_counts.get("housekeeping_control", 0))},
        {"metric": "status_reconstructed", "value": str(status_counts.get("reconstructed", 0))},
        {"metric": "status_stub_sequence_emitted", "value": str(status_counts.get("stub_sequence_emitted", 0))},
        {"metric": "status_planned", "value": str(status_counts.get("planned", 0))},
        {"metric": "status_missing_data", "value": str(status_counts.get("missing_data", 0))},
    ]
    write_tsv(outdir / "recovery_summary.tsv", summary_rows, ["metric", "value"])


if __name__ == "__main__":
    main()
