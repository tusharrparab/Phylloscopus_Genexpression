#!/usr/bin/env python3

import argparse
import csv
import re
from pathlib import Path
from typing import Dict, List, Optional


ASSEMBLY_LEVEL_RANK = {
    "Complete Genome": 0,
    "Chromosome": 1,
    "Scaffold": 2,
    "Contig": 3,
    "": 99,
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build a reference_manifest.tsv from a discovery species snapshot."
    )
    parser.add_argument("--species-manifest", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument(
        "--max-references",
        type=int,
        default=6,
        help="Maximum number of references to emit, including the primary reference.",
    )
    parser.add_argument(
        "--primary-species",
        default="",
        help="Optional scientific name or species_id to force as primary reference.",
    )
    parser.add_argument(
        "--busco-lineage",
        default="",
        help="Optional BUSCO lineage dataset to record in the manifest. Leave blank to use auto-lineage.",
    )
    return parser.parse_args()


def read_tsv(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: List[Dict[str, str]], fieldnames: List[str]):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def slugify(value: str) -> str:
    value = value.strip().lower()
    value = re.sub(r"[^a-z0-9]+", "_", value)
    return value.strip("_")


def sort_key(row: Dict[str, str]):
    has_local_annotation = bool((row.get("annotation_gtf") or "").strip())
    has_local_transcriptome = bool(
        (row.get("transcriptome_fasta") or row.get("transcript_fasta") or "").strip()
    )
    has_local_assembly = bool((row.get("assembly_fasta") or "").strip())
    level = (row.get("assembly_level") or "").strip()
    level_rank = ASSEMBLY_LEVEL_RANK.get(level, 99)
    assembly_count = int((row.get("ncbi_assembly_count") or "0").strip() or 0)
    return (
        0 if has_local_annotation and (has_local_assembly or has_local_transcriptome) else 1,
        0 if has_local_transcriptome else 1,
        0 if has_local_assembly else 1,
        level_rank,
        -assembly_count,
        row.get("scientific_name", "").strip(),
    )


def choose_primary(candidates: List[Dict[str, str]], preferred: str) -> Optional[Dict[str, str]]:
    if not candidates:
        return None

    if preferred:
        preferred_norm = preferred.strip().lower()
        for row in candidates:
            if row.get("scientific_name", "").strip().lower() == preferred_norm:
                return row
            if row.get("species_id", "").strip().lower() == preferred_norm:
                return row

    ordered = sorted(candidates, key=sort_key)
    return ordered[0]


def main():
    args = parse_args()
    species_rows = read_tsv(Path(args.species_manifest))

    candidates = [
        row
        for row in species_rows
        if (row.get("assembly_accession") or "").strip() or (row.get("assembly_fasta") or "").strip()
    ]
    if not candidates:
        raise SystemExit("No assembly-backed species were found in the species manifest.")

    primary = choose_primary(candidates, args.primary_species)
    if primary is None:
        raise SystemExit("Unable to determine a primary reference.")

    ordered = [row for row in sorted(candidates, key=sort_key) if row != primary]
    selected = [primary, *ordered[: max(args.max_references - 1, 0)]]

    rows = []
    for index, row in enumerate(selected):
        scientific_name = row.get("scientific_name", "").strip()
        reference_role = "primary" if index == 0 else "secondary"
        species_slug = slugify(scientific_name)
        reference_id = f"{species_slug}_{reference_role}"

        notes = (row.get("notes") or "").strip()
        source_note = f"derived from snapshot species manifest; assembly_accession={row.get('assembly_accession', '').strip() or 'local-only'}"
        notes = f"{source_note}; {notes}" if notes else source_note

        rows.append(
            {
                "reference_id": reference_id,
                "scientific_name": scientific_name,
                "reference_role": reference_role,
                "reference_quality": (
                    "annotation_ready"
                    if (row.get("assembly_fasta") or "").strip() and (row.get("annotation_gtf") or "").strip()
                    else "assembly_only"
                ),
                "data_provenance": (row.get("data_provenance") or "").strip() or "snapshot_manifest",
                "assembly_accession": (row.get("assembly_accession") or "").strip(),
                "assembly_level": (row.get("assembly_level") or "").strip(),
                "assembly_name": (row.get("assembly_name") or "").strip(),
                "assembly_fasta": (row.get("assembly_fasta") or "").strip(),
                "annotation_gtf": (row.get("annotation_gtf") or "").strip(),
                "protein_fasta": "",
                "transcript_fasta": (row.get("transcriptome_fasta") or row.get("transcript_fasta") or "").strip(),
                "reference_bed12": "",
                "reference_twobit": "",
                "query_chain": "",
                "busco_lineage": args.busco_lineage,
                "analysis_notes": (row.get("analysis_suitability") or "").strip(),
                "notes": notes,
            }
        )

    write_tsv(
        Path(args.out),
        rows,
        [
            "reference_id",
            "scientific_name",
            "reference_role",
            "reference_quality",
            "data_provenance",
            "assembly_accession",
            "assembly_level",
            "assembly_name",
            "assembly_fasta",
            "annotation_gtf",
            "protein_fasta",
            "transcript_fasta",
            "reference_bed12",
            "reference_twobit",
            "query_chain",
            "busco_lineage",
            "analysis_notes",
            "notes",
        ],
    )


if __name__ == "__main__":
    main()
