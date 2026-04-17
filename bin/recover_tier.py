#!/usr/bin/env python3

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path


SAFE_CODONS = [
    "GCT", "GCC", "GCA", "GCG",
    "CCT", "CCC", "CCA", "CCG",
    "ACT", "ACC", "ACA", "ACG",
    "GTT", "GTC", "GTA", "GTG",
    "GAA", "GAG", "CAA", "CAG",
    "AAT", "AAC", "TAT", "TAC",
    "CAT", "CAC", "AAA", "AAG",
    "CGT", "CGC", "CGA", "CGG",
    "AGA", "AGG", "TCT", "TCC",
    "TCA", "TCG", "AGT", "AGC",
]

CONFIDENCE_BY_TIER = {
    "A": "high",
    "B": "medium",
    "C": "medium",
    "D": "low",
}

METHOD_BY_TIER = {
    "A": "projected_from_annotated_assembly",
    "B": "reference_guided_projection",
    "C": "transcript_or_rnaseq_recovery",
    "D": "targeted_consensus_recovery",
}


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tier-label", required=True)
    parser.add_argument("--accepted-tiers", required=True)
    parser.add_argument("--species-plan", required=True)
    parser.add_argument("--ortholog-targets", required=True)
    parser.add_argument("--reference-manifest", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--mode", choices=["stub", "contract"], default="stub")
    return parser.parse_args()


def read_tsv(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows, fieldnames):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def chunk_string(text, width):
    return [text[index:index + width] for index in range(0, len(text), width)]


def make_mock_sequence(species_id: str, gene_id: str, cds_length: int):
    codon_count = max(30, cds_length // 3)
    digest = hashlib.sha256(f"{species_id}:{gene_id}".encode("utf-8")).digest()
    codons = ["ATG"]
    for index in range(max(codon_count - 1, 1)):
        codons.append(SAFE_CODONS[digest[index % len(digest)] % len(SAFE_CODONS)])
    return "".join(codons[:codon_count])


def main():
    args = parse_args()
    accepted_tiers = {value.strip() for value in args.accepted_tiers.split(",") if value.strip()}
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    sequence_root = outdir / "ortholog_sequences"
    sequence_root.mkdir(parents=True, exist_ok=True)

    species_rows = read_tsv(Path(args.species_plan))
    target_rows = read_tsv(Path(args.ortholog_targets))

    status_rows = []
    tier_summary = []
    tier_counts = Counter()

    for species in species_rows:
        tier = species["evidence_tier"].strip()
        if tier not in accepted_tiers:
            continue

        species_id = species["species_id"].strip()
        scientific_name = species["scientific_name"].strip()
        tier_counts[tier] += 1

        species_dir = sequence_root / species_id
        species_dir.mkdir(parents=True, exist_ok=True)
        species_fasta_path = species_dir / f"{species_id}.cds.fasta"
        fasta_lines = []

        for target in target_rows:
            gene_id = target["gene_id"].strip()
            cds_length = int(target["cds_length"])

            if args.mode == "stub":
                sequence = make_mock_sequence(species_id, gene_id, cds_length)
                header = f">{species_id}|{gene_id}|tier={tier}|mode=stub"
                fasta_lines.extend([header, *chunk_string(sequence, 60)])
                reconstruction_status = "reconstructed"
                sequence_length = len(sequence)
                notes = "mock sequence emitted for scaffold validation"
            else:
                reconstruction_status = "planned"
                sequence_length = 0
                notes = "contract mode only; replace module command with real tool chain"

            status_rows.append(
                {
                    "species_id": species_id,
                    "scientific_name": scientific_name,
                    "gene_id": gene_id,
                    "gene_symbol": target["gene_symbol"].strip(),
                    "evidence_tier": tier,
                    "reconstruction_status": reconstruction_status,
                    "sequence_length": str(sequence_length),
                    "confidence_tier": CONFIDENCE_BY_TIER[tier],
                    "method": METHOD_BY_TIER[tier],
                    "notes": notes,
                }
            )

        if fasta_lines:
            species_fasta_path.write_text("\n".join(fasta_lines) + "\n")

        tier_summary.append(
            {
                "species_id": species_id,
                "scientific_name": scientific_name,
                "evidence_tier": tier,
                "genes_requested": str(len(target_rows)),
                "mode": args.mode,
            }
        )

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
    write_tsv(
        outdir / "tier_summary.tsv",
        tier_summary,
        ["species_id", "scientific_name", "evidence_tier", "genes_requested", "mode"],
    )

    metadata = {
        "tier_label": args.tier_label,
        "accepted_tiers": sorted(accepted_tiers),
        "mode": args.mode,
        "species_count": sum(tier_counts.values()),
        "tier_counts": dict(tier_counts),
    }
    (outdir / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")


if __name__ == "__main__":
    main()

