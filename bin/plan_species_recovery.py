#!/usr/bin/env python3

import argparse
import csv
import json
from collections import Counter
from pathlib import Path


TIER_TO_STRATEGY = {
    "A": "assembly_annotation_projection",
    "B": "assembly_reference_guided_projection",
    "C": "transcript_or_rnaseq_recovery",
    "D": "targeted_consensus_recovery",
    "E": "missing_data_report",
}

DEFAULT_CONFIDENCE_BY_TIER = {
    "A": "high",
    "B": "medium",
    "C": "low",
    "D": "low",
    "E": "none",
}

DEFAULT_SUITABILITY_BY_TIER = {
    "A": "assembly_anchored_comparative_recovery",
    "B": "projection_candidate_pending_annotation",
    "C": "rna_backed_candidate_not_orthology_validated",
    "D": "wgs_backed_candidate_not_orthology_validated",
    "E": "not_analyzable_with_current_inputs",
}


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--species-manifest", required=True)
    parser.add_argument("--ortholog-targets", required=True)
    parser.add_argument("--reference-manifest", required=True)
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


def classify_tier(row):
    hint = (row.get("evidence_hint") or "").strip().upper()
    if hint in {"A", "B", "C", "D", "E"}:
        return hint

    assembly_fasta = (row.get("assembly_fasta") or "").strip()
    annotation_gtf = (row.get("annotation_gtf") or "").strip()
    assembly_accession = (row.get("assembly_accession") or "").strip()
    transcriptome_fasta = (row.get("transcriptome_fasta") or "").strip()
    rna_sra = (row.get("rna_sra_accessions") or "").strip()
    wgs_sra = (row.get("wgs_sra_accessions") or "").strip()

    if assembly_fasta and annotation_gtf:
        return "A"
    if assembly_fasta or assembly_accession:
        return "B"
    if transcriptome_fasta or rna_sra:
        return "C"
    if wgs_sra:
        return "D"
    return "E"


def default_provenance(row):
    tags = []
    if (row.get("gbif_species_key") or "").strip():
        tags.append("gbif")
    if (row.get("assembly_accession") or "").strip() or (row.get("ncbi_assembly_count") or "").strip():
        tags.append("ncbi")
    if (row.get("ena_tax_id") or "").strip() or (row.get("ena_read_run_count") or "").strip():
        tags.append("ena")
    if (row.get("assembly_fasta") or "").strip() or (row.get("annotation_gtf") or "").strip():
        tags.append("local")
    return ";".join(tags) if tags else "unspecified"


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    species_rows = read_tsv(Path(args.species_manifest))
    target_rows = read_tsv(Path(args.ortholog_targets))
    reference_rows = read_tsv(Path(args.reference_manifest))

    target_count = len(target_rows)
    tier_counts = Counter()
    planned_rows = []

    for row in species_rows:
        tier = classify_tier(row)
        tier_counts[tier] += 1
        planned_rows.append(
            {
                "species_id": row["species_id"].strip(),
                "scientific_name": row["scientific_name"].strip(),
                "evidence_tier": tier,
                "evidence_confidence": (row.get("evidence_confidence") or "").strip() or DEFAULT_CONFIDENCE_BY_TIER[tier],
                "data_provenance": (row.get("data_provenance") or "").strip() or default_provenance(row),
                "analysis_suitability": (row.get("analysis_suitability") or "").strip() or DEFAULT_SUITABILITY_BY_TIER[tier],
                "recovery_strategy": TIER_TO_STRATEGY[tier],
                "assembly_source": (row.get("assembly_fasta") or row.get("assembly_accession") or "").strip(),
                "annotation_source": (row.get("annotation_gtf") or "").strip(),
                "transcriptome_source": (row.get("transcriptome_fasta") or row.get("rna_sra_accessions") or "").strip(),
                "wgs_source": (row.get("wgs_sra_accessions") or "").strip(),
                "ortholog_target_count": str(target_count),
                "notes": (row.get("notes") or "").strip(),
            }
        )

    write_tsv(
        outdir / "species_plan.tsv",
        planned_rows,
        [
            "species_id",
            "scientific_name",
            "evidence_tier",
            "evidence_confidence",
            "data_provenance",
            "analysis_suitability",
            "recovery_strategy",
            "assembly_source",
            "annotation_source",
            "transcriptome_source",
            "wgs_source",
            "ortholog_target_count",
            "notes",
        ],
    )

    run_summary = {
        "species_count": len(species_rows),
        "ortholog_target_count": target_count,
        "reference_count": len(reference_rows),
        "tiers": {tier: tier_counts.get(tier, 0) for tier in ["A", "B", "C", "D", "E"]},
    }
    (outdir / "run_summary.json").write_text(json.dumps(run_summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
