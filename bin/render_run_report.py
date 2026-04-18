#!/usr/bin/env python3

import argparse
import csv
import json
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--status-matrix", required=True)
    parser.add_argument("--status-long", required=True)
    parser.add_argument("--sequence-dir", required=True)
    parser.add_argument("--recovery-summary", required=True)
    parser.add_argument("--run-summary", required=True)
    parser.add_argument("--outdir", required=True)
    return parser.parse_args()


def read_tsv(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    status_matrix_rows = read_tsv(Path(args.status_matrix))
    status_long_rows = read_tsv(Path(args.status_long))
    recovery_summary_rows = read_tsv(Path(args.recovery_summary))
    run_summary = json.loads(Path(args.run_summary).read_text())

    summary_lookup = {row["metric"]: row["value"] for row in recovery_summary_rows}
    category_counts = {}
    for key, value in summary_lookup.items():
        if key.startswith("category_"):
            category_counts[key.removeprefix("category_")] = value

    lines = [
        "# Phylloscopus Comparative Ortholog Recovery Report",
        "",
        "## Run Summary",
        "",
        f"- Species in manifest: {run_summary['species_count']}",
        f"- Ortholog targets: {run_summary['ortholog_target_count']}",
        f"- References: {run_summary['reference_count']}",
        "",
        "## Tier Distribution",
        "",
    ]

    for tier in ["A", "B", "C", "D", "E"]:
        lines.append(f"- Tier {tier}: {run_summary['tiers'].get(tier, 0)} species")

    lines.extend(
        [
            "",
            "## Target Panel",
            "",
            f"- Phylogenetic backbone targets: {category_counts.get('phylogenetic_backbone', '0')}",
            f"- Migration-linked candidate targets: {category_counts.get('migration_candidate', '0')}",
            f"- Vocalization/neural candidate targets: {category_counts.get('vocalization_neural_candidate', '0')}",
            f"- Hypoxia/elevation candidate targets: {category_counts.get('hypoxia_elevation_candidate', '0')}",
            f"- Housekeeping controls: {category_counts.get('housekeeping_control', '0')}",
            f"- Single-copy-preferred targets: {summary_lookup.get('targets_single_copy_preferred', '0')}",
            f"- Targets requiring explicit paralog screening: {summary_lookup.get('targets_paralogy_screen_required', '0')}",
            "",
            "## Recovery Totals",
            "",
            f"- Tier A recovered: {sum(1 for row in status_long_rows if row.get('evidence_tier') == 'A' and row.get('reconstruction_status') == 'recovered')}",
            f"- Tier B candidate_sequence_recovered: {sum(1 for row in status_long_rows if row.get('evidence_tier') == 'B' and row.get('reconstruction_status') == 'candidate_sequence_recovered')}",
            f"- Recovered sequences: {summary_lookup.get('status_recovered', '0')}",
            f"- Real reconstructed sequences: {summary_lookup.get('status_reconstructed', '0')}",
            f"- Candidate assembly-backed sequences: {summary_lookup.get('status_candidate_sequence_recovered', '0')}",
            f"- Stub sequences emitted: {summary_lookup.get('status_stub_sequence_emitted', '0')}",
            f"- Planned only: {summary_lookup.get('status_planned', '0')}",
            f"- Missing data: {summary_lookup.get('status_missing_data', '0')}",
            "",
            "## Matrix Preview",
            "",
            "| species_id | scientific_name | evidence_tier |",
            "|---|---|---|",
        ]
    )

    for row in status_matrix_rows[:10]:
        lines.append(
            f"| {row['species_id']} | {row['scientific_name']} | {row['evidence_tier']} |"
        )

    lines.extend(
        [
            "",
            "## Notes",
            "",
            "- This repository is a comparative ortholog-recovery scaffold, not a completed biological inference pipeline.",
            "- `stub` mode proves workflow wiring only. Synthetic FASTA emitted by stub modules is not biological sequence.",
            "- Tier B candidate sequences are assembly-backed locus recoveries only and are not validated ortholog calls.",
            "- Tier C and Tier D outputs remain hypothesis-level placeholders until orthology is validated from real RNA or WGS data.",
            "- Candidate-gene panel categories support hypothesis generation only; they do not by themselves establish causal trait associations.",
            "- Single-copy-preferred targets still require locus-level inspection in the chosen reference and query assemblies.",
            f"- Consolidated ortholog FASTA bundles are written under `{args.sequence_dir}`.",
            "",
        ]
    )

    if status_long_rows:
        unresolved = sum(
            1 for row in status_long_rows if row.get("reconstruction_status") in {"planned", "missing_data"}
        )
        lines.extend(
            [
                "## Interpretation Guardrails",
                "",
                f"- Non-final target rows (planned or missing): {unresolved}",
                "- Treat downstream ASR as defensible only for loci reconstructed from real homologous sequence, not from stub outputs.",
                "- Treat expression outputs as optional RNA-backed side analyses rather than the defining product of this repository.",
                "",
            ]
        )

    (outdir / "pipeline_report.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
