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
    recovery_summary_rows = read_tsv(Path(args.recovery_summary))
    run_summary = json.loads(Path(args.run_summary).read_text())

    summary_lookup = {row["metric"]: row["value"] for row in recovery_summary_rows}

    lines = [
        "# Phylloscopus Ortholog Pipeline Report",
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
            "## Recovery Totals",
            "",
            f"- Reconstructed: {summary_lookup.get('status_reconstructed', '0')}",
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
            "- This scaffold is designed to run end-to-end in `stub` mode.",
            "- Replace tier recovery commands with real tool invocations before treating any emitted FASTA as biological sequence.",
            f"- Consolidated ortholog FASTA bundles are written under `{args.sequence_dir}`.",
            "",
        ]
    )

    (outdir / "pipeline_report.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()

