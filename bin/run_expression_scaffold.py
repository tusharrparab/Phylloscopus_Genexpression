#!/usr/bin/env python3

import argparse
import csv
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List


def parse_args():
    parser = argparse.ArgumentParser(
        description="Scaffold cross-species RNA-seq quantification and expression study setup."
    )
    parser.add_argument("--species-manifest", required=True)
    parser.add_argument("--reference-manifest", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--mode", choices=["plan", "scaffold"], default="scaffold")
    parser.add_argument("--threads", type=int, default=2)
    return parser.parse_args()


def read_tsv(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: List[Dict[str, str]], fieldnames: List[str]):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def run_command(cmd: List[str], cwd: Path, log_path: Path) -> int:
    result = subprocess.run(
        cmd,
        cwd=str(cwd),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    log_path.write_text(result.stdout)
    return result.returncode


def split_accessions(value: str) -> List[str]:
    return [item.strip() for item in (value or "").split(",") if item.strip()]


def resolve_local_path(value: str) -> str:
    value = (value or "").strip()
    if value and Path(value).exists():
        return str(Path(value).resolve())
    return ""


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    ref_dir = outdir / "reference"
    plans_dir = outdir / "plans"
    quant_dir = outdir / "quant"
    outdir.mkdir(parents=True, exist_ok=True)
    ref_dir.mkdir(parents=True, exist_ok=True)
    plans_dir.mkdir(parents=True, exist_ok=True)
    quant_dir.mkdir(parents=True, exist_ok=True)

    species_rows = read_tsv(Path(args.species_manifest))
    reference_rows = read_tsv(Path(args.reference_manifest))
    primary_reference = next(
        (row for row in reference_rows if row.get("reference_role", "").strip() == "primary"),
        reference_rows[0] if reference_rows else None,
    )
    if primary_reference is None:
        raise SystemExit("reference_manifest.tsv must contain at least one reference.")

    gffread_exe = shutil.which("gffread")
    salmon_exe = shutil.which("salmon")
    prefetch_exe = shutil.which("prefetch")
    fasterq_exe = shutil.which("fasterq-dump")

    transcript_fasta = resolve_local_path(primary_reference.get("transcript_fasta", ""))
    assembly_fasta = resolve_local_path(primary_reference.get("assembly_fasta", ""))
    annotation_gtf = resolve_local_path(primary_reference.get("annotation_gtf", ""))

    transcript_build_log = ref_dir / "build_transcripts.log"
    if not transcript_fasta and args.mode == "scaffold" and assembly_fasta and annotation_gtf and gffread_exe:
        transcript_fasta = str((ref_dir / "primary_reference_transcripts.fa").resolve())
        code = run_command(
            [gffread_exe, "-w", transcript_fasta, "-g", assembly_fasta, annotation_gtf],
            cwd=ref_dir,
            log_path=transcript_build_log,
        )
        if code != 0:
            transcript_fasta = ""

    salmon_index_dir = ref_dir / "salmon_index"
    salmon_log = ref_dir / "salmon_index.log"
    salmon_index_status = "planned"
    if transcript_fasta and args.mode == "scaffold" and salmon_exe:
        code = run_command(
            [salmon_exe, "index", "-t", transcript_fasta, "-i", str(salmon_index_dir), "-p", str(max(args.threads, 1))],
            cwd=ref_dir,
            log_path=salmon_log,
        )
        salmon_index_status = "completed" if code == 0 else "failed"
    else:
        if not transcript_fasta:
            salmon_index_status = "missing_transcript_fasta"
        elif not salmon_exe:
            salmon_index_status = "salmon_missing"

    sample_rows = []
    plan_rows = []
    design_rows = []

    for species in species_rows:
        scientific_name = species.get("scientific_name", "").strip()
        species_id = species.get("species_id", "").strip()
        run_accessions = split_accessions(species.get("rna_sra_accessions", ""))
        if not run_accessions:
            continue

        for index, accession in enumerate(run_accessions, start=1):
            sample_id = f"{species_id}__{accession}"
            sample_rows.append(
                {
                    "sample_id": sample_id,
                    "species_id": species_id,
                    "scientific_name": scientific_name,
                    "run_accession": accession,
                    "library_layout": "unknown",
                    "library_type": "A",
                    "condition": species_id,
                    "contrast_group": species_id,
                    "batch": "",
                    "notes": "populate experimental condition and batch metadata before DE analysis",
                }
            )
            design_rows.append(
                {
                    "sample_id": sample_id,
                    "condition": species_id,
                    "batch": "",
                    "contrast_group": species_id,
                    "notes": "edit before DESeq2 or edgeR analysis",
                }
            )

            sample_plan_dir = plans_dir / sample_id
            sample_plan_dir.mkdir(parents=True, exist_ok=True)
            reads_dir = outdir / "reads" / accession
            quant_out = quant_dir / sample_id
            plan_script = sample_plan_dir / "run_expression_quant.sh"

            commands = ["#!/usr/bin/env bash", "set -euo pipefail", ""]
            if prefetch_exe and fasterq_exe:
                commands.append(f"{prefetch_exe} {accession}")
                commands.append(f"{fasterq_exe} {accession} -O {reads_dir} --threads {max(args.threads, 1)}")
            else:
                commands.append(f"# Install SRA Toolkit to run prefetch/fasterq-dump for {accession}")
            commands.append("")

            if salmon_exe and transcript_fasta:
                commands.extend(
                    [
                        "# Update this command after confirming whether the run is single-end or paired-end.",
                        (
                            f"# paired-end example:\n# {salmon_exe} quant -i {salmon_index_dir} -l A "
                            f"-1 {reads_dir / (accession + '_1.fastq')} "
                            f"-2 {reads_dir / (accession + '_2.fastq')} "
                            f"--validateMappings -o {quant_out} -p {max(args.threads, 1)}"
                        ),
                        (
                            f"# single-end example:\n# {salmon_exe} quant -i {salmon_index_dir} -l A "
                            f"-r {reads_dir / (accession + '.fastq')} "
                            f"--validateMappings -o {quant_out} -p {max(args.threads, 1)}"
                        ),
                    ]
                )
            else:
                commands.append("# Install salmon and provide a transcript FASTA to run quantification")

            plan_script.write_text("\n".join(commands) + "\n")
            plan_script.chmod(0o755)

            notes = []
            if not prefetch_exe or not fasterq_exe:
                notes.append("sra_tools_missing")
            notes.append("library_layout_unknown")
            if not transcript_fasta:
                notes.append("transcript_fasta_missing")
            if not salmon_exe:
                notes.append("salmon_missing")

            plan_rows.append(
                {
                    "sample_id": sample_id,
                    "species_id": species_id,
                    "run_accession": accession,
                    "plan_script": str(plan_script),
                    "quant_status": "planned" if notes else "ready",
                    "notes": "; ".join(notes) if notes else "ready_to_quantify",
                }
            )

    reference_rows_out = [
        {
            "reference_id": primary_reference.get("reference_id", "").strip(),
            "scientific_name": primary_reference.get("scientific_name", "").strip(),
            "transcript_fasta": transcript_fasta,
            "assembly_fasta": assembly_fasta,
            "annotation_gtf": annotation_gtf,
            "transcript_build_status": (
                "provided"
                if resolve_local_path(primary_reference.get("transcript_fasta", ""))
                else (
                    "built_with_gffread"
                    if transcript_fasta and transcript_build_log.exists()
                    else ("gffread_missing" if not gffread_exe else "planned")
                )
            ),
            "salmon_index_status": salmon_index_status,
        }
    ]

    summary_rows = [
        {"metric": "rna_sample_count", "value": str(len(sample_rows))},
        {"metric": "species_with_rna_runs", "value": str(len({row['species_id'] for row in sample_rows}))},
        {"metric": "salmon_index_status", "value": salmon_index_status},
    ]

    write_tsv(
        outdir / "expression_samples.tsv",
        sample_rows,
        ["sample_id", "species_id", "scientific_name", "run_accession", "library_layout", "library_type", "condition", "contrast_group", "batch", "notes"],
    )
    write_tsv(
        outdir / "expression_design_template.tsv",
        design_rows,
        ["sample_id", "condition", "batch", "contrast_group", "notes"],
    )
    write_tsv(
        outdir / "expression_plan.tsv",
        plan_rows,
        ["sample_id", "species_id", "run_accession", "plan_script", "quant_status", "notes"],
    )
    write_tsv(
        outdir / "expression_reference.tsv",
        reference_rows_out,
        ["reference_id", "scientific_name", "transcript_fasta", "assembly_fasta", "annotation_gtf", "transcript_build_status", "salmon_index_status"],
    )
    write_tsv(outdir / "expression_summary.tsv", summary_rows, ["metric", "value"])


if __name__ == "__main__":
    main()
