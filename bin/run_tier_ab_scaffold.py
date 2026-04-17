#!/usr/bin/env python3

import argparse
import csv
import json
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from zipfile import ZipFile


def parse_args():
    parser = argparse.ArgumentParser(
        description="Stage assembly-backed Tier A/B species for download, BUSCO, and projection scaffolding."
    )
    parser.add_argument("--species-plan", required=True)
    parser.add_argument("--ortholog-targets", required=True)
    parser.add_argument("--reference-manifest", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument(
        "--mode",
        choices=["plan", "scaffold"],
        default="scaffold",
        help="plan emits command plans only; scaffold executes steps when the required tools are present.",
    )
    parser.add_argument("--cpus", type=int, default=2)
    parser.add_argument("--ncbi-api-key", default="")
    parser.add_argument("--busco-lineage", default="")
    return parser.parse_args()


def read_tsv(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: List[Dict[str, str]], fieldnames: List[str]):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def run_command(cmd: List[str], cwd: Optional[Path], log_path: Path) -> Tuple[int, str]:
    result = subprocess.run(
        cmd,
        cwd=str(cwd) if cwd else None,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    log_path.write_text(result.stdout)
    return result.returncode, result.stdout


def find_first(root: Path, patterns: List[str]) -> str:
    for pattern in patterns:
        matches = sorted(root.rglob(pattern))
        if matches:
            return str(matches[0].resolve())
    return ""


def resolve_local_asset(row: Dict[str, str], key: str) -> str:
    value = (row.get(key) or "").strip()
    if value and Path(value).exists():
        return str(Path(value).resolve())
    return ""


def stage_ncbi_datasets(accession: str, stage_dir: Path, api_key: str) -> Dict[str, str]:
    zip_path = stage_dir / "ncbi_dataset.zip"
    unpack_dir = stage_dir / "ncbi_dataset"
    log_path = stage_dir / "datasets_download.log"
    stage_dir.mkdir(parents=True, exist_ok=True)

    if not shutil.which("datasets"):
        return {
            "status": "datasets_missing",
            "genome_fasta": "",
            "annotation_gtf": "",
            "protein_fasta": "",
            "gff3": "",
            "zip_path": str(zip_path),
            "log_path": str(log_path),
        }

    if not unpack_dir.exists():
        cmd = [
            "datasets",
            "download",
            "genome",
            "accession",
            accession,
            "--include",
            "genome,gff3,gtf,protein,cds",
            "--filename",
            str(zip_path.resolve()),
            "--no-progressbar",
        ]
        if api_key:
            cmd.extend(["--api-key", api_key])
        code = 1
        for attempt in range(1, 4):
            code, _ = run_command(cmd, cwd=stage_dir, log_path=log_path)
            if code == 0:
                break
            zip_path.unlink(missing_ok=True)
            log_path.write_text(log_path.read_text() + f"\n[retry] datasets attempt {attempt} failed\n")
        if code != 0:
            return {
                "status": "datasets_failed",
                "genome_fasta": "",
                "annotation_gtf": "",
                "protein_fasta": "",
                "gff3": "",
                "zip_path": str(zip_path),
                "log_path": str(log_path),
            }
        with ZipFile(zip_path) as archive:
            archive.extractall(unpack_dir)

    return {
        "status": "downloaded",
        "genome_fasta": find_first(unpack_dir, ["*_genomic.fna", "*.fna"]),
        "annotation_gtf": find_first(unpack_dir, ["*.gtf"]),
        "protein_fasta": find_first(unpack_dir, ["*protein.faa", "*.faa"]),
        "gff3": find_first(unpack_dir, ["*.gff3", "*.gff"]),
        "zip_path": str(zip_path),
        "log_path": str(log_path),
    }


def run_busco(genome_fasta: str, species_id: str, stage_dir: Path, cpus: int, busco_lineage: str) -> Dict[str, str]:
    output_root = stage_dir / "busco"
    log_path = stage_dir / "busco.log"

    if not genome_fasta:
        return {"status": "skipped_no_genome", "summary_json": "", "selected_lineage": "", "log_path": str(log_path)}
    if not shutil.which("busco"):
        return {"status": "busco_missing", "summary_json": "", "selected_lineage": "", "log_path": str(log_path)}

    output_root.mkdir(parents=True, exist_ok=True)
    cmd = [
        "busco",
        "-i",
        genome_fasta,
        "-m",
        "genome",
        "-o",
        species_id,
        "--out_path",
        str(output_root),
        "-c",
        str(max(cpus, 1)),
    ]
    if busco_lineage:
        cmd.extend(["-l", busco_lineage])
    else:
        cmd.append("--auto-lineage-euk")

    code, _ = run_command(cmd, cwd=stage_dir, log_path=log_path)
    summary_json = find_first(output_root, ["short_summary*.json"])
    selected_lineage = ""
    if summary_json:
        try:
            data = json.loads(Path(summary_json).read_text())
            lineage_dataset = data.get("lineage_dataset") or data.get("lineage")
            if isinstance(lineage_dataset, dict):
                selected_lineage = lineage_dataset.get("name", "")
            elif isinstance(lineage_dataset, str):
                selected_lineage = lineage_dataset
        except Exception:
            selected_lineage = ""

    return {
        "status": "completed" if code == 0 else "busco_failed",
        "summary_json": summary_json,
        "selected_lineage": selected_lineage,
        "log_path": str(log_path),
    }


def write_projection_plan(
    projection_dir: Path,
    species_id: str,
    query_genome_fasta: str,
    reference_row: Dict[str, str],
    reference_assets: Dict[str, str],
) -> Dict[str, str]:
    projection_dir.mkdir(parents=True, exist_ok=True)

    reference_bed12 = resolve_local_asset(reference_row, "reference_bed12")
    reference_twobit = resolve_local_asset(reference_row, "reference_twobit")
    query_chain = (reference_row.get("query_chain") or "").strip()
    query_twobit = ""
    missing = []

    if not reference_bed12:
        missing.append("reference_bed12")
    if not reference_twobit:
        missing.append("reference_twobit")
    if not query_chain:
        missing.append("query_chain")
    if not query_twobit:
        missing.append("query_twobit")

    command_lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        "",
        "# Fill missing prerequisites before executing TOGA.",
        f"# query species: {species_id}",
        f"# primary reference: {reference_row.get('scientific_name', '').strip()}",
        "",
    ]

    toga_exe = shutil.which("toga.py") or shutil.which("TOGA") or "toga.py"
    command_lines.append(
        f"# {toga_exe} {query_chain or '<CHAIN_FILE>'} {reference_bed12 or '<REFERENCE_BED12>'} "
        f"{reference_twobit or '<REFERENCE_2BIT>'} {query_twobit or '<QUERY_2BIT>'} "
        f"--project_dir {projection_dir / 'toga_project'}"
    )

    plan_path = projection_dir / "run_toga.sh"
    plan_path.write_text("\n".join(command_lines) + "\n")
    plan_path.chmod(0o755)

    metadata = {
        "query_species_id": species_id,
        "query_genome_fasta": query_genome_fasta,
        "reference_id": reference_row.get("reference_id", "").strip(),
        "reference_scientific_name": reference_row.get("scientific_name", "").strip(),
        "reference_annotation_gtf": reference_assets.get("annotation_gtf", ""),
        "reference_protein_fasta": reference_assets.get("protein_fasta", ""),
        "reference_bed12": reference_bed12,
        "reference_twobit": reference_twobit,
        "query_chain": query_chain,
        "query_twobit": query_twobit,
        "missing_prerequisites": missing,
        "toga_executable_found": bool(shutil.which("toga.py") or shutil.which("TOGA")),
    }
    metadata_path = projection_dir / "projection_plan.json"
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")

    return {
        "status": "ready" if not missing else "planned_missing_prerequisites",
        "plan_script": str(plan_path),
        "plan_metadata": str(metadata_path),
        "missing_prerequisites": ",".join(missing),
    }


def stage_asset_bundle(
    row: Dict[str, str],
    stage_dir: Path,
    api_key: str,
    mode: str,
) -> Dict[str, str]:
    assembly_fasta = resolve_local_asset(row, "assembly_fasta")
    annotation_gtf = resolve_local_asset(row, "annotation_gtf")
    protein_fasta = resolve_local_asset(row, "protein_fasta")
    accession = (row.get("assembly_accession") or "").strip()

    if assembly_fasta:
        return {
            "status": "local",
            "genome_fasta": assembly_fasta,
            "annotation_gtf": annotation_gtf,
            "protein_fasta": protein_fasta,
            "gff3": "",
            "log_path": "",
            "zip_path": "",
        }

    if mode == "scaffold" and accession:
        return stage_ncbi_datasets(accession, stage_dir, api_key)

    return {
        "status": "planned_download" if accession else "missing_source",
        "genome_fasta": "",
        "annotation_gtf": "",
        "protein_fasta": "",
        "gff3": "",
        "log_path": "",
        "zip_path": "",
    }


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    assets_dir = outdir / "assets"
    ortholog_sequences_dir = outdir / "ortholog_sequences"
    outdir.mkdir(parents=True, exist_ok=True)
    ortholog_sequences_dir.mkdir(parents=True, exist_ok=True)

    species_rows = read_tsv(Path(args.species_plan))
    target_rows = read_tsv(Path(args.ortholog_targets))
    reference_rows = read_tsv(Path(args.reference_manifest))

    query_species = [row for row in species_rows if row.get("evidence_tier", "").strip() in {"A", "B"}]
    primary_reference = next(
        (row for row in reference_rows if row.get("reference_role", "").strip() == "primary"),
        reference_rows[0] if reference_rows else None,
    )
    if primary_reference is None:
        raise SystemExit("reference_manifest.tsv must contain at least one reference row.")

    reference_assets_rows = []
    status_rows = []
    projection_rows = []

    primary_reference_stage = assets_dir / "references" / primary_reference["reference_id"].strip()
    primary_reference_assets = stage_asset_bundle(primary_reference, primary_reference_stage, args.ncbi_api_key, args.mode)
    reference_assets_rows.append(
        {
            "entity_type": "reference",
            "entity_id": primary_reference["reference_id"].strip(),
            "scientific_name": primary_reference["scientific_name"].strip(),
            "assembly_accession": (primary_reference.get("assembly_accession") or "").strip(),
            "asset_status": primary_reference_assets["status"],
            "genome_fasta": primary_reference_assets["genome_fasta"],
            "annotation_gtf": primary_reference_assets["annotation_gtf"],
            "protein_fasta": primary_reference_assets["protein_fasta"],
            "log_path": primary_reference_assets["log_path"],
        }
    )

    for species in query_species:
        species_id = species["species_id"].strip()
        scientific_name = species["scientific_name"].strip()
        stage_dir = assets_dir / "queries" / species_id
        assembly_source = (species.get("assembly_source") or "").strip()
        assembly_source_is_file = bool(assembly_source) and Path(assembly_source).exists()
        asset_row = {
            "assembly_accession": assembly_source if assembly_source and not assembly_source_is_file else "",
            "assembly_fasta": assembly_source if assembly_source_is_file else "",
            "annotation_gtf": "",
            "protein_fasta": "",
        }
        assets = stage_asset_bundle(asset_row, stage_dir, args.ncbi_api_key, args.mode)
        busco = run_busco(
            assets["genome_fasta"],
            species_id,
            stage_dir,
            args.cpus,
            args.busco_lineage or (primary_reference.get("busco_lineage") or "").strip(),
        ) if args.mode == "scaffold" else {
            "status": "planned",
            "summary_json": "",
            "selected_lineage": "",
            "log_path": str(stage_dir / "busco.log"),
        }
        projection = write_projection_plan(
            stage_dir / "projection",
            species_id,
            assets["genome_fasta"],
            primary_reference,
            primary_reference_assets,
        )

        reference_assets_rows.append(
            {
                "entity_type": "query",
                "entity_id": species_id,
                "scientific_name": scientific_name,
                "assembly_accession": asset_row["assembly_accession"],
                "asset_status": assets["status"],
                "genome_fasta": assets["genome_fasta"],
                "annotation_gtf": assets["annotation_gtf"],
                "protein_fasta": assets["protein_fasta"],
                "log_path": assets["log_path"],
            }
        )

        projection_rows.append(
            {
                "species_id": species_id,
                "scientific_name": scientific_name,
                "reference_id": primary_reference["reference_id"].strip(),
                "projection_status": projection["status"],
                "plan_script": projection["plan_script"],
                "plan_metadata": projection["plan_metadata"],
                "missing_prerequisites": projection["missing_prerequisites"],
            }
        )

        note_parts = [
            f"asset_status={assets['status']}",
            f"busco_status={busco['status']}",
            f"projection_status={projection['status']}",
        ]
        if busco["selected_lineage"]:
            note_parts.append(f"busco_lineage={busco['selected_lineage']}")
        if projection["missing_prerequisites"]:
            note_parts.append(f"missing={projection['missing_prerequisites']}")
        notes = "; ".join(note_parts)

        for target in target_rows:
            status_rows.append(
                {
                    "species_id": species_id,
                    "scientific_name": scientific_name,
                    "gene_id": target["gene_id"].strip(),
                    "gene_symbol": target["gene_symbol"].strip(),
                    "evidence_tier": species["evidence_tier"].strip(),
                    "reconstruction_status": "planned",
                    "sequence_length": "0",
                    "confidence_tier": "medium" if species["evidence_tier"].strip() == "A" else "low",
                    "method": "assembly_download_busco_projection_scaffold",
                    "notes": notes,
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
        outdir / "reference_assets.tsv",
        reference_assets_rows,
        [
            "entity_type",
            "entity_id",
            "scientific_name",
            "assembly_accession",
            "asset_status",
            "genome_fasta",
            "annotation_gtf",
            "protein_fasta",
            "log_path",
        ],
    )
    write_tsv(
        outdir / "projection_plan.tsv",
        projection_rows,
        [
            "species_id",
            "scientific_name",
            "reference_id",
            "projection_status",
            "plan_script",
            "plan_metadata",
            "missing_prerequisites",
        ],
    )


if __name__ == "__main__":
    main()
