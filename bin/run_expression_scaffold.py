#!/usr/bin/env python3

import argparse
import csv
import json
import shlex
import shutil
import ssl
import subprocess
from pathlib import Path
from typing import Dict, List, Optional
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen
from zipfile import ZipFile


ENA_PORTAL_BASE = "https://www.ebi.ac.uk/ena/portal/api/search"
ENA_READ_RUN_FIELDS = [
    "run_accession",
    "library_layout",
    "library_source",
    "library_strategy",
    "library_selection",
    "study_accession",
    "experiment_accession",
    "sample_accession",
    "instrument_platform",
    "instrument_model",
    "fastq_ftp",
    "read_count",
    "base_count",
]
NCBI_EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plan or execute cross-species RNA-seq quantification and expression study setup."
    )
    parser.add_argument("--species-manifest", required=True)
    parser.add_argument("--reference-manifest", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument(
        "--run-metadata",
        default="",
        help="Optional run_metadata.tsv emitted by build_species_manifest.py. Defaults to a sibling file next to species_manifest.tsv when present.",
    )
    parser.add_argument("--mode", choices=["plan", "scaffold", "execute"], default="scaffold")
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--ncbi-api-key", default="")
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


def find_first(root: Path, patterns: List[str]) -> str:
    for pattern in patterns:
        matches = sorted(root.rglob(pattern))
        if matches:
            return str(matches[0].resolve())
    return ""


def stage_ncbi_datasets(accession: str, stage_dir: Path, api_key: str) -> Dict[str, str]:
    zip_path = stage_dir / "ncbi_dataset.zip"
    unpack_dir = stage_dir / "ncbi_dataset"
    log_path = stage_dir / "datasets_download.log"
    stage_dir.mkdir(parents=True, exist_ok=True)

    datasets_exe = shutil.which("datasets")
    if not datasets_exe:
        return {
            "status": "datasets_missing",
            "genome_fasta": "",
            "annotation_gtf": "",
            "annotation_gff3": "",
            "protein_fasta": "",
            "zip_path": str(zip_path),
            "log_path": str(log_path),
        }

    if not unpack_dir.exists():
        cmd = [
            datasets_exe,
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
            code = run_command(cmd, cwd=stage_dir, log_path=log_path)
            if code == 0:
                break
            zip_path.unlink(missing_ok=True)
            log_path.write_text(log_path.read_text() + f"\n[retry] datasets attempt {attempt} failed\n")
        if code != 0:
            return {
                "status": "datasets_failed",
                "genome_fasta": "",
                "annotation_gtf": "",
                "annotation_gff3": "",
                "protein_fasta": "",
                "zip_path": str(zip_path),
                "log_path": str(log_path),
            }
        with ZipFile(zip_path) as archive:
            archive.extractall(unpack_dir)

    return {
        "status": "downloaded",
        "genome_fasta": find_first(unpack_dir, ["*_genomic.fna", "*.fna"]),
        "annotation_gtf": find_first(unpack_dir, ["*.gtf"]),
        "annotation_gff3": find_first(unpack_dir, ["*.gff3", "*.gff"]),
        "protein_fasta": find_first(unpack_dir, ["*protein.faa", "*.faa"]),
        "zip_path": str(zip_path),
        "log_path": str(log_path),
    }


def resolve_run_metadata_path(species_manifest_path: Path, explicit_path: str) -> Optional[Path]:
    if explicit_path:
        candidate = Path(explicit_path)
        return candidate if candidate.exists() else None
    sibling = species_manifest_path.parent / "run_metadata.tsv"
    return sibling if sibling.exists() else None


def load_run_metadata(path: Optional[Path]) -> Dict[str, Dict[str, str]]:
    if path is None or not path.exists():
        return {}
    return {row["run_accession"].strip(): row for row in read_tsv(path) if row.get("run_accession", "").strip()}


def parse_fastq_urls(value: str) -> List[str]:
    urls = []
    for raw in split_accessions(value.replace(";", ",")):
        if raw.startswith("http://") or raw.startswith("https://") or raw.startswith("ftp://"):
            urls.append(raw)
        elif raw:
            urls.append(f"https://{raw}")
    return urls


def shell_join(parts: List[str]) -> str:
    return " ".join(shlex.quote(part) for part in parts)


def build_request(url: str) -> Request:
    return Request(
        url,
        headers={
            "User-Agent": "phylloscopus-comparative-ortholog-scaffold/0.1 (+expression scaffold)",
            "Accept": "application/json, text/plain, */*",
        },
    )


def fetch_text(url: str, ssl_context, attempts: int = 3) -> str:
    last_error = None
    for _ in range(attempts):
        try:
            with urlopen(build_request(url), timeout=60, context=ssl_context) as response:
                return response.read().decode("utf-8")
        except (HTTPError, URLError) as error:
            last_error = error
    if last_error:
        raise last_error
    return ""


def fetch_json(url: str, ssl_context):
    payload = fetch_text(url, ssl_context)
    return json.loads(payload) if payload else None


def normalize_layout(value: str) -> str:
    layout = (value or "").strip().upper()
    if layout in {"SINGLE", "PAIRED"}:
        return layout
    return "unknown"


def merge_metadata(base: Dict[str, str], overlay: Dict[str, str], source_label: str) -> Dict[str, str]:
    merged = dict(base)
    if not overlay:
        return merged

    for key, value in overlay.items():
        if key == "metadata_source":
            continue
        value = (value or "").strip()
        if not value:
            continue
        if key == "library_layout":
            current = normalize_layout(merged.get(key, ""))
            candidate = normalize_layout(value)
            if current == "unknown" and candidate != "unknown":
                merged[key] = candidate
                merged["layout_source"] = source_label
            elif key not in merged or not (merged.get(key) or "").strip():
                merged[key] = candidate
        elif not (merged.get(key) or "").strip():
            merged[key] = value

    sources = [item for item in split_accessions(merged.get("metadata_source", "").replace(";", ",")) if item]
    if source_label not in sources:
        sources.append(source_label)
    merged["metadata_source"] = ",".join(sources)
    merged["library_layout"] = normalize_layout(merged.get("library_layout", ""))
    if merged["library_layout"] == "unknown" and not merged.get("layout_source"):
        merged["layout_source"] = ""
    return merged


def fetch_ena_run_metadata(accession: str, ssl_context) -> Dict[str, str]:
    params = {
        "result": "read_run",
        "fields": ",".join(ENA_READ_RUN_FIELDS),
        "query": f'run_accession="{accession}"',
        "format": "json",
        "limit": 1,
    }
    rows = fetch_json(f"{ENA_PORTAL_BASE}?{urlencode(params)}", ssl_context) or []
    if not rows:
        return {}

    row = rows[0]
    return {
        "run_accession": (row.get("run_accession") or accession).strip(),
        "library_layout": normalize_layout(row.get("library_layout", "")),
        "library_source": (row.get("library_source") or "").strip(),
        "library_strategy": (row.get("library_strategy") or "").strip(),
        "library_selection": (row.get("library_selection") or "").strip(),
        "study_accession": (row.get("study_accession") or "").strip(),
        "experiment_accession": (row.get("experiment_accession") or "").strip(),
        "sample_accession": (row.get("sample_accession") or "").strip(),
        "instrument_platform": (row.get("instrument_platform") or "").strip(),
        "instrument_model": (row.get("instrument_model") or "").strip(),
        "fastq_ftp": (row.get("fastq_ftp") or "").strip(),
        "read_count": (row.get("read_count") or "").strip(),
        "base_count": (row.get("base_count") or "").strip(),
        "metadata_source": "ena",
        "layout_source": "ena" if normalize_layout(row.get("library_layout", "")) != "unknown" else "",
    }


def fetch_ncbi_runinfo_metadata(accession: str, ssl_context) -> Dict[str, str]:
    params = {
        "db": "sra",
        "id": accession,
        "rettype": "runinfo",
        "retmode": "text",
    }
    payload = fetch_text(f"{NCBI_EUTILS_BASE}?{urlencode(params)}", ssl_context)
    reader = csv.DictReader(payload.splitlines())
    row = next(reader, None)
    if not row:
        return {}

    return {
        "run_accession": (row.get("Run") or accession).strip(),
        "library_layout": normalize_layout(row.get("LibraryLayout", "")),
        "library_source": (row.get("LibrarySource") or "").strip(),
        "library_strategy": (row.get("LibraryStrategy") or "").strip(),
        "library_selection": (row.get("LibrarySelection") or "").strip(),
        "study_accession": (row.get("BioProject") or row.get("SRAStudy") or "").strip(),
        "experiment_accession": (row.get("Experiment") or "").strip(),
        "sample_accession": (row.get("BioSample") or row.get("Sample") or "").strip(),
        "instrument_platform": (row.get("Platform") or "").strip(),
        "instrument_model": (row.get("Model") or "").strip(),
        "read_count": (row.get("spots") or "").strip(),
        "base_count": (row.get("bases") or "").strip(),
        "metadata_source": "ncbi_runinfo",
        "layout_source": "ncbi_runinfo" if normalize_layout(row.get("LibraryLayout", "")) != "unknown" else "",
    }


def enrich_run_metadata(run_metadata: Dict[str, Dict[str, str]], accessions: List[str]) -> Dict[str, Dict[str, str]]:
    ssl_context = ssl.create_default_context()
    enriched = dict(run_metadata)

    for accession in accessions:
        seed = dict(enriched.get(accession, {}))
        current = dict(seed)
        current.setdefault("run_accession", accession)
        current["library_layout"] = normalize_layout(current.get("library_layout", ""))
        has_local_metadata = any(
            key != "run_accession" and (value or "").strip()
            for key, value in seed.items()
        )
        if current.get("metadata_source", "").strip() == "" and has_local_metadata:
            current["metadata_source"] = "manifest"
        if current["library_layout"] in {"SINGLE", "PAIRED"}:
            current.setdefault("layout_source", "manifest")
            enriched[accession] = current
            continue

        ena_row = {}
        try:
            ena_row = fetch_ena_run_metadata(accession, ssl_context)
        except Exception:
            ena_row = {}
        current = merge_metadata(current, ena_row, "ena") if ena_row else current

        if normalize_layout(current.get("library_layout", "")) == "unknown":
            ncbi_row = {}
            try:
                ncbi_row = fetch_ncbi_runinfo_metadata(accession, ssl_context)
            except Exception:
                ncbi_row = {}
            current = merge_metadata(current, ncbi_row, "ncbi_runinfo") if ncbi_row else current

        current["library_layout"] = normalize_layout(current.get("library_layout", ""))
        enriched[accession] = current

    return enriched


def ensure_reference_assets(
    primary_reference: Dict[str, str],
    ref_dir: Path,
    mode: str,
    threads: int,
    ncbi_api_key: str,
) -> Dict[str, str]:
    datasets_status = "not_needed"
    transcript_build_status = "planned"
    salmon_index_status = "planned"

    transcript_fasta = resolve_local_path(primary_reference.get("transcript_fasta", ""))
    assembly_fasta = resolve_local_path(primary_reference.get("assembly_fasta", ""))
    annotation_gtf = resolve_local_path(primary_reference.get("annotation_gtf", ""))
    annotation_source = annotation_gtf

    if mode in {"scaffold", "execute"} and (not assembly_fasta or not annotation_source) and primary_reference.get("assembly_accession", "").strip():
        downloaded = stage_ncbi_datasets(
            primary_reference.get("assembly_accession", "").strip(),
            ref_dir / "ncbi_dataset",
            ncbi_api_key,
        )
        datasets_status = downloaded["status"]
        if not assembly_fasta:
            assembly_fasta = downloaded["genome_fasta"]
        if not annotation_gtf:
            annotation_gtf = downloaded["annotation_gtf"]
        annotation_source = annotation_gtf or downloaded["annotation_gff3"]
    elif primary_reference.get("assembly_accession", "").strip():
        datasets_status = "planned"

    gffread_exe = shutil.which("gffread")
    if transcript_fasta:
        transcript_build_status = "provided"
    elif mode in {"scaffold", "execute"} and assembly_fasta and annotation_source and gffread_exe:
        transcript_fasta = str((ref_dir / "primary_reference_transcripts.fa").resolve())
        transcript_build_log = ref_dir / "build_transcripts.log"
        code = run_command(
            [gffread_exe, "-w", transcript_fasta, "-g", assembly_fasta, annotation_source],
            cwd=ref_dir,
            log_path=transcript_build_log,
        )
        if code == 0 and Path(transcript_fasta).exists():
            transcript_build_status = "built_with_gffread"
        else:
            transcript_fasta = ""
            transcript_build_status = "gffread_failed"
    elif not gffread_exe and assembly_fasta and annotation_source:
        transcript_build_status = "gffread_missing"
    elif not assembly_fasta:
        transcript_build_status = "assembly_fasta_missing"
    else:
        transcript_build_status = "annotation_missing"

    salmon_exe = shutil.which("salmon")
    salmon_index_dir = ref_dir / "salmon_index"
    salmon_index_log = ref_dir / "salmon_index.log"
    if transcript_fasta and mode in {"scaffold", "execute"} and salmon_exe:
        code = run_command(
            [salmon_exe, "index", "-t", transcript_fasta, "-i", str(salmon_index_dir), "-p", str(max(threads, 1))],
            cwd=ref_dir,
            log_path=salmon_index_log,
        )
        salmon_index_status = "completed" if code == 0 else "failed"
    elif not transcript_fasta:
        salmon_index_status = "missing_transcript_fasta"
    elif not salmon_exe:
        salmon_index_status = "salmon_missing"

    return {
        "reference_id": primary_reference.get("reference_id", "").strip(),
        "scientific_name": primary_reference.get("scientific_name", "").strip(),
        "assembly_accession": primary_reference.get("assembly_accession", "").strip(),
        "assembly_fasta": assembly_fasta,
        "annotation_gtf": annotation_gtf,
        "annotation_source": annotation_source,
        "transcript_fasta": transcript_fasta,
        "datasets_status": datasets_status,
        "transcript_build_status": transcript_build_status,
        "salmon_index_status": salmon_index_status,
        "salmon_index_dir": str(salmon_index_dir.resolve()) if salmon_index_dir.exists() else str(salmon_index_dir.resolve()),
    }


def build_download_commands(
    accession: str,
    reads_dir: Path,
    fastq_urls: List[str],
    layout: str,
    threads: int,
) -> tuple[List[str], str, List[str]]:
    commands: List[str] = []
    read_paths: List[str] = []
    curl_exe = shutil.which("curl")
    prefetch_exe = shutil.which("prefetch")
    fasterq_exe = shutil.which("fasterq-dump")

    if fastq_urls and curl_exe:
        commands.append(f"mkdir -p {shlex.quote(str(reads_dir))}")
        for url in fastq_urls:
            filename = Path(url).name
            destination = reads_dir / filename
            commands.append(
                shell_join([curl_exe, "-L", "--retry", "3", "-o", str(destination), url])
            )
            read_paths.append(str(destination))
        return commands, "ena_fastq", read_paths

    if prefetch_exe and fasterq_exe:
        sra_root = reads_dir.parent.parent / "sra"
        sra_path = sra_root / accession / f"{accession}.sra"
        commands.append(f"mkdir -p {shlex.quote(str(reads_dir))} {shlex.quote(str(sra_root))}")
        commands.append(shell_join([prefetch_exe, accession, "--output-directory", str(sra_root)]))
        fasterq_cmd = [fasterq_exe, str(sra_path), "-O", str(reads_dir), "--threads", str(max(threads, 1))]
        if layout == "PAIRED":
            fasterq_cmd.append("--split-files")
            read_paths = [
                str(reads_dir / f"{accession}_1.fastq"),
                str(reads_dir / f"{accession}_2.fastq"),
            ]
        elif layout == "SINGLE":
            read_paths = [str(reads_dir / f"{accession}.fastq")]
        commands.append(shell_join(fasterq_cmd))
        return commands, "sra_tools", read_paths

    return commands, "unavailable", read_paths


def build_salmon_quant_command(
    layout: str,
    salmon_index_dir: str,
    read_paths: List[str],
    quant_out: Path,
    threads: int,
) -> Optional[List[str]]:
    salmon_exe = shutil.which("salmon")
    if not salmon_exe:
        return None
    if layout == "PAIRED" and len(read_paths) >= 2:
        return [
            salmon_exe,
            "quant",
            "-i",
            salmon_index_dir,
            "-l",
            "A",
            "-1",
            read_paths[0],
            "-2",
            read_paths[1],
            "--validateMappings",
            "-o",
            str(quant_out),
            "-p",
            str(max(threads, 1)),
        ]
    if layout == "SINGLE" and read_paths:
        return [
            salmon_exe,
            "quant",
            "-i",
            salmon_index_dir,
            "-l",
            "A",
            "-r",
            read_paths[0],
            "--validateMappings",
            "-o",
            str(quant_out),
            "-p",
            str(max(threads, 1)),
        ]
    return None


def main():
    args = parse_args()
    species_manifest_path = Path(args.species_manifest)
    reference_manifest_path = Path(args.reference_manifest)
    outdir = Path(args.outdir)
    ref_dir = outdir / "reference"
    plans_dir = outdir / "plans"
    quant_dir = outdir / "quant"
    outdir.mkdir(parents=True, exist_ok=True)
    ref_dir.mkdir(parents=True, exist_ok=True)
    plans_dir.mkdir(parents=True, exist_ok=True)
    quant_dir.mkdir(parents=True, exist_ok=True)

    species_rows = read_tsv(species_manifest_path)
    reference_rows = read_tsv(reference_manifest_path)
    primary_reference = next(
        (row for row in reference_rows if row.get("reference_role", "").strip() == "primary"),
        reference_rows[0] if reference_rows else None,
    )
    if primary_reference is None:
        raise SystemExit("reference_manifest.tsv must contain at least one reference.")

    run_metadata_path = resolve_run_metadata_path(species_manifest_path, args.run_metadata)
    run_metadata = load_run_metadata(run_metadata_path)
    manifest_run_accessions = sorted(
        {
            accession
            for species in species_rows
            for accession in split_accessions(species.get("rna_sra_accessions", ""))
            if accession
        }
    )
    run_metadata = enrich_run_metadata(run_metadata, manifest_run_accessions)
    reference_assets = ensure_reference_assets(
        primary_reference=primary_reference,
        ref_dir=ref_dir,
        mode=args.mode,
        threads=args.threads,
        ncbi_api_key=args.ncbi_api_key,
    )

    sample_rows = []
    plan_rows = []
    design_rows = []

    for species in species_rows:
        scientific_name = species.get("scientific_name", "").strip()
        species_id = species.get("species_id", "").strip()
        run_accessions = split_accessions(species.get("rna_sra_accessions", ""))
        if not run_accessions:
            continue

        for accession in run_accessions:
            metadata = run_metadata.get(accession, {})
            sample_id = f"{species_id}__{accession}"
            layout = normalize_layout(metadata.get("library_layout", ""))
            library_type = "A" if layout in {"SINGLE", "PAIRED"} else "unknown"
            fastq_urls = parse_fastq_urls(metadata.get("fastq_ftp", ""))

            sample_rows.append(
                {
                    "sample_id": sample_id,
                    "species_id": species_id,
                    "scientific_name": scientific_name,
                    "run_accession": accession,
                    "library_layout": layout,
                    "library_type": library_type,
                    "study_accession": (metadata.get("study_accession") or "").strip(),
                    "experiment_accession": (metadata.get("experiment_accession") or "").strip(),
                    "sample_accession": (metadata.get("sample_accession") or "").strip(),
                    "instrument_platform": (metadata.get("instrument_platform") or "").strip(),
                    "instrument_model": (metadata.get("instrument_model") or "").strip(),
                    "library_selection": (metadata.get("library_selection") or "").strip(),
                    "read_count": (metadata.get("read_count") or "").strip(),
                    "base_count": (metadata.get("base_count") or "").strip(),
                    "metadata_source": (metadata.get("metadata_source") or "").strip(),
                    "layout_source": (metadata.get("layout_source") or "").strip(),
                    "condition": species_id,
                    "contrast_group": species_id,
                    "batch": "",
                    "notes": "populate biological condition and batch metadata before DE analysis",
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
            sample_log = sample_plan_dir / "run_expression_quant.log"

            download_commands, download_method, read_paths = build_download_commands(
                accession=accession,
                reads_dir=reads_dir,
                fastq_urls=fastq_urls,
                layout=layout,
                threads=args.threads,
            )
            salmon_cmd = build_salmon_quant_command(
                layout=layout,
                salmon_index_dir=reference_assets["salmon_index_dir"],
                read_paths=read_paths,
                quant_out=quant_out,
                threads=args.threads,
            )

            commands = ["#!/usr/bin/env bash", "set -euo pipefail", ""]
            if download_commands:
                commands.extend(download_commands)
            else:
                commands.append(f"# No download method available for {accession}")
            commands.append("")
            if salmon_cmd:
                commands.append(shell_join(salmon_cmd))
            else:
                commands.append("# Quantification prerequisites missing or library layout unresolved")
            plan_script.write_text("\n".join(commands) + "\n")
            plan_script.chmod(0o755)

            notes = []
            execution_ready = True
            if download_method == "unavailable":
                notes.append("download_method_missing")
                execution_ready = False
            if layout not in {"SINGLE", "PAIRED"}:
                notes.append("library_layout_unknown")
                execution_ready = False
            if reference_assets["transcript_build_status"] not in {"provided", "built_with_gffread"}:
                notes.append(reference_assets["transcript_build_status"])
                execution_ready = False
            if reference_assets["salmon_index_status"] != "completed":
                notes.append(reference_assets["salmon_index_status"])
                execution_ready = False
            if not salmon_cmd:
                execution_ready = False
                if "salmon_missing" not in notes and shutil.which("salmon") is None:
                    notes.append("salmon_missing")
            if not metadata:
                notes.append("run_metadata_unresolved")
                execution_ready = False

            quant_status = "planned"
            if args.mode == "execute" and execution_ready:
                code = run_command(["bash", str(plan_script.resolve())], cwd=sample_plan_dir, log_path=sample_log)
                quant_status = "completed" if code == 0 else "failed"
                if code != 0:
                    notes.append("execution_failed")
            elif execution_ready:
                quant_status = "ready"

            plan_rows.append(
                {
                    "sample_id": sample_id,
                    "species_id": species_id,
                    "run_accession": accession,
                    "library_layout": layout,
                    "download_method": download_method,
                    "plan_script": str(plan_script),
                    "log_path": str(sample_log),
                    "quant_status": quant_status,
                    "notes": "; ".join(sorted(dict.fromkeys(notes))) if notes else "ready_to_quantify",
                }
            )

    reference_rows_out = [
        {
            "reference_id": reference_assets["reference_id"],
            "scientific_name": reference_assets["scientific_name"],
            "assembly_accession": reference_assets["assembly_accession"],
            "transcript_fasta": reference_assets["transcript_fasta"],
            "assembly_fasta": reference_assets["assembly_fasta"],
            "annotation_gtf": reference_assets["annotation_gtf"],
            "annotation_source": reference_assets["annotation_source"],
            "datasets_status": reference_assets["datasets_status"],
            "transcript_build_status": reference_assets["transcript_build_status"],
            "salmon_index_status": reference_assets["salmon_index_status"],
        }
    ]

    summary_rows = [
        {"metric": "rna_sample_count", "value": str(len(sample_rows))},
        {"metric": "species_with_rna_runs", "value": str(len({row['species_id'] for row in sample_rows}))},
        {"metric": "run_metadata_path", "value": str(run_metadata_path) if run_metadata_path else ""},
        {"metric": "run_metadata_rows_loaded", "value": str(len(run_metadata))},
        {"metric": "datasets_status", "value": reference_assets["datasets_status"]},
        {"metric": "transcript_build_status", "value": reference_assets["transcript_build_status"]},
        {"metric": "salmon_index_status", "value": reference_assets["salmon_index_status"]},
        {"metric": "quant_ready", "value": str(sum(1 for row in plan_rows if row["quant_status"] == "ready"))},
        {"metric": "quant_completed", "value": str(sum(1 for row in plan_rows if row["quant_status"] == "completed"))},
        {
            "metric": "layouts_resolved",
            "value": str(sum(1 for row in sample_rows if row["library_layout"] in {"SINGLE", "PAIRED"})),
        },
    ]

    write_tsv(
        outdir / "expression_samples.tsv",
        sample_rows,
        [
            "sample_id",
            "species_id",
            "scientific_name",
            "run_accession",
            "library_layout",
            "library_type",
            "study_accession",
            "experiment_accession",
            "sample_accession",
            "instrument_platform",
            "instrument_model",
            "library_selection",
            "read_count",
            "base_count",
            "metadata_source",
            "layout_source",
            "condition",
            "contrast_group",
            "batch",
            "notes",
        ],
    )
    write_tsv(
        outdir / "expression_design_template.tsv",
        design_rows,
        ["sample_id", "condition", "batch", "contrast_group", "notes"],
    )
    write_tsv(
        outdir / "expression_plan.tsv",
        plan_rows,
        ["sample_id", "species_id", "run_accession", "library_layout", "download_method", "plan_script", "log_path", "quant_status", "notes"],
    )
    write_tsv(
        outdir / "expression_reference.tsv",
        reference_rows_out,
        [
            "reference_id",
            "scientific_name",
            "assembly_accession",
            "transcript_fasta",
            "assembly_fasta",
            "annotation_gtf",
            "annotation_source",
            "datasets_status",
            "transcript_build_status",
            "salmon_index_status",
        ],
    )
    write_tsv(outdir / "expression_summary.tsv", summary_rows, ["metric", "value"])


if __name__ == "__main__":
    main()
