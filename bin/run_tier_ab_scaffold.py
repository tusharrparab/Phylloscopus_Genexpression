#!/usr/bin/env python3

import argparse
import csv
import gzip
import json
import re
import shutil
import subprocess
import urllib.request
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple
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


def normalize_row(row: Dict[str, str], fieldnames: List[str]) -> Dict[str, str]:
    normalized = {}
    for field in fieldnames:
        value = row.get(field, "")
        normalized[field] = "" if value is None else str(value)
    return normalized


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


def http_get_json(url: str) -> Dict[str, object]:
    with urllib.request.urlopen(url) as response:
        return json.loads(response.read().decode("utf-8"))


def download_file(url: str, destination: Path):
    destination.parent.mkdir(parents=True, exist_ok=True)
    if shutil.which("curl"):
        subprocess.run(["curl", "-L", url, "-o", str(destination)], check=True)
        return
    with urllib.request.urlopen(url) as response, destination.open("wb") as handle:
        shutil.copyfileobj(response, handle)


def stage_ncbi_http_fallback(accession: str, stage_dir: Path) -> Dict[str, str]:
    log_path = stage_dir / "http_download.log"
    stage_dir.mkdir(parents=True, exist_ok=True)

    try:
        search_url = (
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
            f"?db=assembly&term={accession}%5BAssembly%20Accession%5D&retmode=json"
        )
        search_data = http_get_json(search_url)
        id_list = search_data.get("esearchresult", {}).get("idlist", [])
        if not id_list:
            raise ValueError(f"Assembly accession {accession} was not found in NCBI Assembly")

        summary_url = (
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
            f"?db=assembly&id={id_list[0]}&retmode=json"
        )
        summary_data = http_get_json(summary_url)
        result = summary_data.get("result", {}).get(id_list[0], {})
        ftp_path = (result.get("ftppath_genbank") or "").strip()
        if not ftp_path:
            raise ValueError(f"NCBI Assembly summary for {accession} did not include a GenBank FTP path")

        prefix = ftp_path.rstrip("/").split("/")[-1]
        genome_gz = stage_dir / f"{prefix}_genomic.fna.gz"
        genome_fasta = stage_dir / f"{prefix}_genomic.fna"
        https_path = ftp_path.replace("ftp://", "https://")
        genome_url = f"{https_path}/{prefix}_genomic.fna.gz"

        if not genome_fasta.exists():
            if not genome_gz.exists():
                download_file(genome_url, genome_gz)
            with gzip.open(genome_gz, "rb") as src, genome_fasta.open("wb") as dst:
                shutil.copyfileobj(src, dst)

        log_path.write_text(json.dumps({"accession": accession, "ftppath_genbank": ftp_path}, indent=2) + "\n")
        return {
            "status": "downloaded_http_fallback",
            "genome_fasta": str(genome_fasta.resolve()),
            "annotation_gtf": "",
            "protein_fasta": "",
            "gff3": "",
            "zip_path": str(genome_gz.resolve()),
            "log_path": str(log_path.resolve()),
        }
    except Exception as exc:
        log_path.write_text(f"{type(exc).__name__}: {exc}\n")
        return {
            "status": "http_download_failed",
            "genome_fasta": "",
            "annotation_gtf": "",
            "protein_fasta": "",
            "gff3": "",
            "zip_path": "",
            "log_path": str(log_path.resolve()),
        }


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
        return stage_ncbi_http_fallback(accession, stage_dir)

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


def read_fasta_dict(path: Path) -> Dict[str, str]:
    records = {}
    header = None
    chunks: List[str] = []
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records[header] = "".join(chunks).upper()
                header = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if header is not None:
        records[header] = "".join(chunks).upper()
    return records


def fasta_iter(path: Path) -> Iterable[Tuple[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    try:
        with opener(path, "rt") as handle:
            header = None
            chunks: List[str] = []
            for raw_line in handle:
                line = raw_line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if header is not None:
                        yield header, "".join(chunks).upper()
                    header = line[1:].split()[0]
                    chunks = []
                else:
                    chunks.append(line)
            if header is not None:
                yield header, "".join(chunks).upper()
    except (EOFError, gzip.BadGzipFile):
        return


def gene_id_from_transcript_id(transcript_id: str) -> str:
    return re.split(r"[-|]", transcript_id, maxsplit=1)[0].upper()


IUPAC = {
    "A": {"A"},
    "C": {"C"},
    "G": {"G"},
    "T": {"T"},
    "R": {"A", "G"},
    "Y": {"C", "T"},
    "S": {"C", "G"},
    "W": {"A", "T"},
    "K": {"G", "T"},
    "M": {"A", "C"},
    "B": {"C", "G", "T"},
    "D": {"A", "G", "T"},
    "H": {"A", "C", "T"},
    "V": {"A", "C", "G"},
    "N": {"A", "C", "G", "T"},
}

RC_TABLE = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(RC_TABLE)[::-1]


def bases_match(query_base: str, target_base: str) -> bool:
    return target_base in IUPAC.get(query_base, {query_base})


def build_seeds(query: str, k: int = 23, max_seeds: int = 3) -> List[Tuple[int, str]]:
    seeds = []
    upper = query.upper()
    if len(upper) < k:
        return seeds

    span = len(upper) - k
    if max_seeds <= 1 or span == 0:
        candidate_offsets = [0]
    else:
        candidate_offsets = sorted({round(index * span / (max_seeds - 1)) for index in range(max_seeds)})

    for offset in candidate_offsets:
        seed = upper[offset:offset + k]
        if len(seed) == k and set(seed) <= {"A", "C", "G", "T"}:
            seeds.append((offset, seed))
    if not seeds:
        fallback_k = 13
        fallback_span = len(upper) - fallback_k
        fallback_offsets = [0] if fallback_span <= 0 else sorted(
            {round(index * fallback_span / 3) for index in range(4)}
        )
        for offset in fallback_offsets:
            seed = upper[offset:offset + fallback_k]
            if len(seed) == fallback_k and set(seed) <= {"A", "C", "G", "T"}:
                seeds.append((offset, seed))
    return seeds


def extend_ungapped_hit(query: str, target: str, query_start: int, target_start: int, seed_len: int) -> Dict[str, object]:
    query_end = query_start + seed_len
    target_end = target_start + seed_len

    best_left_q = query_start
    best_left_t = target_start
    running = 0
    best_score = 0
    xdrop = 30
    q_index = query_start - 1
    t_index = target_start - 1
    while q_index >= 0 and t_index >= 0:
        running += 2 if bases_match(query[q_index], target[t_index]) else -3
        if running > best_score:
            best_score = running
            best_left_q = q_index
            best_left_t = t_index
        if best_score - running > xdrop:
            break
        q_index -= 1
        t_index -= 1

    best_right_q = query_end
    best_right_t = target_end
    running = 0
    best_score = 0
    q_index = query_end
    t_index = target_end
    while q_index < len(query) and t_index < len(target):
        running += 2 if bases_match(query[q_index], target[t_index]) else -3
        if running > best_score:
            best_score = running
            best_right_q = q_index + 1
            best_right_t = t_index + 1
        if best_score - running > xdrop:
            break
        q_index += 1
        t_index += 1

    aligned_query = query[best_left_q:best_right_q]
    aligned_target = target[best_left_t:best_right_t]
    matches = sum(1 for q_base, t_base in zip(aligned_query, aligned_target) if bases_match(q_base, t_base))
    aligned_length = len(aligned_query)
    mismatches = aligned_length - matches
    score = matches * 2 - mismatches * 3
    identity = matches / aligned_length if aligned_length else 0.0
    coverage = aligned_length / len(query) if query else 0.0
    return {
        "query_start": best_left_q,
        "query_end": best_right_q,
        "target_start": best_left_t,
        "target_end": best_right_t,
        "aligned_length": aligned_length,
        "identity": identity,
        "coverage": coverage,
        "score": score,
    }


def resolve_target_query_path(target_row: Dict[str, str], ortholog_targets_path: Path) -> Optional[Path]:
    raw_path = (target_row.get("candidate_query_fasta") or "").strip()
    if not raw_path:
        return None
    direct = Path(raw_path)
    if direct.exists():
        return direct.resolve()
    relative = (ortholog_targets_path.resolve().parent / raw_path)
    if relative.exists():
        return relative.resolve()
    return None


def recover_candidate_sequence(query_path: Path, genome_fasta: Path) -> Optional[Dict[str, object]]:
    query_records = read_fasta_dict(query_path)
    if not query_records:
        return None

    query_name, query_sequence = next(iter(query_records.items()))
    query_sequence = query_sequence.upper()
    orientations = [
        ("+", query_sequence, build_seeds(query_sequence)),
        ("-", reverse_complement(query_sequence), build_seeds(reverse_complement(query_sequence))),
    ]
    best_hit = None

    for contig_id, contig_sequence in fasta_iter(genome_fasta):
        for strand, oriented_query, seeds in orientations:
            for query_offset, seed in seeds:
                search_start = 0
                while True:
                    position = contig_sequence.find(seed, search_start)
                    if position < 0:
                        break
                    candidate = extend_ungapped_hit(oriented_query, contig_sequence, query_offset, position, len(seed))
                    candidate.update(
                        {
                            "contig_id": contig_id,
                            "strand": strand,
                            "query_name": query_name,
                            "query_length": len(oriented_query),
                        }
                    )
                    if candidate["aligned_length"] >= 500 and candidate["identity"] >= 0.75:
                        if candidate["identity"] >= 0.98 and candidate["coverage"] >= 0.95:
                            candidate_sequence = contig_sequence[candidate["target_start"]:candidate["target_end"]]
                            if strand == "-":
                                candidate_sequence = reverse_complement(candidate_sequence)
                            candidate["sequence"] = candidate_sequence
                            candidate["start_1based"] = int(candidate["target_start"]) + 1
                            candidate["end_1based"] = int(candidate["target_end"])
                            return candidate
                        if best_hit is None or (candidate["score"], candidate["aligned_length"]) > (
                            best_hit["score"],
                            best_hit["aligned_length"],
                        ):
                            best_hit = candidate
                    search_start = position + 1

    if best_hit is None:
        return None

    best_sequence = ""
    for contig_id, contig_sequence in fasta_iter(genome_fasta):
        if contig_id != best_hit["contig_id"]:
            continue
        best_sequence = contig_sequence[best_hit["target_start"]:best_hit["target_end"]]
        if best_hit["strand"] == "-":
            best_sequence = reverse_complement(best_sequence)
        break

    if not best_sequence:
        return None

    best_hit["sequence"] = best_sequence
    best_hit["start_1based"] = int(best_hit["target_start"]) + 1
    best_hit["end_1based"] = int(best_hit["target_end"])
    return best_hit


STATUS_FIELDS = [
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
    "candidate_sequence",
    "assembly_backed",
    "orthology_unvalidated",
    "recovery_route",
    "query_source",
    "source_coordinates",
    "query_coverage",
    "match_identity",
]


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
    primary_transcripts = {}
    primary_transcript_path = resolve_local_asset(primary_reference, "transcript_fasta")
    if primary_transcript_path:
        primary_transcripts = read_fasta_dict(Path(primary_transcript_path))
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
            gene_id = target["gene_id"].strip()
            sequence = ""
            reconstruction_status = "planned"
            method = "assembly_download_busco_projection_scaffold"
            confidence_tier = "medium" if species["evidence_tier"].strip() == "A" else "low"
            candidate_sequence = "no"
            assembly_backed = "yes" if assets["genome_fasta"] else "no"
            orthology_unvalidated = "no" if species["evidence_tier"].strip() == "A" else "yes"
            recovery_route = species["evidence_tier"].strip()
            query_source = ""
            source_coordinates = ""
            query_coverage = ""
            match_identity = ""
            row_notes = notes

            if species["evidence_tier"].strip() == "A":
                transcript_match = next(
                    (
                        transcript_sequence
                        for transcript_id, transcript_sequence in primary_transcripts.items()
                        if gene_id_from_transcript_id(transcript_id) == gene_id.upper()
                    ),
                    "",
                )
                if transcript_match:
                    sequence = transcript_match
                    reconstruction_status = "recovered"
                    method = "reference_backed_transcript_extraction"
                    confidence_tier = "high"
                    orthology_unvalidated = "no"
                    recovery_route = "tier_a_reference_transcript"
                    query_source = primary_transcript_path
                    row_notes = f"{notes}; transcript_extracted_from_primary_reference"
                else:
                    reconstruction_status = "missing_data"
                    row_notes = f"{notes}; transcript_missing_from_primary_reference"

            if species["evidence_tier"].strip() == "B" and assets["genome_fasta"]:
                query_path = resolve_target_query_path(target, Path(args.ortholog_targets))
                if query_path is not None:
                    candidate_hit = recover_candidate_sequence(query_path, Path(assets["genome_fasta"]))
                    if candidate_hit is not None:
                        sequence = str(candidate_hit["sequence"])
                        reconstruction_status = "candidate_sequence_recovered"
                        method = "assembly_backed_candidate_sequence_search"
                        confidence_tier = "low"
                        candidate_sequence = "yes"
                        assembly_backed = "yes"
                        orthology_unvalidated = "yes"
                        recovery_route = "tier_b_candidate_assembly_search"
                        query_source = str(query_path)
                        source_coordinates = (
                            f"{candidate_hit['contig_id']}:{candidate_hit['start_1based']}-{candidate_hit['end_1based']}"
                            f"({candidate_hit['strand']})"
                        )
                        query_coverage = f"{candidate_hit['coverage']:.3f}"
                        match_identity = f"{candidate_hit['identity']:.3f}"
                        row_notes = (
                            f"{notes}; candidate_sequence; assembly_backed; orthology_unvalidated; "
                            f"query={query_path.name}; query_record={candidate_hit['query_name']}"
                        )
                    else:
                        row_notes = f"{notes}; candidate_query_present_but_no_hit"

            if sequence:
                species_dir = ortholog_sequences_dir / species_id
                species_dir.mkdir(parents=True, exist_ok=True)
                species_fasta_path = species_dir / f"{species_id}.cds.fasta"
                header = (
                    f">{species_id}|{gene_id}|tier={species['evidence_tier'].strip()}|status={reconstruction_status}"
                    f"|route={recovery_route}"
                )
                with species_fasta_path.open("a") as handle:
                    handle.write(header + "\n")
                    for index in range(0, len(sequence), 60):
                        handle.write(sequence[index:index + 60] + "\n")

            status_rows.append(
                normalize_row(
                    {
                        "species_id": species_id,
                        "scientific_name": scientific_name,
                        "gene_id": gene_id,
                        "gene_symbol": target["gene_symbol"].strip(),
                        "target_category": target["category"].strip(),
                        "orthology_basis": target["orthology_basis"].strip(),
                        "copy_number_expectation": target["copy_number_expectation"].strip(),
                        "target_rationale": target["rationale"].strip(),
                        "evidence_tier": species["evidence_tier"].strip(),
                        "evidence_confidence": species.get("evidence_confidence", "").strip()
                        or ("high" if species["evidence_tier"].strip() == "A" else "medium"),
                        "reconstruction_status": reconstruction_status,
                        "sequence_length": len(sequence) if sequence else 0,
                        "confidence_tier": confidence_tier,
                        "method": method,
                        "notes": row_notes,
                        "candidate_sequence": candidate_sequence,
                        "assembly_backed": assembly_backed,
                        "orthology_unvalidated": orthology_unvalidated,
                        "recovery_route": recovery_route,
                        "query_source": query_source,
                        "source_coordinates": source_coordinates,
                        "query_coverage": query_coverage,
                        "match_identity": match_identity,
                    },
                    STATUS_FIELDS,
                )
            )

    write_tsv(
        outdir / "ortholog_status.tsv",
        status_rows,
        STATUS_FIELDS,
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
