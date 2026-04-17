#!/usr/bin/env python3

import argparse
import csv
import json
import os
import re
import ssl
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from urllib.error import HTTPError, URLError
from urllib.parse import quote, urlencode
from urllib.request import Request, urlopen


GBIF_BASE = "https://api.gbif.org/v1"
NCBI_EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ENA_PORTAL_BASE = "https://www.ebi.ac.uk/ena/portal/api"
ENA_TAXONOMY_BASE = "https://www.ebi.ac.uk/ena/taxonomy/rest"
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
    "fastq_bytes",
    "read_count",
    "base_count",
]

ASSEMBLY_LEVEL_RANK = {
    "Complete Genome": 0,
    "Chromosome": 1,
    "Scaffold": 2,
    "Contig": 3,
}

TRANSCRIPTOMIC_SOURCES = {"TRANSCRIPTOMIC"}
TRANSCRIPTOMIC_STRATEGIES = {"RNA-Seq", "TRANSCRIPTOME"}
GENOMIC_WGS_STRATEGIES = {"WGS", "WGA", "WXS"}

try:
    import certifi
except ImportError:  # pragma: no cover - environment-dependent
    certifi = None


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build a dated species manifest for a target genus using GBIF, NCBI, and ENA."
    )
    parser.add_argument("--genus", default="Phylloscopus", help="Target genus to snapshot.")
    parser.add_argument(
        "--outdir",
        required=True,
        help="Directory where species_manifest.tsv and metadata outputs will be written.",
    )
    parser.add_argument(
        "--authority",
        choices=["gbif"],
        default="gbif",
        help="Taxonomy authority to use for the species snapshot.",
    )
    parser.add_argument(
        "--species-limit",
        type=int,
        default=0,
        help="Optional limit for testing; 0 means all species returned by the authority.",
    )
    parser.add_argument(
        "--max-runs-per-class",
        type=int,
        default=10,
        help="Maximum ENA run accessions to write into rna_sra_accessions and wgs_sra_accessions.",
    )
    parser.add_argument(
        "--include-ena",
        action="store_true",
        help="Include ENA run classification to distinguish transcriptomic versus genomic runs.",
    )
    parser.add_argument(
        "--ncbi-api-key",
        default=os.environ.get("NCBI_API_KEY", ""),
        help="Optional NCBI API key. Defaults to NCBI_API_KEY from the environment.",
    )
    parser.add_argument(
        "--request-sleep-seconds",
        type=float,
        default=0.34,
        help="Delay between remote requests to stay polite to public APIs.",
    )
    parser.add_argument(
        "--allow-insecure-ssl",
        action="store_true",
        help="Disable certificate verification as a last-resort fallback for broken local Python SSL setups.",
    )
    return parser.parse_args()


def build_request(url: str) -> Request:
    return Request(
        url,
        headers={
            "User-Agent": "phylloscopus-ortholog-pipeline/0.1 (+local manifest builder)",
            "Accept": "application/json, text/plain, */*",
        },
    )


def build_ssl_context(allow_insecure_ssl: bool):
    if allow_insecure_ssl:
        return ssl._create_unverified_context()  # nosec: explicit last-resort flag
    if certifi is not None:
        return ssl.create_default_context(cafile=certifi.where())
    return ssl.create_default_context()


def fetch_bytes(url: str, ssl_context, attempts: int = 3, delay: float = 1.0) -> bytes:
    last_error = None
    for attempt in range(1, attempts + 1):
        try:
            with urlopen(build_request(url), timeout=60, context=ssl_context) as response:
                return response.read()
        except (HTTPError, URLError) as error:
            last_error = error
            if attempt == attempts:
                raise
            time.sleep(delay * attempt)
    raise last_error


def fetch_json(url: str, ssl_context):
    payload = fetch_bytes(url, ssl_context=ssl_context)
    if not payload:
        return None
    return json.loads(payload.decode("utf-8"))


def fetch_text(url: str, ssl_context) -> str:
    payload = fetch_bytes(url, ssl_context=ssl_context)
    return payload.decode("utf-8")


def sleep_if_needed(seconds: float):
    if seconds > 0:
        time.sleep(seconds)


def slugify_species_id(scientific_name: str) -> str:
    value = scientific_name.strip().lower()
    value = re.sub(r"[^a-z0-9]+", "_", value)
    return value.strip("_")


def write_tsv(path: Path, rows: List[Dict[str, str]], fieldnames: List[str]):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def resolve_gbif_genus(genus: str, ssl_context) -> Dict:
    url = f"{GBIF_BASE}/species/match?name={quote(genus)}"
    response = fetch_json(url, ssl_context=ssl_context)
    if not response or response.get("rank") != "GENUS":
        raise SystemExit(f"Could not resolve genus {genus!r} in GBIF.")
    return response


def fetch_gbif_species(genus_key: int, species_limit: int, ssl_context) -> List[Dict]:
    results = []
    offset = 0
    page_size = 100

    while True:
        params = urlencode(
            {
                "rank": "SPECIES",
                "highertaxon_key": genus_key,
                "status": "ACCEPTED",
                "limit": page_size,
                "offset": offset,
            }
        )
        page = fetch_json(f"{GBIF_BASE}/species/search?{params}", ssl_context=ssl_context)
        page_results = page.get("results", []) if page else []
        results.extend(page_results)
        if species_limit and len(results) >= species_limit:
            return results[:species_limit]
        if not page_results or page.get("endOfRecords", True):
            break
        offset += page_size
    return results


def ncbi_esearch(db: str, organism_name: str, api_key: str, ssl_context) -> Dict:
    params = {
        "db": db,
        "term": f"\"{organism_name}\"[Organism]",
        "retmode": "json",
        "retmax": 20,
    }
    if api_key:
        params["api_key"] = api_key
    return fetch_json(f"{NCBI_EUTILS_BASE}/esearch.fcgi?{urlencode(params)}", ssl_context=ssl_context)


def ncbi_assembly_summaries(ids: List[str], api_key: str, ssl_context) -> Dict:
    if not ids:
        return {}
    params = {
        "db": "assembly",
        "id": ",".join(ids),
        "retmode": "json",
    }
    if api_key:
        params["api_key"] = api_key
    return fetch_json(f"{NCBI_EUTILS_BASE}/esummary.fcgi?{urlencode(params)}", ssl_context=ssl_context)


def pick_best_assembly(esummary_payload: Dict) -> Optional[Dict]:
    result = esummary_payload.get("result", {}) if esummary_payload else {}
    uids = result.get("uids", [])
    if not uids:
        return None

    candidates = []
    for uid in uids:
        row = result.get(uid, {})
        assembly_level = row.get("assemblystatus", "")
        level_rank = ASSEMBLY_LEVEL_RANK.get(assembly_level, 99)
        updated = row.get("lastupdatedate", "")
        candidates.append((level_rank, updated, row))

    candidates.sort(key=lambda item: (item[0], item[1]), reverse=False)
    return candidates[0][2] if candidates else None


def ena_taxonomy_lookup(scientific_name: str, ssl_context) -> Optional[Dict]:
    url = f"{ENA_TAXONOMY_BASE}/scientific-name/{quote(scientific_name)}"
    try:
        payload = fetch_json(url, ssl_context=ssl_context)
    except HTTPError as error:
        if error.code == 404:
            return None
        raise
    if not payload:
        return None
    return payload[0]


def ena_read_runs(tax_id: str, ssl_context) -> List[Dict]:
    params = {
        "result": "read_run",
        "fields": ",".join(ENA_READ_RUN_FIELDS),
        "format": "json",
        "limit": 0,
        "query": f"tax_eq({tax_id})",
    }
    payload = fetch_json(f"{ENA_PORTAL_BASE}/search?{urlencode(params)}", ssl_context=ssl_context)
    return payload or []


def is_transcriptomic_run(run: Dict) -> bool:
    library_source = (run.get("library_source") or "").strip().upper()
    library_strategy = (run.get("library_strategy") or "").strip()
    return library_source in TRANSCRIPTOMIC_SOURCES or library_strategy in TRANSCRIPTOMIC_STRATEGIES


def is_genomic_wgs_run(run: Dict) -> bool:
    library_strategy = (run.get("library_strategy") or "").strip()
    return library_strategy in GENOMIC_WGS_STRATEGIES


def summarize_species(
    species: Dict,
    args,
    ssl_context,
    raw_ncbi_dir: Path,
    raw_ena_dir: Path,
) -> Tuple[Dict[str, str], List[Dict[str, str]]]:
    canonical_name = species.get("canonicalName", "").strip()
    scientific_name = canonical_name or species.get("species", "").strip()
    species_id = slugify_species_id(scientific_name)

    assembly_search = ncbi_esearch("assembly", scientific_name, args.ncbi_api_key, ssl_context)
    sleep_if_needed(args.request_sleep_seconds)
    bioproject_search = ncbi_esearch("bioproject", scientific_name, args.ncbi_api_key, ssl_context)
    sleep_if_needed(args.request_sleep_seconds)
    sra_search = ncbi_esearch("sra", scientific_name, args.ncbi_api_key, ssl_context)
    sleep_if_needed(args.request_sleep_seconds)

    assembly_ids = assembly_search.get("esearchresult", {}).get("idlist", [])
    assembly_count = int(assembly_search.get("esearchresult", {}).get("count", "0"))
    bioproject_count = int(bioproject_search.get("esearchresult", {}).get("count", "0"))
    sra_count = int(sra_search.get("esearchresult", {}).get("count", "0"))

    assembly_summary = (
        ncbi_assembly_summaries(assembly_ids, args.ncbi_api_key, ssl_context) if assembly_ids else {}
    )
    if assembly_ids:
        sleep_if_needed(args.request_sleep_seconds)
    best_assembly = pick_best_assembly(assembly_summary)

    ncbi_raw_path = raw_ncbi_dir / f"{species_id}.json"
    ncbi_raw_path.write_text(
        json.dumps(
            {
                "assembly_search": assembly_search,
                "assembly_summary": assembly_summary,
                "bioproject_search": bioproject_search,
                "sra_search": sra_search,
            },
            indent=2,
        )
        + "\n"
    )

    ena_tax = None
    ena_runs = []
    ena_read_run_count = 0
    ena_transcriptomic_runs: List[str] = []
    ena_genomic_runs: List[str] = []
    run_metadata_rows: List[Dict[str, str]] = []

    if args.include_ena:
        ena_tax = ena_taxonomy_lookup(scientific_name, ssl_context)
        sleep_if_needed(args.request_sleep_seconds)
        if ena_tax:
            ena_runs = ena_read_runs(ena_tax["taxId"], ssl_context)
            sleep_if_needed(args.request_sleep_seconds)
            ena_read_run_count = len(ena_runs)
            for run in ena_runs:
                accession = (run.get("run_accession") or "").strip()
                if not accession:
                    continue
                evidence_class = "other"
                if is_transcriptomic_run(run):
                    ena_transcriptomic_runs.append(accession)
                    evidence_class = "rna"
                if is_genomic_wgs_run(run):
                    ena_genomic_runs.append(accession)
                    evidence_class = "wgs"
                run_metadata_rows.append(
                    {
                        "species_id": species_id,
                        "scientific_name": scientific_name,
                        "ena_tax_id": (ena_tax or {}).get("taxId", ""),
                        "run_accession": accession,
                        "evidence_class": evidence_class,
                        "library_layout": (run.get("library_layout") or "").strip(),
                        "library_source": (run.get("library_source") or "").strip(),
                        "library_strategy": (run.get("library_strategy") or "").strip(),
                        "library_selection": (run.get("library_selection") or "").strip(),
                        "study_accession": (run.get("study_accession") or "").strip(),
                        "experiment_accession": (run.get("experiment_accession") or "").strip(),
                        "sample_accession": (run.get("sample_accession") or "").strip(),
                        "instrument_platform": (run.get("instrument_platform") or "").strip(),
                        "instrument_model": (run.get("instrument_model") or "").strip(),
                        "read_count": (run.get("read_count") or "").strip(),
                        "base_count": (run.get("base_count") or "").strip(),
                        "fastq_ftp": (run.get("fastq_ftp") or "").strip(),
                        "fastq_bytes": (run.get("fastq_bytes") or "").strip(),
                    }
                )

        ena_raw_path = raw_ena_dir / f"{species_id}.json"
        ena_raw_path.write_text(
            json.dumps(
                {
                    "taxonomy": ena_tax,
                    "read_runs": ena_runs,
                },
                indent=2,
            )
            + "\n"
        )

    assembly_accession = (best_assembly or {}).get("assemblyaccession", "")
    assembly_level = (best_assembly or {}).get("assemblystatus", "")
    assembly_name = (best_assembly or {}).get("assemblyname", "")

    if assembly_accession:
        evidence_hint = "B"
    elif ena_transcriptomic_runs:
        evidence_hint = "C"
    elif ena_genomic_runs or sra_count > 0:
        evidence_hint = "D"
    else:
        evidence_hint = "E"

    notes = (
        f"GBIF accepted species; ncbi_assemblies={assembly_count}; "
        f"ncbi_bioprojects={bioproject_count}; ncbi_sra={sra_count}; "
        f"ena_runs={ena_read_run_count}; best_assembly_level={assembly_level or 'none'}"
    )

    return {
        "species_id": species_id,
        "scientific_name": scientific_name,
        "gbif_species_key": str(species.get("speciesKey") or species.get("key") or ""),
        "gbif_scientific_name": species.get("scientificName", "").strip(),
        "gbif_authorship": (species.get("authorship") or "").strip(),
        "ena_tax_id": (ena_tax or {}).get("taxId", ""),
        "assembly_accession": assembly_accession,
        "assembly_level": assembly_level,
        "assembly_name": assembly_name,
        "ncbi_assembly_count": str(assembly_count),
        "ncbi_bioproject_count": str(bioproject_count),
        "ncbi_sra_count": str(sra_count),
        "ena_read_run_count": str(ena_read_run_count),
        "ena_transcriptomic_run_count": str(len(ena_transcriptomic_runs)),
        "ena_genomic_run_count": str(len(ena_genomic_runs)),
        "assembly_fasta": "",
        "annotation_gtf": "",
        "transcriptome_fasta": "",
        "rna_sra_accessions": ",".join(ena_transcriptomic_runs[: args.max_runs_per_class]),
        "wgs_sra_accessions": ",".join(ena_genomic_runs[: args.max_runs_per_class]),
        "evidence_hint": evidence_hint,
        "notes": notes,
    }, run_metadata_rows


def build_taxonomy_source(
    args,
    genus_record: Dict,
    species_records: List[Dict],
    manifest_path: Path,
) -> Dict:
    return {
        "snapshot_date": datetime.now(timezone.utc).date().isoformat(),
        "snapshot_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "authority": args.authority,
        "query_genus": args.genus,
        "species_limit": args.species_limit,
        "resolved_genus": {
            "gbif_usage_key": genus_record.get("usageKey"),
            "scientific_name": genus_record.get("scientificName"),
            "canonical_name": genus_record.get("canonicalName"),
            "family": genus_record.get("family"),
        },
        "query_urls": {
            "gbif_genus_match": f"{GBIF_BASE}/species/match?name={quote(args.genus)}",
            "gbif_species_search": (
                f"{GBIF_BASE}/species/search?rank=SPECIES&highertaxon_key={genus_record.get('usageKey')}"
                "&status=ACCEPTED"
            ),
            "ncbi_eutils_base": NCBI_EUTILS_BASE,
            "ena_portal_base": ENA_PORTAL_BASE,
            "ena_taxonomy_base": ENA_TAXONOMY_BASE,
        },
        "species_count": len(species_records),
        "manifest_path": str(manifest_path),
        "include_ena": args.include_ena,
    }


def build_discovery_summary(rows: List[Dict[str, str]]) -> List[Dict[str, str]]:
    tier_counts = {"B": 0, "C": 0, "D": 0, "E": 0}
    for row in rows:
        tier_counts[row["evidence_hint"]] = tier_counts.get(row["evidence_hint"], 0) + 1

    return [
        {"metric": "species_count", "value": str(len(rows))},
        {"metric": "tier_B", "value": str(tier_counts.get("B", 0))},
        {"metric": "tier_C", "value": str(tier_counts.get("C", 0))},
        {"metric": "tier_D", "value": str(tier_counts.get("D", 0))},
        {"metric": "tier_E", "value": str(tier_counts.get("E", 0))},
        {
            "metric": "species_with_assembly_accession",
            "value": str(sum(1 for row in rows if row["assembly_accession"])),
        },
        {
            "metric": "species_with_rna_runs",
            "value": str(sum(1 for row in rows if row["rna_sra_accessions"])),
        },
        {
            "metric": "species_with_wgs_runs",
            "value": str(sum(1 for row in rows if row["wgs_sra_accessions"])),
        },
    ]


def main():
    args = parse_args()
    ssl_context = build_ssl_context(args.allow_insecure_ssl)
    outdir = Path(args.outdir)
    raw_dir = outdir / "raw"
    raw_ncbi_dir = raw_dir / "ncbi"
    raw_ena_dir = raw_dir / "ena"
    outdir.mkdir(parents=True, exist_ok=True)
    raw_ncbi_dir.mkdir(parents=True, exist_ok=True)
    if args.include_ena:
        raw_ena_dir.mkdir(parents=True, exist_ok=True)

    genus_record = resolve_gbif_genus(args.genus, ssl_context)
    sleep_if_needed(args.request_sleep_seconds)
    species_records = fetch_gbif_species(genus_record["usageKey"], args.species_limit, ssl_context)
    (raw_dir / "gbif_species.json").write_text(json.dumps(species_records, indent=2) + "\n")

    manifest_rows = []
    run_metadata_rows = []
    for species in species_records:
        manifest_row, species_run_metadata = summarize_species(species, args, ssl_context, raw_ncbi_dir, raw_ena_dir)
        manifest_rows.append(manifest_row)
        run_metadata_rows.extend(species_run_metadata)

    manifest_rows.sort(key=lambda row: row["scientific_name"])

    manifest_path = outdir / "species_manifest.tsv"
    write_tsv(
        manifest_path,
        manifest_rows,
        [
            "species_id",
            "scientific_name",
            "gbif_species_key",
            "gbif_scientific_name",
            "gbif_authorship",
            "ena_tax_id",
            "assembly_accession",
            "assembly_level",
            "assembly_name",
            "ncbi_assembly_count",
            "ncbi_bioproject_count",
            "ncbi_sra_count",
            "ena_read_run_count",
            "ena_transcriptomic_run_count",
            "ena_genomic_run_count",
            "assembly_fasta",
            "annotation_gtf",
            "transcriptome_fasta",
            "rna_sra_accessions",
            "wgs_sra_accessions",
            "evidence_hint",
            "notes",
        ],
    )

    snapshot_date = datetime.now(timezone.utc).date().isoformat()
    (outdir / "snapshot_date").write_text(snapshot_date + "\n")

    taxonomy_source = build_taxonomy_source(args, genus_record, species_records, manifest_path)
    (outdir / "taxonomy_source.json").write_text(json.dumps(taxonomy_source, indent=2) + "\n")

    summary_rows = build_discovery_summary(manifest_rows)
    write_tsv(outdir / "discovery_summary.tsv", summary_rows, ["metric", "value"])

    if args.include_ena:
        write_tsv(
            outdir / "run_metadata.tsv",
            sorted(run_metadata_rows, key=lambda row: (row["species_id"], row["run_accession"])),
            [
                "species_id",
                "scientific_name",
                "ena_tax_id",
                "run_accession",
                "evidence_class",
                "library_layout",
                "library_source",
                "library_strategy",
                "library_selection",
                "study_accession",
                "experiment_accession",
                "sample_accession",
                "instrument_platform",
                "instrument_model",
                "read_count",
                "base_count",
                "fastq_ftp",
                "fastq_bytes",
            ],
        )


if __name__ == "__main__":
    main()
