#!/usr/bin/env python3

import argparse
import csv
import importlib.util
from pathlib import Path
from typing import Dict, List, Optional, Tuple


def parse_args():
    parser = argparse.ArgumentParser(
        description="Write narrow Tier B sanity tables into an existing validation archive."
    )
    parser.add_argument("--archive-dir", required=True)
    parser.add_argument("--assembly-fasta", required=True)
    return parser.parse_args()


def load_tier_ab_module():
    module_path = Path(__file__).resolve().parent / "run_tier_ab_scaffold.py"
    spec = importlib.util.spec_from_file_location("run_tier_ab_scaffold", module_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def read_tsv(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: List[Dict[str, str]], fieldnames: List[str]):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def read_fasta(path: Path) -> Dict[str, str]:
    records: Dict[str, str] = {}
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
                header = line[1:]
                chunks = []
            else:
                chunks.append(line)
    if header is not None:
        records[header] = "".join(chunks).upper()
    return records


def ambiguous_base_count(sequence: str) -> int:
    return sum(1 for base in sequence.upper() if base not in {"A", "C", "G", "T"})


def smith_waterman(query: str, target: str) -> Dict[str, float]:
    if not query or not target:
        return {
            "aligned_span": 0,
            "identity": 0.0,
            "query_start_1based": 0,
            "query_end_1based": 0,
            "target_start_1based": 0,
            "target_end_1based": 0,
        }

    match_score = 2
    mismatch_penalty = -3
    gap_penalty = -4
    rows = len(query) + 1
    cols = len(target) + 1
    scores = [[0] * cols for _ in range(rows)]
    traceback: List[List[Optional[Tuple[int, int]]]] = [[None] * cols for _ in range(rows)]
    best_score = 0
    best_position = (0, 0)

    for i in range(1, rows):
        query_base = query[i - 1]
        for j in range(1, cols):
            target_base = target[j - 1]
            diag_score = scores[i - 1][j - 1] + (match_score if query_base == target_base else mismatch_penalty)
            up_score = scores[i - 1][j] + gap_penalty
            left_score = scores[i][j - 1] + gap_penalty
            score = max(0, diag_score, up_score, left_score)
            scores[i][j] = score
            if score == 0:
                traceback[i][j] = None
            elif score == diag_score:
                traceback[i][j] = (i - 1, j - 1)
            elif score == up_score:
                traceback[i][j] = (i - 1, j)
            else:
                traceback[i][j] = (i, j - 1)
            if score > best_score:
                best_score = score
                best_position = (i, j)

    aligned_query: List[str] = []
    aligned_target: List[str] = []
    i, j = best_position
    end_i = i
    end_j = j
    while i > 0 and j > 0 and scores[i][j] > 0:
        prev = traceback[i][j]
        if prev is None:
            break
        prev_i, prev_j = prev
        if prev_i == i - 1 and prev_j == j - 1:
            aligned_query.append(query[i - 1])
            aligned_target.append(target[j - 1])
        elif prev_i == i - 1 and prev_j == j:
            aligned_query.append(query[i - 1])
            aligned_target.append("-")
        else:
            aligned_query.append("-")
            aligned_target.append(target[j - 1])
        i, j = prev_i, prev_j

    aligned_query.reverse()
    aligned_target.reverse()
    aligned_span = sum(1 for query_base, target_base in zip(aligned_query, aligned_target) if query_base != "-" and target_base != "-")
    matches = sum(1 for query_base, target_base in zip(aligned_query, aligned_target) if query_base == target_base and query_base != "-")
    identity = matches / aligned_span if aligned_span else 0.0
    return {
        "aligned_span": aligned_span,
        "identity": identity,
        "query_start_1based": i + 1 if aligned_span else 0,
        "query_end_1based": end_i if aligned_span else 0,
        "target_start_1based": j + 1 if aligned_span else 0,
        "target_end_1based": end_j if aligned_span else 0,
    }


def separation_call(best_hit: Dict[str, object], second_hit: Optional[Dict[str, object]]) -> Tuple[str, str, str, str]:
    if second_hit is None:
        return ("yes", "no second hit passed the reporting threshold", "NA", "NA")

    score_delta = int(best_hit["score"]) - int(second_hit["score"])
    identity_delta = float(best_hit["identity"]) - float(second_hit["identity"])
    coverage_delta = float(best_hit["coverage"]) - float(second_hit["coverage"])
    clearly_separated = "yes" if (
        score_delta >= 100 and (identity_delta >= 0.05 or coverage_delta >= 0.05)
    ) else "no"
    explanation = (
        f"best-vs-second score_delta={score_delta}; "
        f"identity_delta={identity_delta:.3f}; coverage_delta={coverage_delta:.3f}"
    )
    return (clearly_separated, explanation, str(score_delta), f"{identity_delta:.3f}")


def main():
    args = parse_args()
    archive_dir = Path(args.archive_dir).resolve()
    summary_dir = archive_dir / "summary_tables"
    assembly_fasta = Path(args.assembly_fasta).resolve()

    tier_ab = load_tier_ab_module()
    status_rows = read_tsv(archive_dir / "ortholog_status_long.tsv")
    target_rows = read_tsv(archive_dir / "manifests_used" / "ortholog_targets.tsv")
    target_lookup = {row["gene_id"].strip(): row for row in target_rows}
    fasta_records = read_fasta(archive_dir / "ortholog_sequences" / "RAG1.fna")

    tier_b_row = next(
        row
        for row in status_rows
        if row["species_id"].strip() == "phylloscopus_collybita" and row["gene_id"].strip() == "RAG1"
    )
    tier_a_row = next(
        row
        for row in status_rows
        if row["species_id"].strip() == "phylloscopus_trochilus" and row["gene_id"].strip() == "RAG1"
    )
    target_row = target_lookup["RAG1"]
    query_path = Path(tier_b_row["query_source"]).resolve()
    query_records = tier_ab.read_fasta_dict(query_path)
    query_name, query_sequence = next(iter(query_records.items()))
    selected_hit = tier_ab.recover_candidate_sequence(query_path, assembly_fasta)
    if selected_hit is None:
        raise SystemExit("No Tier B candidate hit could be recovered from the supplied assembly.")
    ranked_hits = tier_ab.rank_candidate_hits(query_path, assembly_fasta)
    selected_key = (
        str(selected_hit["contig_id"]),
        str(selected_hit["strand"]),
        int(selected_hit["start_1based"]),
        int(selected_hit["end_1based"]),
    )
    ordered_hits: List[Dict[str, object]] = []
    selected_rank = "not_ranked"
    for hit in ranked_hits:
        hit_key = (
            str(hit["contig_id"]),
            str(hit["strand"]),
            int(hit["start_1based"]),
            int(hit["end_1based"]),
        )
        if hit_key == selected_key:
            selected_rank = str(hit["hit_rank"])
            selected_hit["hit_rank"] = hit["hit_rank"]
            selected_hit["aligned_length"] = hit["aligned_length"]
            selected_hit["recovered_length"] = hit["recovered_length"]
            selected_hit["score"] = hit["score"]
            selected_hit["identity"] = hit["identity"]
            selected_hit["coverage"] = hit["coverage"]
        else:
            ordered_hits.append(hit)
    reported_hits = [selected_hit, *ordered_hits[:4]]
    second_hit = ordered_hits[0] if ordered_hits else None
    clearly_separated, explanation, score_delta, identity_delta = separation_call(selected_hit, second_hit)

    candidate_summary_rows = [
        {
            "gene_id": "RAG1",
            "species_id": "phylloscopus_collybita",
            "label": "candidate assembly-backed locus sequence",
            "assembly_backed": tier_b_row["assembly_backed"],
            "orthology_unvalidated": tier_b_row["orthology_unvalidated"],
            "contig": selected_hit["contig_id"],
            "start_1based": str(selected_hit["start_1based"]),
            "end_1based": str(selected_hit["end_1based"]),
            "strand": str(selected_hit["strand"]),
            "sequence_length": tier_b_row["sequence_length"],
            "query_length_used_for_retrieval": str(len(query_sequence)),
            "nominal_target_cds_length": target_row["cds_length"].strip(),
            "query_coverage": f"{float(selected_hit['coverage']):.3f}",
            "match_identity": f"{float(selected_hit['identity']):.3f}",
            "interpretation_note": (
                "assembly hit spans the full search query length used for retrieval; "
                "the recovered region is still shorter than the nominal target CDS length; "
                "orthology and structure remain unvalidated"
            ),
        }
    ]
    write_tsv(
        summary_dir / "tier_b_candidate_summary.tsv",
        candidate_summary_rows,
        [
            "gene_id",
            "species_id",
            "label",
            "assembly_backed",
            "orthology_unvalidated",
            "contig",
            "start_1based",
            "end_1based",
            "strand",
            "sequence_length",
            "query_length_used_for_retrieval",
            "nominal_target_cds_length",
            "query_coverage",
            "match_identity",
            "interpretation_note",
        ],
    )

    second_best_rows = []
    for hit in reported_hits:
        hit_rank_label = "selected_best_hit" if hit is selected_hit else str(hit["hit_rank"])
        second_best_rows.append(
            {
                "gene_id": "RAG1",
                "species_id": "phylloscopus_collybita",
                "query_record": query_name,
                "query_length": str(len(query_sequence)),
                "hit_rank": hit_rank_label,
                "contig": str(hit["contig_id"]),
                "start_1based": str(hit["start_1based"]),
                "end_1based": str(hit["end_1based"]),
                "strand": str(hit["strand"]),
                "aligned_length": str(hit["aligned_length"]),
                "recovered_length": str(hit["recovered_length"]),
                "percent_identity": f"{float(hit['identity']) * 100:.2f}",
                "query_coverage": f"{float(hit['coverage']) * 100:.2f}",
                "score": str(hit["score"]),
                "selected_hit_rank_among_all_hits": selected_rank,
                "score_delta_vs_best": "0" if hit is selected_hit else str(int(selected_hit["score"]) - int(hit["score"])),
                "identity_delta_vs_best": (
                    "0.000"
                    if hit is selected_hit
                    else f"{float(selected_hit['identity']) - float(hit['identity']):.3f}"
                ),
                "coverage_delta_vs_best": (
                    "0.000"
                    if hit is selected_hit
                    else f"{float(selected_hit['coverage']) - float(hit['coverage']):.3f}"
                ),
                "best_hit_clearly_separated_from_second": clearly_separated,
                "separation_note": explanation,
            }
        )
    write_tsv(
        summary_dir / "second_best_hit_sanity.tsv",
        second_best_rows,
        [
            "gene_id",
            "species_id",
            "query_record",
            "query_length",
            "hit_rank",
            "contig",
            "start_1based",
            "end_1based",
            "strand",
            "aligned_length",
            "recovered_length",
            "percent_identity",
            "query_coverage",
            "score",
            "selected_hit_rank_among_all_hits",
            "score_delta_vs_best",
            "identity_delta_vs_best",
            "coverage_delta_vs_best",
            "best_hit_clearly_separated_from_second",
            "separation_note",
        ],
    )

    tier_a_header = next(header for header in fasta_records if header.startswith("phylloscopus_trochilus|RAG1"))
    tier_b_header = next(header for header in fasta_records if header.startswith("phylloscopus_collybita|RAG1"))
    tier_a_sequence = fasta_records[tier_a_header]
    tier_b_sequence = fasta_records[tier_b_header]
    pairwise = smith_waterman(tier_a_sequence, tier_b_sequence)
    pairwise_rows = [
        {
            "gene_id": "RAG1",
            "species_id": "phylloscopus_trochilus",
            "evidence_tier": tier_a_row["evidence_tier"],
            "label_or_status": tier_a_row["reconstruction_status"],
            "sequence_length": str(len(tier_a_sequence)),
            "ambiguous_base_count": str(ambiguous_base_count(tier_a_sequence)),
            "pairwise_partner_species": "phylloscopus_collybita",
            "pairwise_aligned_span": str(pairwise["aligned_span"]),
            "pairwise_identity": f"{pairwise['identity']:.3f}",
            "pairwise_query_start_1based": str(pairwise["query_start_1based"]),
            "pairwise_query_end_1based": str(pairwise["query_end_1based"]),
            "pairwise_target_start_1based": str(pairwise["target_start_1based"]),
            "pairwise_target_end_1based": str(pairwise["target_end_1based"]),
            "note": "Tier A sequence is the bundled mock reference transcript retained only to anchor this proof-of-feasibility archive.",
        },
        {
            "gene_id": "RAG1",
            "species_id": "phylloscopus_collybita",
            "evidence_tier": tier_b_row["evidence_tier"],
            "label_or_status": tier_b_row["reconstruction_status"],
            "sequence_length": str(len(tier_b_sequence)),
            "ambiguous_base_count": str(ambiguous_base_count(tier_b_sequence)),
            "pairwise_partner_species": "phylloscopus_trochilus",
            "pairwise_aligned_span": str(pairwise["aligned_span"]),
            "pairwise_identity": f"{pairwise['identity']:.3f}",
            "pairwise_query_start_1based": str(pairwise["target_start_1based"]),
            "pairwise_query_end_1based": str(pairwise["target_end_1based"]),
            "pairwise_target_start_1based": str(pairwise["query_start_1based"]),
            "pairwise_target_end_1based": str(pairwise["query_end_1based"]),
            "note": "Pairwise sanity only; candidate sequence remains assembly-backed and does not validate orthology, exon structure, or completeness.",
        },
    ]
    pairwise_fieldnames = [
        "gene_id",
        "species_id",
        "evidence_tier",
        "label_or_status",
        "sequence_length",
        "ambiguous_base_count",
        "pairwise_partner_species",
        "pairwise_aligned_span",
        "pairwise_identity",
        "pairwise_query_start_1based",
        "pairwise_query_end_1based",
        "pairwise_target_start_1based",
        "pairwise_target_end_1based",
        "note",
    ]
    write_tsv(summary_dir / "pairwise_tierA_tierB_RAG1_sanity.tsv", pairwise_rows, pairwise_fieldnames)
    write_tsv(summary_dir / "sequence_sanity_check.tsv", pairwise_rows, pairwise_fieldnames)


if __name__ == "__main__":
    main()
