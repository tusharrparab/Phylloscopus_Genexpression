#!/usr/bin/env python3

import argparse
import csv
import shlex
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple


def parse_args():
    parser = argparse.ArgumentParser(
        description="Scaffold ancestral sequence reconstruction from ortholog FASTA files."
    )
    parser.add_argument("--sequence-dir", required=True, help="Directory containing per-locus FASTA files.")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--species-tree", default="", help="Optional user-supplied Newick tree.")
    parser.add_argument("--mode", choices=["plan", "scaffold"], default="scaffold")
    parser.add_argument("--min-taxa", type=int, default=4)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--sequence-type", choices=["DNA", "AA", "CODON"], default="DNA")
    return parser.parse_args()


def read_fasta(path: Path) -> List[Tuple[str, str]]:
    records = []
    header = None
    chunks: List[str] = []
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks)))
                header = line[1:]
                chunks = []
            else:
                chunks.append(line)
    if header is not None:
        records.append((header, "".join(chunks)))
    return records


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


def main():
    args = parse_args()
    sequence_dir = Path(args.sequence_dir)
    outdir = Path(args.outdir)
    align_dir = outdir / "alignments"
    plan_dir = outdir / "plans"
    tree_dir = outdir / "trees"
    outdir.mkdir(parents=True, exist_ok=True)
    align_dir.mkdir(parents=True, exist_ok=True)
    plan_dir.mkdir(parents=True, exist_ok=True)
    tree_dir.mkdir(parents=True, exist_ok=True)

    mafft_exe = shutil.which("mafft")
    iqtree_exe = shutil.which("iqtree2") or shutil.which("iqtree")

    locus_rows = []
    asr_rows = []
    summary_rows = []

    fasta_paths = sorted(
        [
            path
            for pattern in ("*.fna", "*.fa", "*.fasta")
            for path in sequence_dir.glob(pattern)
        ]
    )

    for fasta_path in fasta_paths:
        locus_id = fasta_path.stem
        records = read_fasta(fasta_path)
        taxa = [header.split("|")[0] for header, _ in records]
        taxa_count = len(taxa)
        site_lengths = [len(sequence) for _, sequence in records]
        min_len = min(site_lengths) if site_lengths else 0
        max_len = max(site_lengths) if site_lengths else 0
        locus_rows.append(
            {
                "locus_id": locus_id,
                "input_fasta": str(fasta_path.resolve()),
                "taxa_count": str(taxa_count),
                "min_length": str(min_len),
                "max_length": str(max_len),
            }
        )

        locus_plan_dir = plan_dir / locus_id
        locus_plan_dir.mkdir(parents=True, exist_ok=True)
        alignment_path = align_dir / f"{locus_id}.aligned.fasta"
        tree_prefix = tree_dir / locus_id
        plan_script = locus_plan_dir / "run_asr.sh"
        plan_log = locus_plan_dir / "asr.log"
        state_path = tree_dir / f"{locus_id}.state"
        treefile_path = tree_dir / f"{locus_id}.treefile"

        if taxa_count < args.min_taxa:
            status = "too_few_taxa"
            notes = f"requires at least {args.min_taxa} taxa for ASR planning"
        else:
            status = "planned"
            notes = []
            analysis_input = fasta_path.resolve()
            if args.mode == "scaffold" and mafft_exe:
                code = run_command([mafft_exe, "--auto", str(fasta_path.resolve())], cwd=locus_plan_dir, log_path=plan_log)
                if code == 0 and plan_log.exists():
                    alignment_path.write_text(plan_log.read_text())
                    analysis_input = alignment_path.resolve()
                    notes.append("alignment_completed_with_mafft")
                else:
                    notes.append("alignment_failed_fallback_to_input_fasta")
            else:
                notes.append("mafft_missing" if not mafft_exe else "alignment_planned")

            cmd = [iqtree_exe or "iqtree2", "-s", str(analysis_input)]
            cmd.extend(["-st", args.sequence_type, "-m", "MFP", "-asr", "-T", "AUTO"])
            if args.species_tree:
                cmd.extend(["-te", args.species_tree])

            plan_script.write_text(
                "#!/usr/bin/env bash\nset -euo pipefail\n" + " ".join(shlex.quote(part) for part in cmd) + "\n"
            )
            plan_script.chmod(0o755)

            if args.mode == "scaffold" and iqtree_exe:
                code = run_command(cmd, cwd=tree_dir, log_path=plan_log)
                if code == 0:
                    notes.append("iqtree_asr_completed")
                    status = "completed"
                else:
                    notes.append("iqtree_asr_failed")
                    status = "iqtree_failed"
            else:
                notes.append("iqtree_missing" if not iqtree_exe else "iqtree_planned")

            notes = "; ".join(notes)

        asr_rows.append(
            {
                "locus_id": locus_id,
                "status": status,
                "taxa_count": str(taxa_count),
                "alignment_path": str(alignment_path) if alignment_path.exists() else "",
                "treefile_path": str(treefile_path) if treefile_path.exists() else "",
                "state_path": str(state_path) if state_path.exists() else "",
                "plan_script": str(plan_script) if plan_script.exists() else "",
                "notes": notes,
            }
        )

    summary_rows.append({"metric": "locus_count", "value": str(len(locus_rows))})
    summary_rows.append(
        {
            "metric": "loci_with_enough_taxa",
            "value": str(sum(1 for row in asr_rows if row["status"] != "too_few_taxa")),
        }
    )
    summary_rows.append(
        {
            "metric": "completed_asr_runs",
            "value": str(sum(1 for row in asr_rows if row["status"] == "completed")),
        }
    )

    write_tsv(outdir / "locus_manifest.tsv", locus_rows, ["locus_id", "input_fasta", "taxa_count", "min_length", "max_length"])
    write_tsv(
        outdir / "asr_plan.tsv",
        asr_rows,
        ["locus_id", "status", "taxa_count", "alignment_path", "treefile_path", "state_path", "plan_script", "notes"],
    )
    write_tsv(outdir / "asr_summary.tsv", summary_rows, ["metric", "value"])


if __name__ == "__main__":
    main()
