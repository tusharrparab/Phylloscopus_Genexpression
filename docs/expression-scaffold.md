# Expression Scaffold

The expression scaffold is implemented in [run_expression_scaffold.py](/Users/sam/Documents/New%20project/bin/run_expression_scaffold.py).

It converts RNA-backed species in the species manifest into a sample sheet and quantification plan.

## What It Does

1. extracts RNA run accessions from `rna_sra_accessions`
2. chooses the primary reference from `reference_manifest.tsv`
3. uses a provided transcript FASTA when available
4. otherwise tries to build transcripts from genome + GTF using `gffread -w ... -g ...` ([gffread docs](https://ccb.jhu.edu/software/stringtie/gff.shtml))
5. resolves run-level metadata from `run_metadata.tsv`, ENA, and then NCBI SRA RunInfo so single-end versus paired-end runs do not stay `library_layout_unknown` when public metadata are available
6. plans or runs transcript quantification with `salmon index` and `salmon quant` ([Salmon docs](https://salmon.readthedocs.io/en/latest/salmon.html))
7. writes a design template for downstream DESeq2 or edgeR analysis

For SRA download planning, the scaffold uses `prefetch` + `fasterq-dump`, which the SRA Toolkit documents as the standard fast path for FASTQ extraction ([SRA Toolkit wiki](https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump)).

## Typical Usage

```bash
python3 bin/run_expression_scaffold.py \
  --species-manifest examples/species_manifest.tsv \
  --reference-manifest examples/reference_manifest.tsv \
  --outdir results/expression \
  --mode scaffold
```

To execute quantification instead of only writing runnable plans:

```bash
python3 bin/run_expression_scaffold.py \
  --species-manifest snapshots/phylloscopus_2026-04-17/species_manifest.tsv \
  --reference-manifest snapshots/phylloscopus_2026-04-17/reference_manifest.tsv \
  --run-metadata snapshots/phylloscopus_2026-04-17/run_metadata.tsv \
  --outdir results/expression \
  --mode execute
```

## Outputs

- `expression/expression_samples.tsv`
- `expression/expression_design_template.tsv`
- `expression/expression_plan.tsv`
- `expression/expression_reference.tsv`
- `expression/expression_summary.tsv`

If a local `run_metadata.tsv` row is incomplete, the scaffold now attempts live ENA lookup by `run_accession` and then falls back to NCBI SRA RunInfo before leaving the sample unresolved.

If the primary reference lacks `transcript_fasta` and cannot build one from `assembly_fasta` + `annotation_gtf`, the scaffold records `transcript_fasta_missing` and leaves quantification in planning mode.
