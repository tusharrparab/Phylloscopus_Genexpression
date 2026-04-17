# Expression Scaffold

The expression scaffold is implemented in [run_expression_scaffold.py](/Users/sam/Documents/New%20project/bin/run_expression_scaffold.py).

It converts RNA-backed species in the species manifest into a sample sheet and quantification plan.

## What It Does

1. extracts RNA run accessions from `rna_sra_accessions`
2. chooses the primary reference from `reference_manifest.tsv`
3. uses a provided transcript FASTA when available
4. otherwise tries to build transcripts from genome + GTF using `gffread -w ... -g ...` ([gffread docs](https://ccb.jhu.edu/software/stringtie/gff.shtml))
5. plans transcript quantification with `salmon index` and `salmon quant` ([Salmon docs](https://salmon.readthedocs.io/en/latest/salmon.html))
6. writes a design template for downstream DESeq2 or edgeR analysis

For SRA download planning, the scaffold uses `prefetch` + `fasterq-dump`, which the SRA Toolkit documents as the standard fast path for FASTQ extraction ([SRA Toolkit wiki](https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump)).

## Typical Usage

```bash
python3 bin/run_expression_scaffold.py \
  --species-manifest examples/species_manifest.tsv \
  --reference-manifest examples/reference_manifest.tsv \
  --outdir results/expression \
  --mode scaffold
```

## Outputs

- `expression/expression_samples.tsv`
- `expression/expression_design_template.tsv`
- `expression/expression_plan.tsv`
- `expression/expression_reference.tsv`
- `expression/expression_summary.tsv`

If the scaffold cannot determine library layout, it does not emit a misleading paired-end command. Instead it writes commented single-end and paired-end `salmon quant` examples and marks the sample as `library_layout_unknown` until metadata are confirmed.

If the primary reference lacks `transcript_fasta` and cannot build one from `assembly_fasta` + `annotation_gtf`, the scaffold records `transcript_fasta_missing` and leaves quantification in planning mode.
