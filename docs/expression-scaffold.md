# Expression Scaffold

The expression module exists to support RNA-backed taxa and optional downstream analyses. It is not the repository's central biological product.

## Implemented Now

- selection of RNA-backed samples from the species manifest
- reference transcript FASTA resolution or construction with `gffread`
- run-level metadata enrichment from local metadata, ENA, and NCBI SRA RunInfo
- `salmon` index construction
- `salmon quant` plan or execution
- sample-sheet and design-template generation

## Still Missing

- ortholog-level count aggregation across species
- transcript-to-gene harmonization across references
- batch/tissue/season design enforcement
- downstream statistical testing
- phylogenetically aware comparative expression models

## Use It For

- preparing RNA-backed species for later ortholog-aware analysis
- staging transcript evidence for Tier C taxa
- local quantification experiments within a single reference context

## Do Not Use It For Yet

- generic cross-species expression claims
- automatic differential expression across taxa
- treating this repository as primarily an RNA-seq project

## Example

```bash
python3 bin/run_expression_scaffold.py \
  --species-manifest snapshots/phylloscopus_2026-04-17/species_manifest.tsv \
  --reference-manifest snapshots/phylloscopus_2026-04-17/reference_manifest.tsv \
  --run-metadata snapshots/phylloscopus_2026-04-17/run_metadata.tsv \
  --outdir results/expression \
  --mode execute
```
