# Expression Scope

## Status

Expression is optional and secondary in this repository.

The codebase can already stage or execute:

- run-level metadata enrichment from local metadata, ENA, and NCBI SRA RunInfo
- transcript FASTA generation with `gffread`
- `salmon` index construction
- `salmon quant` planning or execution

## What This Does Not Mean

These capabilities do not make the repository a comparative expression pipeline.

Missing pieces include:

- ortholog-level count aggregation across species
- tissue harmonization
- seasonal and developmental covariate handling
- batch correction strategy
- transcript-to-gene comparability across references
- phylogenetically aware comparative expression models

## Defensible Use Right Now

Right now, expression should be used only as:

- optional support for RNA-backed Tier C taxa
- a way to stage transcript evidence and sample metadata
- a path toward later ortholog-aware transcript quantification

## Not Yet Defensible

- direct cross-species differential expression claims
- treating `salmon` output from different references as inherently comparable
- treating candidate housekeeping loci as sufficient normalization strategy across the genus
