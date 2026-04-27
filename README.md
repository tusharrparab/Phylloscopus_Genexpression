# Phylloscopus Tier-Aware Ortholog Recovery Pipeline
This repository contains a tier-aware ortholog recovery framework for `Phylloscopus` species under uneven data availability. It is organized as a pipeline repository with retained validation artifacts for the currently implemented recovery paths.

## Implemented Functionality

- manifest validation for species, targets, and references
- evidence-tier assignment from staged assets and accessions
- Tier A reference-backed transcript extraction when a staged annotation contains the requested target
- Tier A/B assembly staging, including optional NCBI `datasets` download, optional `BUSCO`, and explicit TOGA prerequisite planning
- merged ortholog status tables and per-locus FASTA bundles
- optional downstream ASR staging with `mafft` and `IQ-TREE` when installed
- optional RNA-backed expression staging with `gffread`, `salmon`, and run-metadata enrichment when installed

## Partial Proof Of Concept

- Tier A has real non-stub recovery for locally staged annotated references
- Tier B has one minimal archived candidate recovery for `RAG1` in `phylloscopus_collybita`
- the Tier B result is an assembly-backed candidate sequence recovery only; it is not a validated ortholog call
- the retained validation package is under [validation/minimal_tier_b_candidate_recovery](/Users/sam/Documents/New project/validation/minimal_tier_b_candidate_recovery)

## Not Implemented

- general Tier B ortholog validation from unannotated assemblies
- real Tier C transcript-based ortholog recovery
- real Tier D WGS-based ortholog recovery
- paralog disambiguation, orthogroup inference, or gene-tree reconciliation
- full TOGA execution with chains and `2bit` assets
- comparative inference built on unvalidated candidate sequences
- cross-species expression normalization or comparative differential expression analysis

## Validation Artifacts

- `validation/minimal_tier_b_candidate_recovery` retains the current narrow validation package
- retained artifacts include merged status tables, merged `RAG1` FASTA output, Tier B candidate summary tables, and the second-best-hit sanity table
- the Tier A vs Tier B pairwise check is retained only as [validation/minimal_tier_b_candidate_recovery/summary_tables/non_informative_pairwise_check.tsv](/Users/sam/Documents/New project/validation/minimal_tier_b_candidate_recovery/summary_tables/non_informative_pairwise_check.tsv) because the Tier A comparator is a short mock reference transcript and the result is not biologically informative

## Technical Documentation

- [docs/implementation.md](/Users/sam/Documents/New project/docs/implementation.md)
- [docs/manifest-schemas.md](/Users/sam/Documents/New project/docs/manifest-schemas.md)
- [docs/tier-ab-scaffold.md](/Users/sam/Documents/New project/docs/tier-ab-scaffold.md)
- [docs/asr-scaffold.md](/Users/sam/Documents/New project/docs/asr-scaffold.md)
- [docs/expression-scaffold.md](/Users/sam/Documents/New project/docs/expression-scaffold.md)
- [docs/discovery-builder.md](/Users/sam/Documents/New project/docs/discovery-builder.md)
- [docs/limitations.md](/Users/sam/Documents/New project/docs/limitations.md)

## Input Manifests

- [examples/species_manifest.tsv](/Users/sam/Documents/New project/examples/species_manifest.tsv)
- [examples/ortholog_targets.tsv](/Users/sam/Documents/New project/examples/ortholog_targets.tsv)
- [examples/reference_manifest.tsv](/Users/sam/Documents/New project/examples/reference_manifest.tsv)

## Quick Start

Run the workflow in `stub` mode:

```bash
./nextflow run . -profile local -ansi-log false
```

This validates manifests, assigns tiers, emits stub recovery outputs for unimplemented routes, merges status tables, and writes a run report. Stub FASTA output is synthetic and should be treated as workflow-validation output only.

## Repository Layout

```text
.
├── bin/              Python entrypoints and helper scripts
├── docs/             technical documentation
├── envs/             optional runtime definition
├── examples/         example manifests and toy Tier A reference inputs
├── modules/local/    Nextflow process wrappers
├── validation/       retained validation artifacts
├── main.nf           top-level workflow
└── nextflow.config   runtime configuration
```
