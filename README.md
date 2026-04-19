# Phylloscopus Comparative Ortholog Recovery Scaffold

It is a scaffold that separates:

- what is already executable,
- what is only staged or planned,
- what still requires real orthology validation and comparative analysis.

## What The Repository Actually Implements

- manifest validation for species, targets, and references
- tier assignment across species with uneven data availability
- assembly-backed Tier A/B staging, including NCBI `datasets` download when available, `BUSCO` execution when installed, and explicit TOGA handoff plans
- stub-mode tier recovery that proves workflow plumbing end to end
- merged ortholog status matrices and per-locus FASTA bundles
- optional ASR staging, with real `mafft` and `IQ-TREE` execution if those tools are installed
- optional RNA-backed expression staging, with real `gffread`, `salmon`, and SRA/ENA download planning when tools and metadata are available

## What The Repository Does Not Yet Implement

- validated ortholog inference across all tiers
- paralog disambiguation or gene-tree reconciliation
- real Tier C transcript reconstruction
- real Tier D targeted consensus recovery
- automatic TOGA execution with chains and `2bit` assets
- comparative hypothesis testing, trait association, or publication-grade phylogenomic inference
- defensible cross-species expression normalization or differential expression analysis

Those gaps are intentional and are documented explicitly in:

- [docs/repo_audit.md](/Users/sam/Documents/New%20project/docs/repo_audit.md)
- [docs/biological_scope.md](/Users/sam/Documents/New%20project/docs/biological_scope.md)
- [docs/ortholog_strategy.md](/Users/sam/Documents/New%20project/docs/ortholog_strategy.md)
- [docs/expression_scope.md](/Users/sam/Documents/New%20project/docs/expression_scope.md)
- [docs/limitations.md](/Users/sam/Documents/New%20project/docs/limitations.md)

## Biological Framing

The repository is organized around a biologically conservative idea:

- recover or stage ortholog evidence across a clade with uneven genomic resources
- keep evidence strength explicit by tier
- separate backbone markers from trait-linked candidate genes
- treat migration, vocalization, and elevation physiology as hypothesis-linked use cases, not completed findings

The example target panel is now structured into:

- `phylogenetic_backbone`
- `migration_candidate`
- `vocalization_neural_candidate`
- `hypoxia_elevation_candidate`
- `housekeeping_control`

Each target includes a short rationale plus an orthology-risk annotation. Candidate genes are included as hypotheses only. They are not treated as proof of trait causation.

## Evidence Tiers

- `A`: assembly and annotation already staged locally
- `B`: assembly-backed, but annotation or projection assets are incomplete
- `C`: transcriptome or RNA-seq evidence only
- `D`: WGS evidence only
- `E`: no usable public sequence evidence currently staged

These tiers are useful for planning. They are not equivalent evidence classes for downstream inference.

## Input Manifests

- [examples/species_manifest.tsv](/Users/sam/Documents/New%20project/examples/species_manifest.tsv)
- [examples/ortholog_targets.tsv](/Users/sam/Documents/New%20project/examples/ortholog_targets.tsv)
- [examples/reference_manifest.tsv](/Users/sam/Documents/New%20project/examples/reference_manifest.tsv)

The example manifests are designed to be honest rather than polished. They include:

- evidence confidence
- provenance
- analysis suitability notes
- ortholog target category and rationale
- explicit orthology expectations

Schema details are in [docs/manifest-schemas.md](/Users/sam/Documents/New%20project/docs/manifest-schemas.md).

## Quick Start

Run the full scaffold in `stub` mode:

```bash
./nextflow run . -profile local -ansi-log false
```

This validates the manifests, assigns tiers, emits stub recovery outputs, merges status tables, and writes a markdown report. Stub FASTA output is synthetic and is only for workflow validation.

## Real And Partially Real Components

### Discovery

Build a dated species snapshot from GBIF, NCBI, and optionally ENA:

```bash
python3 bin/build_species_manifest.py \
  --genus Phylloscopus \
  --include-ena \
  --outdir snapshots/phylloscopus_2026-04-17
```

### Reference Selection

Promote assembly-backed taxa into a reference manifest:

```bash
python3 bin/build_reference_manifest.py \
  --species-manifest snapshots/phylloscopus_2026-04-17/species_manifest.tsv \
  --out snapshots/phylloscopus_2026-04-17/reference_manifest.tsv
```

### Tier A/B Staging

Tier A/B can already stage:

- `datasets` assembly downloads
- `BUSCO` genome runs
- TOGA prerequisite planning

It still does not execute full projection.

### ASR

ASR can already use `mafft` and `IQ-TREE` if they are installed:

```bash
python3 bin/run_asr_scaffold.py \
  --sequence-dir examples/asr/ortholog_sequences \
  --species-tree examples/asr/phylloscopus_example_tree.nwk \
  --outdir results/asr \
  --mode scaffold
```

Interpret ASR outputs biologically only when the input loci are real homologous sequences. Stub-generated FASTA is not valid input for inference.

### Expression

Expression is optional and secondary. It is intended for RNA-backed taxa once ortholog recovery already has a comparative context.

The expression scaffold can already:

- resolve run metadata from local `run_metadata.tsv`, ENA, and NCBI SRA RunInfo
- build transcript FASTA from a genome plus annotation with `gffread`
- build a `salmon` index
- plan or execute `salmon quant`

It still does not solve cross-species normalization or comparative expression interpretation.

```bash
python3 bin/run_expression_scaffold.py \
  --species-manifest snapshots/phylloscopus_2026-04-17/species_manifest.tsv \
  --reference-manifest snapshots/phylloscopus_2026-04-17/reference_manifest.tsv \
  --run-metadata snapshots/phylloscopus_2026-04-17/run_metadata.tsv \
  --outdir results/expression \
  --mode execute
```

## Repository Layout

```text
.
├── bin/                  Python entrypoints and helper scripts
├── docs/                 biological framing, audit, limitations, and module notes
├── envs/                 optional conda runtime definition
├── examples/             example manifests and toy inputs
├── modules/local/        Nextflow process wrappers
├── main.nf               top-level workflow
└── nextflow.config       runtime configuration
```

## Design Principles

- be explicit about evidence quality
- prefer single-copy orthologs for backbone inference
- treat trait-linked candidate genes as candidates, not claims
- keep discovery snapshots frozen and auditable
- preserve stub executability without pretending stub outputs are biology

## Further Reading

- [docs/repo_audit.md](/Users/sam/Documents/New%20project/docs/repo_audit.md)
- [docs/implementation.md](/Users/sam/Documents/New%20project/docs/implementation.md)
- [docs/tier-ab-scaffold.md](/Users/sam/Documents/New%20project/docs/tier-ab-scaffold.md)
- [docs/asr-scaffold.md](/Users/sam/Documents/New%20project/docs/asr-scaffold.md)
- [docs/expression-scaffold.md](/Users/sam/Documents/New%20project/docs/expression-scaffold.md)
