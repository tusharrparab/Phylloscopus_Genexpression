# Phylloscopus_Genexpression

Integrated local scaffold for `tusharrparab/Phylloscopus_Genexpression`, focused on tiered ortholog recovery and sequence reconstruction across `Phylloscopus`.

This repository contains a concrete `Nextflow` scaffold for a tiered ortholog-recovery pipeline across the genus `Phylloscopus`.

The current implementation is intentionally executable in `stub` mode so the project can be validated end-to-end without external bioinformatics tools. The module boundaries, manifests, outputs, and upgrade path are laid out so the scaffold can be hardened into a real production pipeline.

## What This Scaffold Does

- validates the species, target, and reference manifests
- classifies each species into evidence tiers `A` through `E`
- routes recovery through tier-specific modules
- emits ortholog status tables and mock FASTA bundles in `stub` mode
- merges all tier outputs into a final status matrix and report

## Evidence Tiers

- `A`: assembly + annotation available
- `B`: assembly available, annotation incomplete or missing
- `C`: transcriptome or RNA-seq only
- `D`: WGS reads only
- `E`: no usable public sequence input

## Repository Layout

```text
.
├── bin/
├── docs/
├── examples/
├── modules/local/
├── main.nf
└── nextflow.config
```

## Build A Real Species Snapshot

The missing upstream step in the scaffold is now implemented as a manifest builder that freezes the accepted `Phylloscopus` species list and inventories public sequence evidence.

Small validation run:

```bash
python3 bin/build_species_manifest.py \
  --genus Phylloscopus \
  --species-limit 5 \
  --include-ena \
  --outdir snapshots/phylloscopus_test_2026-04-17
```

Full genus snapshot:

```bash
python3 bin/build_species_manifest.py \
  --genus Phylloscopus \
  --include-ena \
  --outdir snapshots/phylloscopus_2026-04-17
```

See [docs/discovery-builder.md](/Users/sam/Documents/New%20project/docs/discovery-builder.md) for the builder outputs and downstream handoff.

## Quick Start

1. Bootstrap `Nextflow` locally if it is not already installed.
2. Run the scaffold in `stub` mode with the bundled example manifests.

```bash
./nextflow run . -profile local
```

Outputs will be published under `results/`.

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `params.species_manifest` | `examples/species_manifest.tsv` | species inventory plus available input metadata |
| `params.ortholog_targets` | `examples/ortholog_targets.tsv` | ortholog panel to recover |
| `params.reference_manifest` | `examples/reference_manifest.tsv` | reference assemblies and annotations |
| `params.outdir` | `results` | published output directory |
| `params.execution_mode` | `stub` | `stub` emits mock sequences; `contract` emits status-only contracts |

## Primary Outputs

- `results/validation/validation_summary.tsv`
- `results/planning/species_plan.tsv`
- `results/tier_ab/`
- `results/tier_c/`
- `results/tier_d/`
- `results/tier_e/`
- `results/merged/ortholog_status_long.tsv`
- `results/merged/ortholog_status_matrix.tsv`
- `results/merged/ortholog_sequences/`
- `results/reports/pipeline_report.md`

## Upgrading This Scaffold To Real Recovery

This repository already has the right split points for a full implementation:

- `RECOVER_TIER_AB`: add assembly download, BUSCO QC, Cactus alignment, TOGA projection, CDS extraction, OrthoFinder audit
- `RECOVER_TIER_C`: add read download, RNA-seq QC, transcript assembly or reference-guided CDS recovery, ORF calling, orthogroup assignment
- `RECOVER_TIER_D`: add read download, mapping to closest reference ortholog loci, consensus generation, depth and allele-balance filters
- `REPORT_TIER_E`: retain explicit missing-data accounting

The detailed mapping from scaffold modules to real tools is in [docs/implementation.md](/Users/sam/Documents/New%20project/docs/implementation.md).

## Example Run Contract

The bundled manifests intentionally span all evidence tiers so the example run produces:

- at least one `A` tier species
- multiple `B` tier species
- one `C` tier species
- one `D` tier species
- one `E` tier species

That makes it useful for verifying merge logic and downstream consumers before wiring in heavy external tools.
