# Taxonomy And Discovery Builder

The recovery scaffold expects a frozen `species_manifest.tsv`. This builder creates that manifest from live public sources.

## What It Does

1. Resolves the target genus in GBIF.
2. Downloads the accepted species list for that genus.
3. For each species, queries NCBI Assembly, BioProject, and SRA.
4. Optionally queries ENA to distinguish transcriptomic versus genomic runs.
5. Writes a dated `species_manifest.tsv` plus reproducibility metadata.

## Outputs

- `species_manifest.tsv`
- `taxonomy_source.json`
- `snapshot_date`
- `discovery_summary.tsv`
- `raw/gbif_species.json`
- `raw/ncbi/*.json`
- `raw/ena/*.json` when `--include-ena` is enabled

## Usage

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

## Tier Rules Used By The Builder

- `B`: at least one NCBI assembly accession exists
- `C`: no assembly, but ENA transcriptomic runs were found
- `D`: no assembly, but ENA WGS-like genomic runs were found, or NCBI SRA has records with unresolved library typing
- `E`: no assembly and no recoverable public run evidence was found

The builder does **not** auto-assign `A`. Tier `A` requires a usable assembly plus an annotation path that the recovery pipeline can directly consume. In practice that usually means you stage or generate annotation files after the snapshot step.

## Feeding The Manifest Into Recovery

Once you have a snapshot, point the pipeline at it:

```bash
./nextflow run . \
  -profile local \
  --species_manifest snapshots/phylloscopus_2026-04-17/species_manifest.tsv \
  --ortholog_targets examples/ortholog_targets.tsv \
  --reference_manifest examples/reference_manifest.tsv
```

## Recommended Next Manual Step

After building the species snapshot:

1. Inspect `species_manifest.tsv`.
2. Promote one or more strong assembly species into `reference_manifest.tsv`.
3. Replace the example ortholog panel with your real single-copy target panel.
4. Only then swap the `stub` tier recovery modules for real download, QC, and inference commands.

