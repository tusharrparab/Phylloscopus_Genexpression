# Tier A/B Scaffold

Tier `A/B` now has a concrete scaffold in [run_tier_ab_scaffold.py](/Users/sam/Documents/New%20project/bin/run_tier_ab_scaffold.py).

## What It Does

For assembly-backed species:

1. resolves local assembly assets when they already exist
2. downloads NCBI genome packages with `datasets` when an assembly accession is present and the CLI is installed
3. runs `BUSCO` in genome mode when `busco` is installed
4. writes a TOGA projection handoff plan for each query species

## Current Behavior

- If `datasets` is available, genome packages are downloaded with `datasets download genome accession ... --include genome,gff3,gtf,protein,cds` ([NCBI docs](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/datasets/download/genome/)).
- If `busco` is available, the scaffold runs genome mode with either the configured lineage or `--auto-lineage-euk` ([BUSCO user guide](https://busco.ezlab.org/busco_userguide)).
- If TOGA prerequisites are incomplete, the pipeline still writes a per-species `run_toga.sh` and `projection_plan.json` so the missing pieces are explicit instead of implicit.

## Important Limitation

TOGA itself is not executed automatically yet. That is deliberate.

According to the TOGA README, the real mandatory inputs are:

- a reference `bed12`
- a chain file
- reference and query genomes in `2bit`

Those assets are not derivable from an assembly accession alone. The scaffold now prepares everything up to that boundary and writes the exact missing prerequisites into `projection_plan.tsv` and each species `projection_plan.json` ([TOGA README](https://github.com/hillerlab/TOGA)).

## Outputs

Under `results/tier_ab/`:

- `ortholog_status.tsv`
- `reference_assets.tsv`
- `projection_plan.tsv`
- `assets/references/...`
- `assets/queries/...`
- `assets/queries/<species>/projection/run_toga.sh`

