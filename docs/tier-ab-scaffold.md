# Tier A/B Scaffold

Tier A/B is the most biologically mature part of the repository because assembly-backed taxa can be staged against standard tools.

## Implemented Now

- local assembly asset resolution
- NCBI `datasets` download when an assembly accession is available
- optional `BUSCO` genome mode runs
- explicit TOGA prerequisite planning per species

## Not Yet Implemented

- chain generation
- `2bit` conversion workflow
- full TOGA execution
- projected CDS extraction
- orthogroup validation after projection

## Why This Matters

Tier A/B is where the repository begins to transition from workflow scaffold to potentially useful comparative genomics staging. It is still not enough to claim validated ortholog recovery.

## Current Outputs

- `ortholog_status.tsv`
- `reference_assets.tsv`
- `projection_plan.tsv`
- downloaded or staged asset directories
- per-species `run_toga.sh` plans

## Interpretation

The presence of a TOGA handoff plan means the prerequisite accounting is explicit. It does not mean projection has been performed.
