# Implementation Map

This file describes the real behavior of the repository as it exists now.

## Module Status Map

| Module | Implemented now | Still scaffolded |
|---|---|---|
| `VALIDATE_INPUTS` | schema checks, uniqueness checks, ortholog target metadata checks | no deep semantic validation of biological plausibility |
| `PLAN_RECOVERY` | evidence-tier planning and per-species recovery strategy | no automated locus-specific readiness assessment |
| `RECOVER_TIER_AB` | assembly staging, optional `datasets`, optional `BUSCO`, TOGA prerequisite planning | no chain building, no TOGA execution, no projected CDS extraction |
| `RECOVER_TIER_C` | stub or contract outputs only | no real RNA-based ortholog recovery |
| `RECOVER_TIER_D` | stub or contract outputs only | no real WGS-based ortholog recovery |
| `REPORT_TIER_E` | explicit missing-data accounting | nothing major missing for its narrow purpose |
| `MERGE_RECOVERIES` | status merging and per-locus FASTA merging | no occupancy filtering or orthology QC summaries |
| `RENDER_REPORT` | honest markdown report | no richer QC dashboard |
| `RUN_ASR` | optional `mafft`, optional `IQ-TREE`, plan-script generation | no codon-aware alignment or tree conflict checks |
| `RUN_EXPRESSION` | run metadata enrichment, transcript build, `salmon` planning or execution | no cross-species expression interpretation |

## Why Tier A/B Is Ahead Of Tier C/D

Assembly-backed taxa are the only part of the repository where the required assets are close enough to standard tool expectations that partial real execution is already practical. RNA-only and WGS-only tiers still need the hard biological work:

- ortholog assignment under incomplete evidence
- paralog control
- locus-specific QC
- reference bias management

That is why Tier C and Tier D remain scaffolded instead of pretending to be solved.

## Real Tools Already Wired

- `datasets`
- `BUSCO`
- `gffread`
- `mafft`
- `IQ-TREE`
- `prefetch`
- `fasterq-dump`
- `salmon`

Their presence improves staging or execution. Their presence does not by itself make the repository biologically validated.

## Honest Reading Of Outputs

- `ortholog_status.tsv` files track planning and evidence, not guaranteed orthology
- merged FASTA from stub mode is workflow test material only
- ASR outputs are only meaningful when the upstream loci are real validated homologs
- expression outputs are optional support products, not the central biological deliverable
