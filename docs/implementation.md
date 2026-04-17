# Implementation Notes

This scaffold is intentionally split so each evidence tier can be upgraded independently.

## Current Module Map

| Module | Purpose now | Real implementation target |
|---|---|---|
| `VALIDATE_INPUTS` | schema and uniqueness checks | retain as-is |
| `PLAN_RECOVERY` | classify evidence tiers from manifest fields | retain as-is, optionally enrich with automated discovery |
| `RECOVER_TIER_AB` | emit contract outputs for assembly-backed species | add assembly fetch, BUSCO, Cactus, TOGA, CDS extraction, OrthoFinder audit |
| `RECOVER_TIER_C` | emit contract outputs for RNA-backed species | add read fetch, QC, transcript assembly or reference-guided CDS recovery, ORF calling |
| `RECOVER_TIER_D` | emit contract outputs for WGS-only species | add targeted mapping, variant calling, consensus generation, coverage filters |
| `REPORT_TIER_E` | emit explicit missing-data rows | retain as-is |
| `MERGE_RECOVERIES` | merge tier outputs into long, wide, and FASTA contracts | retain as-is |
| `RENDER_REPORT` | markdown summary | retain as-is or replace with richer reporting |

## Recommended Real Tooling

### Tier A and Tier B

Suggested chain:

1. `datasets` or direct FTP fetch for assemblies and metadata
2. `BUSCO` on genome and projected annotation
3. `Progressive Cactus` for whole-genome alignment
4. `TOGA` for ortholog projection
5. CDS and protein extraction
6. `OrthoFinder` for orthogroup validation

### Tier C

Suggested chain:

1. `prefetch` or ENA fetch
2. `fastp` for read QC
3. `rnaSPAdes` or another transcript assembly path, or direct reference-guided CDS recovery
4. `TransDecoder` or equivalent ORF recovery
5. orthogroup assignment against the reference panel

### Tier D

Suggested chain:

1. WGS read fetch
2. mapping to the closest `Phylloscopus` reference ortholog panel
3. `bcftools mpileup` and `bcftools call`
4. consensus generation
5. per-locus depth and allele-balance filtering

## Why The Pipeline Starts From A Manifest

The scaffold treats taxonomy and public-data discovery as an upstream refresh step. That is deliberate:

- taxonomy changes over time
- public accessions are unstable and can expand mid-project
- reconstruction should run on a frozen inventory for reproducibility

In production, a separate discovery builder should update the manifest on a schedule and commit dated snapshots.

