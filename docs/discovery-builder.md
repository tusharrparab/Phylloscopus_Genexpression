# Taxonomy And Discovery Builder

The recovery workflow should not operate on an unfrozen view of public databases. `build_species_manifest.py` exists to produce dated species snapshots and expose evidence heterogeneity before any recovery attempt.

## What It Does

1. resolves the target genus in GBIF
2. fetches accepted species for that genus
3. queries NCBI Assembly, BioProject, and SRA
4. optionally queries ENA read-run metadata
5. writes a dated `species_manifest.tsv` plus reproducibility metadata

## What It Does Not Do

- prove that a taxon is ready for biological inference
- auto-assign Tier `A`
- validate orthology
- guarantee that all accessions are biologically usable

## Tier Rules Used By The Builder

- `B`: at least one assembly accession exists
- `C`: no assembly, but transcriptomic runs exist
- `D`: no assembly, but WGS-like or unresolved run evidence exists
- `E`: no public evidence useful to this scaffold was detected

Tier `A` still requires local assembly plus annotation assets staged in a way the recovery workflow can actually consume.

## Output Files

- `species_manifest.tsv`
- `taxonomy_source.json`
- `snapshot_date`
- `discovery_summary.tsv`
- `run_metadata.tsv` when `--include-ena` is enabled
- `raw/gbif_species.json`
- `raw/ncbi/*.json`
- `raw/ena/*.json` when `--include-ena` is enabled

## Notes On Honesty

- `evidence_confidence` and `analysis_suitability` are planning fields, not proof fields
- `run_metadata.tsv` improves run-level planning but does not solve biological comparability
- ENA and SRA metadata can be incomplete or contradictory

## Typical Usage

```bash
python3 bin/build_species_manifest.py \
  --genus Phylloscopus \
  --include-ena \
  --outdir snapshots/phylloscopus_2026-04-17
```

## Recommended Manual Next Steps

1. inspect the species snapshot critically
2. choose one or more assembly-backed taxa for the reference manifest
3. decide whether the ortholog target panel is appropriate for the biological question
4. only then run recovery in stub or staging mode
