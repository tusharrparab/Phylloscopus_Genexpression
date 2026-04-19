# Phylloscopus Comparative Ortholog Recovery Report

## Run Summary

- Species in manifest: 2
- Ortholog targets: 1
- References: 1

## Tier Distribution

- Tier A: 1 species
- Tier B: 1 species
- Tier C: 0 species
- Tier D: 0 species
- Tier E: 0 species

## Target Panel

- Phylogenetic backbone targets: 1
- Migration-linked candidate targets: 0
- Vocalization/neural candidate targets: 0
- Hypoxia/elevation candidate targets: 0
- Housekeeping controls: 0
- Single-copy-preferred targets: 1
- Targets requiring explicit paralog screening: 0

## Recovery Totals

- Tier A recovered: 1
- Tier B candidate_sequence_recovered: 1
- Recovered sequences: 1
- Real reconstructed sequences: 0
- Candidate assembly-backed sequences: 1
- Stub sequences emitted: 0
- Planned only: 0
- Missing data: 0

## Matrix Preview

| species_id | scientific_name | evidence_tier |
|---|---|---|
| phylloscopus_collybita | Phylloscopus collybita | B |
| phylloscopus_trochilus | Phylloscopus trochilus | A |

## Notes

- This repository is a comparative ortholog-recovery scaffold, not a completed biological inference pipeline.
- `stub` mode proves workflow wiring only. Synthetic FASTA emitted by stub modules is not biological sequence.
- Tier B candidate sequences are assembly-backed locus recoveries only and are not validated ortholog calls.
- Tier C and Tier D outputs remain hypothesis-level placeholders until orthology is validated from real RNA or WGS data.
- Candidate-gene panel categories support hypothesis generation only; they do not by themselves establish causal trait associations.
- Single-copy-preferred targets still require locus-level inspection in the chosen reference and query assemblies.
- Consolidated ortholog FASTA bundles are written under `/tmp/phyllo-tierb-validation-local-20260418/merged/ortholog_sequences`.

## Interpretation Guardrails

- Non-final target rows (planned or missing): 0
- Treat downstream ASR as defensible only for loci reconstructed from real homologous sequence, not from stub outputs.
- Treat expression outputs as optional RNA-backed side analyses rather than the defining product of this repository.

