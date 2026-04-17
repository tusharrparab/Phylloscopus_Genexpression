# Ortholog Strategy

## Core Principle

A target list is not an ortholog callset.

This repository currently defines a candidate panel and recovery plan. It does not yet prove one-to-one orthology across all taxa. That distinction matters most for Tier C and Tier D taxa, where sequence evidence is incomplete and paralog confusion is a realistic failure mode.

## Current Strategy

1. Define a biologically motivated target panel.
2. Annotate each target with an orthology expectation.
3. Prefer single-copy loci for phylogenetic backbone analyses.
4. Keep candidate-gene classes explicit when single-copy status is less certain.
5. Track missingness and evidence tier across species.

## Current Orthology Metadata

`examples/ortholog_targets.tsv` now includes:

- `category`
- `orthology_basis`
- `copy_number_expectation`
- `rationale`

These fields are not cosmetic. They are there to prevent the panel from being treated as if every locus carries the same inferential weight.

## Recommended Validation Path Before Biological Inference

### Backbone Loci

- confirm reference annotation quality
- confirm expected exon structure
- verify one-to-one orthology with reciprocal best hits plus synteny where possible
- inspect alignments for frame integrity and unexpected indels

### Candidate Loci

- inspect copy number in the primary reference
- inspect whether close paralogs exist in other avian genomes
- prefer gene-tree or orthogroup-based validation before sequence-evolution claims
- require stronger manual review for loci marked `screen_for_paralogs`

## Why Single-Copy Preference Matters

Backbone phylogenetic inference is much more stable when loci are conservative and low-copy. That is why the example panel distinguishes `single_copy_preferred` loci from candidates that need explicit paralog checks.

Historical avian marker use is part of the motivation, not proof of universal safety. `RAG1`, for example, is attractive because of its conservative evolutionary behavior and broad phylogenetic use in birds, but even classical markers still need reference-specific inspection in a modern comparative pipeline: [Groth and Barrowclough 1999](https://pubmed.ncbi.nlm.nih.gov/10381315/).

## What Is Still Missing

- orthogroup inference from proteomes
- reciprocal-hit automation
- synteny-aware validation
- gene-tree/species-tree reconciliation
- locus occupancy filters and orthology QC dashboards
