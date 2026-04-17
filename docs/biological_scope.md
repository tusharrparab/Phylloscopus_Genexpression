# Biological Scope

## What This Scaffold Is For

This repository is designed to support a future comparative genomics workflow across `Phylloscopus` species with uneven data availability. The immediate goal is not to prove biological hypotheses. The immediate goal is to stage ortholog evidence honestly and traceably so later comparative analyses are possible.

## Biological Questions The Design Can Eventually Support

- reconstruction of an ortholog matrix with explicit missingness
- phylogenetic backbone inference from conservative nuclear loci
- ancestral sequence reconstruction on validated loci
- exploratory tests of whether trait-linked candidate genes show unusual sequence evolution in focal lineages
- optional RNA-backed expression quantification for taxa that only have transcriptomic evidence

## Biological Questions It Cannot Yet Support

- causal claims about migration genetics
- causal claims about vocal learning or acoustic divergence
- causal claims about high-elevation adaptation
- cross-species differential expression claims
- robust phylogenomic conclusions from current stub outputs

## Why The Target Panel Is Structured

The example panel mixes several target classes on purpose:

- `phylogenetic_backbone`: conservative loci with long use in avian systematics
- `migration_candidate`: genes previously discussed in avian migration candidate-gene literature
- `vocalization_neural_candidate`: genes with known relevance to songbird neural or vocal systems
- `hypoxia_elevation_candidate`: genes in hypoxia-response pathways repeatedly examined in altitude adaptation work
- `housekeeping_control`: genes useful as expression or assembly sanity checks, not as primary phylogenetic evidence

This structure is meant to keep hypothesis generation separated from ortholog reliability.

## Candidate-Gene Caution

Candidate genes are included because they are plausible biological entry points, not because they prove trait evolution. For migration in particular, published associations are mixed and context dependent. The panel therefore treats migration, neural, and hypoxia loci as targets for later validation and comparative screening, not as privileged truth.

## Literature Anchors For The Example Panel

- RAG1 has long-standing use as a conservative avian nuclear marker: [Groth and Barrowclough 1999](https://pubmed.ncbi.nlm.nih.gov/10381315/)
- Migration candidates such as `ADCYAP1` and `CLOCK` have been studied in passerines, but effects are not uniformly replicated: [Mueller et al. 2011](https://pubmed.ncbi.nlm.nih.gov/21325325/), [Peterson et al. 2013](https://pubmed.ncbi.nlm.nih.gov/24627781/)
- `FOXP2` is an established vocal-learning and songbird neural candidate: [Xiao et al. 2021](https://www.nature.com/articles/s41467-021-22918-2)
- HIF-pathway genes such as `EGLN1` and `EPAS1` are recurrent altitude-adaptation candidates in birds and other vertebrates: [Graham and McCracken 2019](https://pubmed.ncbi.nlm.nih.gov/30631144/)

## Appropriate Downstream Uses After Real Ortholog Validation

- occupancy-filtered locus matrices
- codon or amino-acid alignments
- branch-wise dN/dS or convergence screens
- validated ASR on clade subsets
- trait-linked comparative analyses with explicit phylogenetic control
