# Phylloscopus Ortholog Pipeline Feasibility

## Executive Summary

Yes, we can build a highly automated pipeline to infer orthologous loci and reconstruct ortholog sequences across the genus `Phylloscopus`, but not a truly complete push-button pipeline for **all** species using public data alone. The limiting factor is not the workflow engine; it is source-data coverage.

As of April 16, 2026, GBIF returns **74 accepted species** under `Phylloscopus`, while NCBI assembly search returns **13 assemblies representing 7 unique species** (`P. collybita`, `P. ibericus`, `P. plumbeitarsus`, `P. schwarzi`, `P. trochiloides`, `P. trochilus`, `P. whistleri`). NCBI also returns **51 BioProjects** and **586 SRA records** for the genus, which means some additional species likely have raw data, but not genus-wide, analysis-ready genomes. That makes a fully automated reconstruction pipeline feasible for a **data-available subset**, and feasible for the full genus only as a **tiered system** that reports missing species instead of pretending to reconstruct them ([GBIF genus query](https://api.gbif.org/v1/species/search?rank=SPECIES&highertaxon_key=2493047&status=ACCEPTED&limit=200), [NCBI assembly query](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=assembly&term=Phylloscopus%5BOrganism%5D&retmode=json), [NCBI BioProject query](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=bioproject&term=Phylloscopus%5BOrganism%5D&retmode=json), [NCBI SRA query](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=Phylloscopus%5BOrganism%5D&retmode=json)).

## Bottom Line

If by "create orthologs" you mean **infer orthologous genes/loci and recover their sequences**, then:

- **Yes** for species with a public genome assembly.
- **Usually yes** for species with transcriptomes or sufficiently deep raw reads.
- **Sometimes** for species with only sparse locus-level data.
- **No** for species with no recoverable public sequence data.

So the right engineering goal is not "fully automated for all species no matter what." The right goal is: **fully automated data discovery, QC, ortholog inference, sequence recovery, and confidence scoring, with explicit failure states for missing taxa**.

## Recommended Technical Design

### 1. Freeze the species list

Build a reproducible taxonomy snapshot step that writes a dated species manifest from GBIF or a chosen avian authority. This matters because the genus is not a static target.

Outputs:

- `species_manifest.tsv`
- `taxonomy_source.json`
- `snapshot_date`

### 2. Inventory available sequence data per species

For each species in the manifest, automatically query:

- NCBI Assembly
- NCBI BioProject
- NCBI SRA
- ENA as a secondary mirror if needed

Classify each species into evidence tiers:

| Tier | Input data | What the pipeline can do |
|---|---|---|
| A | Annotated or high-quality genome assembly | Genome-wide ortholog inference and sequence extraction |
| B | Genome assembly without good annotation | Reference-guided annotation plus ortholog inference |
| C | Transcriptome or RNA-seq | Transcript assembly or direct CDS recovery for many genes |
| D | Raw WGS only | Consensus recovery for targeted loci, lower confidence |
| E | No usable public data | Report as unavailable; no reconstruction |

### 3. Use a reference-based ortholog recovery backbone

For genome-bearing species, the best core is:

1. Choose one or more high-quality bird references.
2. Build pairwise or multi-genome alignments.
3. Project coding genes and infer orthologous loci.

`TOGA` is the strongest fit here because it combines gene projection and orthology inference from whole-genome alignments, scales to hundreds of genomes, and already demonstrated use on **501 bird assemblies** ([TOGA paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC10193443/)). For alignment, `Progressive Cactus` is the right large-scale whole-genome alignment backbone ([Nature paper](https://www.nature.com/articles/s41586-020-2871-y)).

### 4. Add an orthogroup layer for cross-checking

After projected gene models are produced, run a protein-level orthology layer such as `OrthoFinder` on predicted proteomes. This gives:

- orthogroups
- gene trees
- duplication calls
- cross-checking against projection-based orthology

This is useful because projection methods are strongest for conserved loci and synteny-aware inference, while orthogroup clustering is useful for auditing unexpected duplications and split models ([OrthoFinder 2019](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1832-y)).

### 5. Recover sequences by evidence tier

Recommended recovery logic:

- **Tier A/B genomes**: TOGA projection -> CDS/protein extraction -> ortholog table -> codon alignments.
- **Tier C transcriptomes/RNA-seq**: transcript assembly or direct mapping to reference ortholog CDS -> ORF calling -> orthogroup assignment.
- **Tier D raw WGS**: map reads to closest `Phylloscopus` reference ortholog loci -> call consensus only where depth and allele balance pass thresholds.
- **Tier E no data**: emit `missing_data`, not an imputed sequence.

Important constraint: do **not** silently impute full sequences for absent species from sister taxa. That produces something, but it is not a biologically observed ortholog sequence.

### 6. Build QC into every stage

QC should be mandatory, not optional:

- BUSCO completeness on assemblies, annotations, and transcriptomes ([BUSCO user guide](https://busco.ezlab.org/busco_userguide), [BUSCO paper](https://academic.oup.com/bioinformatics/article/31/19/3210/211866))
- stop-codon and frameshift checks
- exon boundary consistency
- synteny support for projected loci
- read-depth thresholds for consensus reconstructions
- ortholog class support (`1:1`, `1:many`, `many:1`, `many:many`)
- per-sequence confidence score

Final outputs should always carry provenance:

- source database accession
- reconstruction method
- QC status
- confidence tier

## What "Fully Automated" Should Mean Here

A credible fully automated pipeline should automate:

- taxonomy snapshotting
- public-data discovery
- download and metadata normalization
- assembly/transcriptome QC
- ortholog inference
- sequence extraction or reconstruction
- alignment generation
- confidence scoring
- report generation

It should **not** claim to automate away biological absence. For some `Phylloscopus` species, the automated result today will correctly be: "no recoverable nuclear ortholog sequence because no suitable public data were found."

## Practical Recommendation

The right first version is:

1. Target a **single-copy coding ortholog panel**, not the entire gene space.
2. Support all five evidence tiers.
3. Produce a genus-wide status matrix before attempting any large reconstruction run.
4. Start with the species that already have assemblies, then expand to transcriptome/read-backed species.

That gives a scientifically defensible pipeline quickly, instead of over-promising whole-genome ortholog catalogs for taxa that currently have no usable inputs.

## Proof-of-Concept Status

The current scaffold now includes a minimal empirical proof-of-concept for Tier B: one assembly-backed candidate sequence recovery has been demonstrated for a single conservative backbone locus from an unannotated assembly. In the archived `RAG1` example, the selected assembly hit spans the full search query length used for retrieval, but the recovered region is still shorter than the nominal target CDS length and remains a candidate assembly-backed locus recovery with unvalidated orthology and structure. This is intentionally narrow and should be interpreted only as a feasibility check for candidate locus recovery, not as a validated orthology result or a general Tier B ortholog workflow.

## Decision

**Recommendation:** build it, but define the product as a **tiered automated ortholog recovery system** rather than a universal all-species reconstruction machine.

If you want, the next step should be a concrete implementation spec in either `Nextflow` or `Snakemake`, including:

- directory layout
- manifest schemas
- module boundaries
- container/tool choices
- output contracts
- failure and retry rules

## Sources

- GBIF accepted-species query for `Phylloscopus`: [https://api.gbif.org/v1/species/search?rank=SPECIES&highertaxon_key=2493047&status=ACCEPTED&limit=200](https://api.gbif.org/v1/species/search?rank=SPECIES&highertaxon_key=2493047&status=ACCEPTED&limit=200)
- NCBI Assembly search for `Phylloscopus`: [https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=assembly&term=Phylloscopus%5BOrganism%5D&retmode=json](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=assembly&term=Phylloscopus%5BOrganism%5D&retmode=json)
- NCBI BioProject search for `Phylloscopus`: [https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=bioproject&term=Phylloscopus%5BOrganism%5D&retmode=json](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=bioproject&term=Phylloscopus%5BOrganism%5D&retmode=json)
- NCBI SRA search for `Phylloscopus`: [https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=Phylloscopus%5BOrganism%5D&retmode=json](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=Phylloscopus%5BOrganism%5D&retmode=json)
- TOGA: [https://pmc.ncbi.nlm.nih.gov/articles/PMC10193443/](https://pmc.ncbi.nlm.nih.gov/articles/PMC10193443/)
- Progressive Cactus: [https://www.nature.com/articles/s41586-020-2871-y](https://www.nature.com/articles/s41586-020-2871-y)
- OrthoFinder: [https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1832-y](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1832-y)
- BUSCO user guide: [https://busco.ezlab.org/busco_userguide](https://busco.ezlab.org/busco_userguide)
- BUSCO paper: [https://academic.oup.com/bioinformatics/article/31/19/3210/211866](https://academic.oup.com/bioinformatics/article/31/19/3210/211866)
