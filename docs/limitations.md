# Limitations

## Methodological Limits

- Tier C and Tier D are not real recovery pipelines yet.
- Stub FASTA output is synthetic and must never be analyzed as biological sequence.
- Ortholog targets are design intents, not validated one-to-one orthologs.
- The workflow does not yet include paralog filtering, orthogroup inference, or gene-tree QC.
- Missing data are tracked, but locus occupancy thresholds are not yet enforced.

## Biological Limits

- Trait-linked target categories are hypothesis-linked only.
- Candidate genes for migration, vocalization, or hypoxia are not evidence of adaptation on their own.
- Expression outputs are not yet suitable for comparative interpretation across species.
- ASR is only as defensible as the orthology and alignment quality of the input loci.

## Data Limits

- public metadata can be incomplete or wrong
- ENA/SRA library layout may still be unresolved for some runs
- assembly accessions do not guarantee usable annotation assets
- different species will inevitably enter the workflow with different evidence quality

## Publication-Grade Requirements Still Missing

- validated locus selection
- explicit paralog rejection
- high-quality comparative alignments
- species-tree support and conflict assessment
- reproducible comparative statistics tied to a biological study design
