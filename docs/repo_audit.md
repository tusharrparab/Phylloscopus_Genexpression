# Repository Audit

## Bottom Line

This repository is scientifically more honest as a `Phylloscopus` comparative ortholog-recovery scaffold than as a "gene expression" repository. The current codebase can validate manifests, plan tier-specific recovery, stage some real assembly-backed work, run optional ASR tooling, and stage optional RNA quantification. It cannot yet support publication-grade orthology calls, clade-wide comparative inference, or trait-linked conclusions.

## Already Implemented

- manifest validation and uniqueness checks
- frozen-species discovery from GBIF, NCBI, and ENA
- tier assignment from staged evidence
- Tier A/B assembly staging with `datasets`, `BUSCO`, and explicit TOGA prerequisite plans
- merged per-gene recovery matrices and merged FASTA bundles
- ASR alignment and `IQ-TREE` execution when tools are present
- expression metadata enrichment, transcript construction, and `salmon` quant planning or execution when tools are present

## Scaffolded Only

- Tier C ortholog recovery from RNA evidence
- Tier D ortholog recovery from WGS evidence
- full TOGA execution
- ortholog vetting beyond target-list intent
- clade-wide ortholog matrix curation
- cross-species expression interpretation

## Biologically Implied But Not Yet Supported

- that any recovered tier output is a validated ortholog
- that candidate genes are causal for migration, vocalization, or elevation adaptation
- that ASR on stub-generated or weakly validated loci is meaningful
- that RNA-backed quantification is comparable across species without ortholog-level aggregation and normalization design
- that backbone markers alone are sufficient for robust species-tree inference across the full genus

## Main Scientific Problems Corrected In This Revision

- The repository no longer presents itself primarily as an expression project.
- The target panel is no longer an unstructured gene list.
- Ortholog targets now carry category, rationale, and orthology-risk metadata.
- Tier evidence is now framed as evidence strength, not as proof of orthology.
- Documentation now separates implemented behavior from planned upgrades.

## Remaining Hard Problems

- no gene-tree or synteny-based ortholog validation
- no paralog filtering workflow
- no codon-aware comparative alignment path
- no species-tree estimation or ortholog occupancy filtering pipeline
- no publication-ready treatment of batch, tissue, season, or developmental heterogeneity for expression
