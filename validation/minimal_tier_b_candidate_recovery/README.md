# Minimal Tier B Candidate Recovery

This validation package captures one narrow Tier B validation run in which a conservative backbone locus was recovered as an assembly-backed candidate sequence from an unannotated assembly and propagated through the existing merge and reporting layers.

Recovered locus

- target locus: `RAG1`
- Tier A species: `phylloscopus_trochilus`
- Tier B species: `phylloscopus_collybita`
- Tier B assembly cache used for the run: `GCA_056576225.1_ASM5657622v1_genomic.clean.fna.gz`

How it was recovered

- a real `Phylloscopus collybita` `RAG1` query sequence from GenBank record `AY319997.1` was used as the narrow reference query
- the Tier B scaffold searched the unannotated assembly for a best-hit contiguous region and extracted the best matching segment
- the recovered Tier B segment maps to `JBQOOF010106519.1:6876447-6879318(+)`
- recovered candidate length: `2872 bp`
- query coverage: `1.000`
- match identity: `0.996`
- the selected assembly hit spans the full search query length used for retrieval, but the recovered region is still shorter than the nominal target CDS length in the manifest
- under the current reporting thresholds, no second assembly hit passed the same sanity filter for this query
- the second-best-hit table is a sanity artifact only and does not convert this candidate into a validated ortholog

Interpretation guardrail

- the Tier B sequence is labeled as a `candidate assembly-backed locus sequence`
- it is marked `candidate_sequence=yes`
- it is marked `assembly_backed=yes`
- it is marked `orthology_unvalidated=yes`
- it must not be treated as a validated ortholog call

Package contents

- original validation manifests: [manifests](/Users/sam/Documents/New project/validation/minimal_tier_b_candidate_recovery/manifests)
- exact manifests used for the successful cached run: [manifests_used](/Users/sam/Documents/New project/validation/minimal_tier_b_candidate_recovery/manifests_used)
- merged status tables: [ortholog_status_long.tsv](/Users/sam/Documents/New project/validation/minimal_tier_b_candidate_recovery/ortholog_status_long.tsv), [ortholog_status_matrix.tsv](/Users/sam/Documents/New project/validation/minimal_tier_b_candidate_recovery/ortholog_status_matrix.tsv)
- sequence outputs: [ortholog_sequences/RAG1.fna](/Users/sam/Documents/New project/validation/minimal_tier_b_candidate_recovery/ortholog_sequences/RAG1.fna)
- summary tables: [summary_tables/tier_b_candidate_summary.tsv](/Users/sam/Documents/New project/validation/minimal_tier_b_candidate_recovery/summary_tables/tier_b_candidate_summary.tsv), [summary_tables/second_best_hit_sanity.tsv](/Users/sam/Documents/New project/validation/minimal_tier_b_candidate_recovery/summary_tables/second_best_hit_sanity.tsv), [summary_tables/non_informative_pairwise_check.tsv](/Users/sam/Documents/New project/validation/minimal_tier_b_candidate_recovery/summary_tables/non_informative_pairwise_check.tsv)
- run-level report: [pipeline_report.md](/Users/sam/Documents/New project/validation/minimal_tier_b_candidate_recovery/pipeline_report.md)

Pairwise note

- `summary_tables/non_informative_pairwise_check.tsv` is retained only as a record that the attempted Tier A vs Tier B comparison is not biologically informative in this archive because the Tier A comparator is a short mock transcript rather than a comparable recovered CDS.
