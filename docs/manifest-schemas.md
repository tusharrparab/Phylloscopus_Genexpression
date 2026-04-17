# Manifest Schemas

## `species_manifest.tsv`

Required columns:

- `species_id`
- `scientific_name`

Supported evidence and planning columns:

- `assembly_accession`
- `assembly_level`
- `assembly_name`
- `assembly_fasta`
- `annotation_gtf`
- `transcriptome_fasta`
- `rna_sra_accessions`
- `wgs_sra_accessions`
- `evidence_hint`
- `evidence_confidence`
- `data_provenance`
- `analysis_suitability`
- `notes`

Discovery-oriented optional columns:

- `gbif_species_key`
- `gbif_scientific_name`
- `gbif_authorship`
- `ena_tax_id`
- `ncbi_assembly_count`
- `ncbi_bioproject_count`
- `ncbi_sra_count`
- `ena_read_run_count`
- `ena_transcriptomic_run_count`
- `ena_genomic_run_count`

Tier rules:

- `A`: `assembly_fasta` and `annotation_gtf`
- `B`: `assembly_fasta` or `assembly_accession`
- `C`: `transcriptome_fasta` or `rna_sra_accessions`
- `D`: `wgs_sra_accessions`
- `E`: none of the above

`evidence_hint` overrides automatic classification when present.

Recommended interpretation of `evidence_confidence`:

- `high`: local assembly and annotation are already staged
- `medium`: assembly-backed but annotation or projection assets remain incomplete
- `low`: evidence exists but orthology will require much more work
- `none`: no current route to recovery

## `ortholog_targets.tsv`

Required columns:

- `gene_id`
- `gene_symbol`
- `category`
- `orthology_basis`
- `copy_number_expectation`
- `ref_species`
- `cds_length`
- `rationale`

Supported optional columns:

- `analysis_notes`
- `ref_transcript_id`
- `ref_protein_id`
- `reference_cds_fasta`
- `reference_protein_fasta`

Allowed `category` values:

- `phylogenetic_backbone`
- `migration_candidate`
- `vocalization_neural_candidate`
- `hypoxia_elevation_candidate`
- `housekeeping_control`

Allowed `orthology_basis` values:

- `historical_single_locus_marker`
- `single_copy_preferred`
- `candidate_gene_requires_validation`
- `expression_control_not_for_phylogeny`

Allowed `copy_number_expectation` values:

- `single_copy_preferred`
- `screen_for_paralogs`
- `control_not_interpreted_as_phylogeny`

Interpretation:

- backbone markers are the most suitable starting point for comparative sequence inference
- candidate loci are biologically motivated but require stronger orthology checks
- housekeeping controls are not primary phylogenetic markers

## `reference_manifest.tsv`

Required columns:

- `reference_id`
- `scientific_name`
- `reference_role`

Supported optional columns:

- `reference_quality`
- `data_provenance`
- `assembly_fasta`
- `annotation_gtf`
- `protein_fasta`
- `transcript_fasta`
- `assembly_accession`
- `assembly_level`
- `assembly_name`
- `reference_bed12`
- `reference_twobit`
- `query_chain`
- `busco_lineage`
- `analysis_notes`
- `notes`

`reference_role` must contain exactly one `primary` row.
