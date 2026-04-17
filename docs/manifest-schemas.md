# Manifest Schemas

## `species_manifest.tsv`

Required columns:

- `species_id`
- `scientific_name`

Supported optional columns:

- `gbif_species_key`
- `gbif_scientific_name`
- `gbif_authorship`
- `ena_tax_id`
- `assembly_accession`
- `assembly_level`
- `assembly_name`
- `ncbi_assembly_count`
- `ncbi_bioproject_count`
- `ncbi_sra_count`
- `ena_read_run_count`
- `ena_transcriptomic_run_count`
- `ena_genomic_run_count`
- `assembly_fasta`
- `annotation_gtf`
- `transcriptome_fasta`
- `rna_sra_accessions`
- `wgs_sra_accessions`
- `evidence_hint`
- `notes`

Classification rules:

- `A`: `assembly_fasta` and `annotation_gtf`
- `B`: `assembly_fasta` or `assembly_accession`
- `C`: `transcriptome_fasta` or `rna_sra_accessions`
- `D`: `wgs_sra_accessions`
- `E`: none of the above

`evidence_hint` overrides automatic classification when present.

The live discovery builder populates the discovery-oriented columns from GBIF, NCBI, and ENA, while leaving `assembly_fasta`, `annotation_gtf`, and `transcriptome_fasta` blank until local resources are staged.

## `ortholog_targets.tsv`

Required columns:

- `gene_id`
- `gene_symbol`
- `ref_species`
- `cds_length`

Useful future columns:

- `ref_transcript_id`
- `ref_protein_id`
- `reference_cds_fasta`
- `reference_protein_fasta`

## `reference_manifest.tsv`

Required columns:

- `reference_id`
- `scientific_name`
- `reference_role`

Supported optional columns:

- `assembly_fasta`
- `annotation_gtf`
- `protein_fasta`
- `assembly_accession`
- `assembly_level`
- `assembly_name`
- `reference_bed12`
- `reference_twobit`
- `query_chain`
- `busco_lineage`
- `notes`
