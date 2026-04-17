process VALIDATE_INPUTS {
    label 'light'
    publishDir "${params.outdir}/validation", mode: params.publish_mode

    input:
    path species_manifest
    path ortholog_targets
    path reference_manifest

    output:
    path 'validated/species_manifest.tsv', emit: species_manifest
    path 'validated/ortholog_targets.tsv', emit: ortholog_targets
    path 'validated/reference_manifest.tsv', emit: reference_manifest
    path 'validated/validation_summary.tsv', emit: summary

    script:
    """
    mkdir -p validated

    validate_inputs.py \
      --species-manifest ${species_manifest} \
      --ortholog-targets ${ortholog_targets} \
      --reference-manifest ${reference_manifest} \
      --outdir validated
    """
}

