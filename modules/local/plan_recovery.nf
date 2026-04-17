process PLAN_RECOVERY {
    label 'light'
    publishDir "${params.outdir}/planning", mode: params.publish_mode

    input:
    path species_manifest
    path ortholog_targets
    path reference_manifest

    output:
    path 'planning/species_plan.tsv', emit: species_plan
    path 'planning/run_summary.json', emit: run_summary

    script:
    """
    mkdir -p planning

    plan_species_recovery.py \
      --species-manifest ${species_manifest} \
      --ortholog-targets ${ortholog_targets} \
      --reference-manifest ${reference_manifest} \
      --outdir planning
    """
}

