process RECOVER_TIER_C {
    label 'medium'
    publishDir "${params.outdir}/tier_c", mode: params.publish_mode

    input:
    path species_plan
    path ortholog_targets
    path reference_manifest

    output:
    path 'tier_c', emit: tier_dir
    path 'tier_c/ortholog_status.tsv', emit: status_tsv

    script:
    """
    recover_tier.py \
      --tier-label tier_c \
      --accepted-tiers C \
      --species-plan ${species_plan} \
      --ortholog-targets ${ortholog_targets} \
      --reference-manifest ${reference_manifest} \
      --outdir tier_c \
      --mode ${params.execution_mode}
    """
}

