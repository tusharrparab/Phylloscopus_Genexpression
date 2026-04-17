process RECOVER_TIER_D {
    label 'medium'
    publishDir "${params.outdir}/tier_d", mode: params.publish_mode

    input:
    path species_plan
    path ortholog_targets
    path reference_manifest

    output:
    path 'tier_d', emit: tier_dir
    path 'tier_d/ortholog_status.tsv', emit: status_tsv

    script:
    """
    recover_tier.py \
      --tier-label tier_d \
      --accepted-tiers D \
      --species-plan ${species_plan} \
      --ortholog-targets ${ortholog_targets} \
      --reference-manifest ${reference_manifest} \
      --outdir tier_d \
      --mode ${params.execution_mode}
    """
}

