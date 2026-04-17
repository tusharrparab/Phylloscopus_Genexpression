process RECOVER_TIER_AB {
    label 'medium'
    publishDir "${params.outdir}/tier_ab", mode: params.publish_mode

    input:
    path species_plan
    path ortholog_targets
    path reference_manifest

    output:
    path 'tier_ab', emit: tier_dir
    path 'tier_ab/ortholog_status.tsv', emit: status_tsv

    script:
    """
    recover_tier.py \
      --tier-label tier_ab \
      --accepted-tiers A,B \
      --species-plan ${species_plan} \
      --ortholog-targets ${ortholog_targets} \
      --reference-manifest ${reference_manifest} \
      --outdir tier_ab \
      --mode ${params.execution_mode}
    """
}

