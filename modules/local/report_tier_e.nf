process REPORT_TIER_E {
    label 'light'
    publishDir "${params.outdir}/tier_e", mode: params.publish_mode

    input:
    path species_plan
    path ortholog_targets

    output:
    path 'tier_e', emit: tier_dir
    path 'tier_e/ortholog_status.tsv', emit: status_tsv
    path 'tier_e/missing_species.tsv', emit: missing_tsv

    script:
    """
    report_missing_tier_e.py \
      --species-plan ${species_plan} \
      --ortholog-targets ${ortholog_targets} \
      --outdir tier_e
    """
}

