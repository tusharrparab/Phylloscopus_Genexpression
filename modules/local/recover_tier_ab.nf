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
    run_tier_ab_scaffold.py \
      --species-plan ${species_plan} \
      --ortholog-targets ${ortholog_targets} \
      --reference-manifest ${reference_manifest} \
      --outdir tier_ab \
      --mode "${params.tier_ab_mode}" \
      --cpus ${task.cpus} \
      --ncbi-api-key "${params.ncbi_api_key}" \
      --busco-lineage "${params.busco_lineage}"
    """
}
