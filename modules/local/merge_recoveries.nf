process MERGE_RECOVERIES {
    label 'light'
    publishDir "${params.outdir}/merged", mode: params.publish_mode

    input:
    path tier_ab_dir
    path tier_ab_status
    path tier_c_dir
    path tier_c_status
    path tier_d_dir
    path tier_d_status
    path tier_e_dir
    path tier_e_status
    path species_plan
    path ortholog_targets

    output:
    path 'merged/ortholog_status_long.tsv', emit: status_long
    path 'merged/ortholog_status_matrix.tsv', emit: status_matrix
    path 'merged/ortholog_sequences', emit: sequence_dir
    path 'merged/recovery_summary.tsv', emit: summary_tsv

    script:
    """
    merge_recoveries.py \
      --species-plan ${species_plan} \
      --ortholog-targets ${ortholog_targets} \
      --tier-dir ${tier_ab_dir} \
      --tier-status ${tier_ab_status} \
      --tier-dir ${tier_c_dir} \
      --tier-status ${tier_c_status} \
      --tier-dir ${tier_d_dir} \
      --tier-status ${tier_d_status} \
      --tier-dir ${tier_e_dir} \
      --tier-status ${tier_e_status} \
      --outdir merged
    """
}

