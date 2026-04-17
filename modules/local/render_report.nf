process RENDER_REPORT {
    label 'light'
    publishDir "${params.outdir}/reports", mode: params.publish_mode

    input:
    path status_matrix
    path status_long
    path sequence_dir
    path recovery_summary
    path run_summary

    output:
    path 'reports/pipeline_report.md', emit: report_md

    script:
    """
    render_run_report.py \
      --status-matrix ${status_matrix} \
      --status-long ${status_long} \
      --sequence-dir ${sequence_dir} \
      --recovery-summary ${recovery_summary} \
      --run-summary ${run_summary} \
      --outdir reports
    """
}

