process RUN_EXPRESSION {
    label 'medium'
    publishDir "${params.outdir}/expression", mode: params.publish_mode

    input:
    path species_manifest
    path reference_manifest

    output:
    path 'expression', emit: expression_dir

    script:
    """
    run_expression_scaffold.py \
      --species-manifest ${species_manifest} \
      --reference-manifest ${reference_manifest} \
      --outdir expression \
      --run-metadata "${params.run_metadata}" \
      --mode "${params.expression_mode}" \
      --threads ${task.cpus} \
      --ncbi-api-key "${params.ncbi_api_key}"
    """
}
