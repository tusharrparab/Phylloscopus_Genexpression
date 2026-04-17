process RUN_ASR {
    label 'medium'
    publishDir "${params.outdir}/asr", mode: params.publish_mode

    input:
    path sequence_dir

    output:
    path 'asr', emit: asr_dir

    script:
    """
    run_asr_scaffold.py \
      --sequence-dir ${sequence_dir} \
      --outdir asr \
      --species-tree "${params.asr_species_tree}" \
      --mode "${params.asr_mode}" \
      --min-taxa ${params.asr_min_taxa} \
      --threads ${task.cpus} \
      --sequence-type "${params.asr_sequence_type}"
    """
}

