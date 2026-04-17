nextflow.enable.dsl = 2

include { VALIDATE_INPUTS } from './modules/local/validate_inputs'
include { PLAN_RECOVERY } from './modules/local/plan_recovery'
include { RECOVER_TIER_AB } from './modules/local/recover_tier_ab'
include { RECOVER_TIER_C } from './modules/local/recover_tier_c'
include { RECOVER_TIER_D } from './modules/local/recover_tier_d'
include { REPORT_TIER_E } from './modules/local/report_tier_e'
include { MERGE_RECOVERIES } from './modules/local/merge_recoveries'
include { RENDER_REPORT } from './modules/local/render_report'
include { RUN_ASR } from './modules/local/run_asr'
include { RUN_EXPRESSION } from './modules/local/run_expression'

def helpMessage() {
    return """
    Phylloscopus Comparative Ortholog Recovery Scaffold

    Required inputs:
      --species_manifest   TSV describing species and available sequence evidence
      --ortholog_targets   TSV defining a hypothesis-aware ortholog target panel
      --reference_manifest TSV listing reference assemblies and annotations

    Optional parameters:
      --outdir             Output directory for published results
      --execution_mode     stub or contract

    Example:
      ./nextflow run . -profile local
    """.stripIndent()
}

if (params.help) {
    log.info helpMessage()
    System.exit(0)
}

workflow {
    species_manifest_ch = Channel.fromPath(params.species_manifest, checkIfExists: true)
    ortholog_targets_ch = Channel.fromPath(params.ortholog_targets, checkIfExists: true)
    reference_manifest_ch = Channel.fromPath(params.reference_manifest, checkIfExists: true)

    validated = VALIDATE_INPUTS(species_manifest_ch, ortholog_targets_ch, reference_manifest_ch)

    planned = PLAN_RECOVERY(
        validated.species_manifest,
        validated.ortholog_targets,
        validated.reference_manifest
    )

    tier_ab = RECOVER_TIER_AB(
        planned.species_plan,
        validated.ortholog_targets,
        validated.reference_manifest
    )

    tier_c = RECOVER_TIER_C(
        planned.species_plan,
        validated.ortholog_targets,
        validated.reference_manifest
    )

    tier_d = RECOVER_TIER_D(
        planned.species_plan,
        validated.ortholog_targets,
        validated.reference_manifest
    )

    tier_e = REPORT_TIER_E(
        planned.species_plan,
        validated.ortholog_targets
    )

    merged = MERGE_RECOVERIES(
        tier_ab.tier_dir,
        tier_ab.status_tsv,
        tier_c.tier_dir,
        tier_c.status_tsv,
        tier_d.tier_dir,
        tier_d.status_tsv,
        tier_e.tier_dir,
        tier_e.status_tsv,
        planned.species_plan,
        validated.ortholog_targets
    )

    RENDER_REPORT(
        merged.status_matrix,
        merged.status_long,
        merged.sequence_dir,
        merged.summary_tsv,
        planned.run_summary
    )

    if (params.enable_asr) {
        RUN_ASR(merged.sequence_dir)
    }

    if (params.enable_expression) {
        RUN_EXPRESSION(validated.species_manifest, validated.reference_manifest)
    }
}
