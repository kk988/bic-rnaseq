/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Print even though we only need
// DE specific configs

include { paramsSummaryLog } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Validate DE specific configs only

if (!params.sample_groups) {
    exit 1, "--sample_groups (tab delim file) required for differential expression!"
}

if (!params.group_comparisons) {
    exit 1, "--group_comparisons (tab delim file) required for differential expression!"
}

if (params.diff_expression_only && !params.counts_file) {
    exit 1, "--counts_file required for differential expression!"
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Take in default config for DE scripts
// Then check for config from params
ch_de_config        = Channel.fromPath("$projectDir/assets/de_config.yml", checkIfExists: true)
ch_de_custom_config = params.de_config ? Channel.fromPath(params.de_config) : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { DIFFERENTIAL_EXPRESSION as RUN_DE from './subworkflows/local/differential_expression' }

workflow DIFF_EXPRESSION {
    RUN_DE(
        ch_de_config,
        ch_de_custom_config.collect().ifEmpty([]),
        params.sample_groups,
        params.group_comparisons,
        params.counts_file
    )
}