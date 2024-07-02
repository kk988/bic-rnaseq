 
include { CLUSTERING           } from '../../../modules/local/rnaseq_analysis_modules/clustering'
include { COUNTS_NORMALIZATION } from '../../../modules/local/rnaseq_analysis_modules/counts_normalization'
include { HEATMAP              } from '../../../modules/local/rnaseq_analysis_modules/heatmap'
include { RUN_DE               } from '../../../modules/local/rnaseq_analysis_modules/run_de'
include { VOLCANO              } from '../../../modules/local/rnaseq_analysis_modules/volcano'

workflow DIFFERENTIAL_EXPRESSION {
    take:
    default_config    // channel: /path/to/config.yml
    custom_config     // channel: /path/to/config.yml
    sample_groups     // channel: /path/to/sample/groups.txt
    group_comparisons // channel: /path/to/group/comparisons.txt
    counts_file       // channel: /path/to/counts/file.tsv

    main:

    //
    // Load default config
    // If a custom config is present overwrite the default values with custom vals
    // 

    
}