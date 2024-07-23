
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
    def yaml = new org.yaml.snakeyaml.Yaml();
    def de_run_conf = yaml.load(default_config);
    def conf2 = custom_config != null ? yaml.load(custom_config) : [:];
    de_run_conf.putAll(conf2);

    // Run steps that are reqested in config
    ch_scaled_counts = Channel.empty()

    if de_run_conf.normalization.run == true {
        COUNTS_NORMALIZATION {
            counts_file,
            de_run_conf.normalization
        }
        ch_scaled_counts = COUNTS_NORMALIZATION.out.scaled_counts
    } else {
        ch_scaled_counts = counts_file
    }

    ch_all_results = Channel.empty()
    ch_de_results = Channel.empty()
    if de_run_conf.run_de.run == true {
        RUN_DE {
            counts_file,
            sample_groups,
            group_comparisons,
            de_run_conf.run_de
        }
        ch_all_results = RUN_DE.out.all_gene_results
        ch_de_results = RUN_DE.out.deseq2_results
    }

    if de_run_conf.clustering.run == true {
        CLUSTERING {
            counts_file,
            de_run_conf.clustering
        }
    }

    //
    // Prepare metadata for heatmap plots
    //
    ch_de_results.map{
        filename ->
            def comp_matcher = filename =~ /"ALLResDESeq_"(.+)_vs_(.+).xlsx/
            def meta = {
                "comparison": comp_matcher[0][1] + "_vs_" + comp_matcher[0][2],
                "compb": comp_matcher[0][1],
                "compa": comp_matcher[0][2]
            }
            return [meta, filename]
    }.set{ ch_de_results_meta }

    if de_run_conf.heatmap.run == true {
        HEATMAP {
            ch_de_results_meta,
            ch_scaled_counts,
            sample_groups,
            de_run_conf.heatmap
        }
    }

    if de_run_conf.volcano.run == true {
        VOLCANO {
            ch_all_results.first(),
            group_comparisons,
            de_run_conf.volcano
        }
    }

}
