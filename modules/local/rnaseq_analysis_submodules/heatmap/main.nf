process HEATMAP {
    label "process_medium"

    container "docker://kk988/bic_rnaseq_modules:2.0.1"

    input:
    tuple(meta, de_results)
    path norm_counts
    path sample_groups
    def conf_info

    output:
    path "*.pdf"        , emit: heatmap_plot
    path "versions.yml" , emit: versions

    // grab the configuration info
    def conf_args = conf_info.annotate_samples == TRUE ? "--annotate_samples TRUE" : ""

    script:
    """
    Rscript /rnaseq_analysis_modules/make_de_heatmap.R \
    --norm_counts_file $norm_counts \
    --de_file $de_results \
    --sample_key $sample_groups \
    --conditions $meta.compa,$meta.compb \
    $conf_args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-deseq2: \$(Rscript -e "library(DESeq2); cat(as.character(packageVersion('DESeq2')))")
    END_VERSIONS
    """
}
