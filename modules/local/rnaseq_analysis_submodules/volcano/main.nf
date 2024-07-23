process VOLCANO {
    label "process_medium"

    container "docker://kk988/bic_rnaseq_modules:2.0.0"

    input:
    path all_results
    path group_comparisons
    val conf_info

    output:
    path "*.pdf"        , emit: volcano_plot
    path "versions.yml" , emit: versions

    // grab the configuration info

    script:
    """
    Rscript /rnaseq_analysis_modules/make_volcano_plot.R \
    --all_results_dir . \
    --comparisons_file $group_comparisons

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-deseq2: \$(Rscript -e "library(DESeq2); cat(as.character(packageVersion('DESeq2')))")
    END_VERSIONS
    """
}
