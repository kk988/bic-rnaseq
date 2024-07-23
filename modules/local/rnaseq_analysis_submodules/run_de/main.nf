process RUN_DE {
    label "process_medium"

    container "docker://kk988/bic_rnaseq_modules:2.0.0"

    input:
    path counts_file
    path sample_groups
    path group_comparisons
    def conf_info

    output:
    path "AllRes*xlsx" , emit: all_gene_results
    path "ResDE*xlsx"  , emit: deseq2_results
    path "versions.yml", emit: versions

    def args  = task.ext.args  ?: ''

    script:
    """
    Rscript /rnaseq_analysis_modules/run_de.R \
    --htseq_file $counts_file \
    --key_file $sample_groups \
    --unfiltered_dir . \
    --filtered_dir . \
    --comparisons_file $group_comparisons \
    $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-deseq2: \$(Rscript -e "library(DESeq2); cat(as.character(packageVersion('DESeq2')))")
    END_VERSIONS
    """

}
