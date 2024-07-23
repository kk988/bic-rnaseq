process COUNTS_NORMALIZATION {
    label "process_medium"

    container "docker://kk988/bic_rnaseq_modules:2.0.0"

    input:
    path counts_file
    path sample_groups
    def conf_info

    output:
    path "*.xlsx"      , emit: scaled_counts
    path "versions.yml", emit: versions


    def norm_counts = counts_file.replace(".txt", "_norm_counts.xlsx")

    // grab the configuration info
    def min_tot_reads = conf_info.min_tot_reads ? "--min_tot_reads ${conf_info.min_tot_reads}" : ""

    script:
    """
    Rscript /rnaseq_analysis_modules/normalize_counts.R \
    --htseq_file $counts_file \
    --norm_counts_file $norm_counts \
    --key_file $sample_groups \
    $min_tot_reads

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-deseq2: \$(Rscript -e "library(DESeq2); cat(as.character(packageVersion('DESeq2')))")
    END_VERSIONS
    """
}
