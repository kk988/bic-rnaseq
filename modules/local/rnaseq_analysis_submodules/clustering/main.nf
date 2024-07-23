process CLUSTERING {
    label "process_medium"

    container "docker://kk988/bic_rnaseq_modules:2.0.0"

    input:
    path counts_file
    path sample_groups
    def conf_info

    output:
    path "*_MDS.pdf"   , optional:true, emit: mds_plot
    path "*_PCA.pdf"   , optional:true, emit: pca_plot
    path "*_hclust.pdf", optional:true, emit: hclust_plot
    path "versions.yml", emit: versions

    // grab the configuration info
    def pca = conf_info.pca.run && conf_info.pca.run == true ? "PCA" : ""
    def mds = conf_info.pca.run && conf_info.pca.run == true ? "MDS" : ""
    def hclust = conf_info.hclust.run && conf_info.hclust.run == true ? "hclust" : ""
    def method = "${pca} ${mds} ${hclust}".trim()
    def pca_label_samps = conf_info.pca.label_samples ? "--pca_label_samps ${conf_info.pca_label_samples}" : ""
    def mds_label_samps = conf_info.mds.label_samples ? "--mds_label_samps ${conf_info.mds_label_samples}" : ""
    def min_tot_reads = conf_info.min_total_reads ? "--min_total_reads ${conf_info.min_total_reads}" : ""

    script:
    """
    Rscript /rnaseq_analysis_modules/clustering.R --htseq_file $counts_file --key_file $sample_groups --method $method $pca_label_samps $mds_label_samps $min_tot_reads

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-deseq2: \$(Rscript -e "library(DESeq2); cat(as.character(packageVersion('DESeq2')))")
    END_VERSIONS
    """
}
