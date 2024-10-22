process HTSEQ_MERGE_COUNTS {
    label "process_medium"

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    val(count_files)

    output:
    path "htseq.merged.counts.tsv", emit: merged_counts
    path "versions.yml", emit: versions

    script:
    """
    #!/bin/bash

    mkdir -p tmp_genes/

    echo "GeneID\tGeneSymbol" > tmp_genes/gene_ids.txt
    cut -f 1,2 \$(echo ${count_files[0]}) >> tmp_genes/gene_ids.txt

    for f in \$(echo ${count_files.join(' ')}); do
        sample_name=\$(basename \$f _htseq_counts.txt)

        echo \$sample_name >> tmp_genes/\$sample_name.counts.tsv
        cut -f 3 \$f >> tmp_genes/\$sample_name.counts.tsv
    done

    paste tmp_genes/gene_ids.txt tmp_genes/*.counts.tsv > htseq.merged.counts.tsv

    rm -r tmp_genes/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(echo \$(sed --version 2>&1) | sed 's/^.*GNU sed) //; s/ .*\$//')
        paste: \$(echo \$(paste --version 2>&1) | sed 's/^.*GNU coreutils) //; s/ .*\$//')
    END_VERSIONS
    """
}
