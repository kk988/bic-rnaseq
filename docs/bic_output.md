# nf-core/rnaseq: Output

## Introduction

This document describes the output produced by the pipeline that BIC will run. Most of the plots are taken from the MultiQC report generated from the [full-sized test dataset](https://github.com/nf-core/test-datasets/tree/rnaseq#full-test-dataset-origin) utilized in nf-core/rnaseq.


The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Preprocessing](#preprocessing)
  - [cat](#cat) - Merge re-sequenced FastQ files
  - [FastQC](#fastqc) - Raw read QC
- [Alignment and quantification](#alignment-and-quantification)
  - [STAR and HTSeq](#star-and-htseq) - alignment and gene level quantification
- [Alignment post-processing](#alignment-post-processing)
  - [SAMtools](#samtools) - Sort and index alignments
  - [picard MarkDuplicates](#picard-markduplicates) - Duplicate read marking
- [Differential Expression](#differential-expression)
- [Other steps](#other-steps)
  - [StringTie](#stringtie) - Transcript assembly and quantification
  - [BEDTools and bedGraphToBigWig](#bedtools-and-bedgraphtobigwig) - Create bigWig coverage files
- [Quality control](#quality-control)
  - [RSeQC](#rseqc) - Various RNA-seq QC metrics
  - [Qualimap](#qualimap) - Various RNA-seq QC metrics
  - [dupRadar](#dupradar) - Assessment of technical / biological read duplication
  - [DESeq2](#deseq2) - PCA plot and sample pairwise distance heatmap and dendrogram
  - [MultiQC](#multiqc) - Present QC for raw reads, alignment, read counting and sample similiarity

## Preprocessing

### cat

If multiple libraries/runs have been provided for the same sample in the input samplesheet (e.g. to increase sequencing depth) then these will be merged at the very beginning of the pipeline in order to have consistent sample naming throughout the pipeline. Please refer to the [usage documentation](https://nf-co.re/rnaseq/usage#samplesheet-input) to see how to specify these samples in the input samplesheet.

### FastQC

<details markdown="1">
<summary>Output files</summary>

- `fastqc/`
  - `*_fastqc.html`: FastQC report containing quality metrics.
  - `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

> **NB:** The FastQC plots in this directory are generated relative to the raw, input reads. They may contain adapter sequence and regions of low quality. To see how your reads look after adapter and quality trimming please refer to the FastQC reports in the `trimgalore/fastqc/` directory.
</details>
.

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

![MultiQC - FastQC sequence counts plot](images/mqc_fastqc_counts.png)

![MultiQC - FastQC mean quality scores plot](images/mqc_fastqc_quality.png)

![MultiQC - FastQC adapter content plot](images/mqc_fastqc_adapter.png)

## Alignment and quantification

### STAR and HTSeq

<details markdown="1">
<summary>Output files</summary>

- `star_htseq/log/`
  - `*.SJ.out.tab`: File containing filtered splice junctions detected after mapping the reads.
  - `*.Log.final.out`: STAR alignment report containing the mapping results summary.
  - `*.Log.out` and `*.Log.progress.out`: STAR log files containing detailed information about the run. Typically only useful for debugging purposes.
- `star_htseq/htseq`
  - `*_htseq_counts.txt`: individual count results from htseq
  - `htseq.merged.counts.tsv`: tab separated file containing counts for all samples 

> **NB:** Star alignment output is not saved by default. A downstream sorted alignment with duplicates marked (not removed) is stored by default. See [picard MarkDuplicates](#picard-markduplicates) for that information.


</details>

[STAR](https://github.com/alexdobin/STAR) is a read aligner designed for splice aware mapping typical of RNA sequencing data. STAR stands for *S*pliced *T*ranscripts *A*lignment to a *R*eference, and has been shown to have high accuracy and outperforms other aligners by more than a factor of 50 in mapping speed, but it is memory intensive. Using `--aligner star_htseq` is the default alignment and quantification option.

The STAR section of the MultiQC report shows a bar plot with alignment rates: good samples should have most reads as _Uniquely mapped_ and few _Unmapped_ reads.

![MultiQC - STAR alignment scores plot](images/mqc_star.png)

[`HTSeq`](https://htseq.readthedocs.io/en/master/) is utilized by running `htseq-count`, which counts reads within features. The output is a tab delimited file with counts for each feature in the gtf (as set by command line arguments) for the given sample. This is then merged into a matrix which has counts for each feature in the gtf for each sample.


## Alignment post-processing

The pipeline has been written in a way where all the files generated downstream of the alignment are placed in the same directory as specified by `--aligner` e.g. if `--aligner star_htseq` is specified then all the downstream results will be placed in the `star_htseq/` directory. This helps with organising the directory structure. It also means that results won't be overwritten when resuming the pipeline and can be used for benchmarking between alignment algorithms if required. BIC pipeline will only run star_htseq by default.

### SAMtools

<details markdown="1">
<summary>Output files</summary>

- `<ALIGNER>/samtools_stats/`
  - `<SAMPLE>.sorted.bam.flagstat`
  - `<SAMPLE>.sorted.bam.idxstats`
  - `<SAMPLE>.sorted.bam.stats` 

</details>

The original BAM files generated by the selected alignment algorithm are further processed with [SAMtools](http://samtools.sourceforge.net/) to sort them by coordinate, for indexing, as well as to generate read mapping statistics using the duplicate marked alignment files.

![MultiQC - SAMtools alignment scores plot](images/mqc_samtools_mapped.png)

![MultiQC - SAMtools mapped reads per contig plot](images/mqc_samtools_idxstats.png)

### picard MarkDuplicates

<details markdown="1">
<summary>Output files</summary>

- `<ALIGNER>/alignment/`
  - `<SAMPLE>.markdup.sorted.bam`: Coordinate sorted BAM file after duplicate marking. This is the final post-processed BAM file and so will be saved by default in the results directory.
  - `<SAMPLE>.markdup.sorted.bam.bai`: BAI index file for coordinate sorted BAM file after duplicate marking. This is the final post-processed BAM index file and so will be saved by default in the results directory.
- `<ALIGNER>/picard_metrics/`
  - `<SAMPLE>.markdup.sorted.MarkDuplicates.metrics.txt`: Metrics file from MarkDuplicates.

</details>

Unless you are using [UMIs](https://emea.illumina.com/science/sequencing-method-explorer/kits-and-arrays/umi.html) it is not possible to establish whether the fragments you have sequenced from your sample were derived via true biological duplication (i.e. sequencing independent template fragments) or as a result of PCR biases introduced during the library preparation. By default, the pipeline uses [picard MarkDuplicates](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates) to _mark_ the duplicate reads identified amongst the alignments to allow you to guage the overall level of duplication in your samples. However, for RNA-seq data it is not recommended to physically remove duplicate reads from the alignments (unless you are using UMIs) because you expect a significant level of true biological duplication that arises from the same fragments being sequenced from for example highly expressed genes. This step will be skipped automatically when using the `--with_umi` option or explicitly via the `--skip_markduplicates` parameter.

![MultiQC - Picard MarkDuplicates metrics plot](images/mqc_picard_markduplicates.png)

## Differential Expression

<details markdown="1">
<summary>Output files</summary>

- `differentialExpression_gene/`
</details>

Differential expression is an optional output, and will only be provided if key and comparison information was provided at the time of requesting analysis. Differential expression utilizes another pipeline called [bic-differentialabundance](https://github.com/cBio-MSKCC/bic-differentialabundance).

## Other steps

### StringTie

<details markdown="1">
<summary>Output files</summary>

- `<ALIGNER>/stringtie/`
  - `*.coverage.gtf`: GTF file containing transcripts that are fully covered by reads.
  - `*.transcripts.gtf`: GTF file containing all of the assembled transcipts from StringTie.
  - `*.gene_abundance.txt`: Text file containing gene aboundances and FPKM values.
- `<ALIGNER>/stringtie/<SAMPLE>.ballgown/`: Ballgown output directory.

</details>

[StringTie](https://ccb.jhu.edu/software/stringtie/) is a fast and highly efficient assembler of RNA-Seq alignments into potential transcripts. It uses a novel network flow algorithm as well as an optional de novo assembly step to assemble and quantitate full-length transcripts representing multiple splice variants for each gene locus. In order to identify differentially expressed genes between experiments, StringTie's output can be processed by specialized software like [Ballgown](https://github.com/alyssafrazee/ballgown), [Cuffdiff](http://cole-trapnell-lab.github.io/cufflinks/cuffdiff/index.html) or other programs ([DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html), [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html), etc.).

### BEDTools and bedGraphToBigWig

<details markdown="1">
<summary>Output files</summary>

- `<ALIGNER>/bigwig/`
  - `*.forward.bigWig`: bigWig coverage file relative to genes on the forward DNA strand.
  - `*.reverse.bigWig`: bigWig coverage file relative to genes on the reverse DNA strand.

</details>

The [bigWig](https://genome.ucsc.edu/goldenpath/help/bigWig.html) format is an indexed binary format useful for displaying dense, continuous data in Genome Browsers such as the [UCSC](https://genome.ucsc.edu/cgi-bin/hgTracks) and [IGV](http://software.broadinstitute.org/software/igv/). This mitigates the need to load the much larger BAM files for data visualisation purposes which will be slower and result in memory issues. The bigWig format is also supported by various bioinformatics software for downstream processing such as meta-profile plotting.

## Quality control

### RSeQC

[RSeQC](<(http://rseqc.sourceforge.net/)>) is a package of scripts designed to evaluate the quality of RNA-seq data. This pipeline runs several, but not all RSeQC scripts. By default the pipeline will run these following rseqc scripts: `bam_stat.py`, `inner_distance.py`, `infer_experiment.py`, `junction_annotation.py`, `junction_saturation.py`,`read_distribution.py` and `read_duplication.py`.

The majority of RSeQC scripts generate output files which can be plotted and summarised in the MultiQC report.

#### Infer experiment

<details markdown="1">
<summary>Output files</summary>

- `<ALIGNER>/rseqc/infer_experiment/`
  - `*.infer_experiment.txt`: File containing fraction of reads mapping to given strandedness configurations.

</details>

This script predicts the "strandedness" of the protocol (i.e. unstranded, sense or antisense) that was used to prepare the sample for sequencing by assessing the orientation in which aligned reads overlay gene features in the reference genome. The strandedness of each sample has to be provided to the pipeline in the input samplesheet (see [usage docs](https://nf-co.re/rnaseq/usage#samplesheet-input)). However, this information is not always available, especially for public datasets. As a result, additional features have been incorporated into this pipeline to auto-detect whether you have provided the correct information in the samplesheet, and if this is not the case then the affected libraries will be flagged in the table under 'Strandedness Checks' elsewhere in the report. If required, this will allow you to correct the input samplesheet and rerun the pipeline with the accurate strand information. Note, it is important to get this information right because it can affect the final results.

RSeQC documentation: [infer_experiment.py](http://rseqc.sourceforge.net/#infer-experiment-py)

![MultiQC - Strand check table](images/mqc_strand_check.png)

![MultiQC - RSeQC infer experiment plot](images/mqc_rseqc_inferexperiment.png)

#### Read distribution

<details markdown="1">
<summary>Output files</summary>

- `<ALIGNER>/rseqc/read_distribution/`
  - `*.read_distribution.txt`: File containing fraction of reads mapping to genome feature e.g. CDS exon, 5’UTR exon, 3’ UTR exon, Intron, Intergenic regions etc.

</details>

This tool calculates how mapped reads are distributed over genomic features. A good result for a standard RNA-seq experiments is generally to have as many exonic reads as possible (`CDS_Exons`). A large amount of intronic reads could be indicative of DNA contamination in your sample but may be expected for a total RNA preparation.

RSeQC documentation: [read_distribution.py](http://rseqc.sourceforge.net/#read-distribution-py)

![MultiQC - RSeQC read distribution plot](images/mqc_rseqc_readdistribution.png)

#### Junction annotation

<details markdown="1">
<summary>Output files</summary>

- `<ALIGNER>/rseqc/junction_annotation/bed/`
  - `*.junction.bed`: BED file containing splice junctions.
  - `*.junction.Interact.bed`: BED file containing interacting splice junctions.
- `<ALIGNER>/rseqc/junction_annotation/log/`
  - `*.junction_annotation.log`: Log file generated by the program.
- `<ALIGNER>/rseqc/junction_annotation/pdf/`
  - `*.splice_events.pdf`: PDF file containing splicing events plot.
  - `*.splice_junction.pdf`: PDF file containing splice junctions plot.
- `<ALIGNER>/rseqc/junction_annotation/rscript/`
  - `*.junction_plot.r`: R script used to generate pdf plots above.
- `<ALIGNER>/rseqc/junction_annotation/xls/`
  - `*.junction.xls`: Excel spreadsheet with junction information.

</details>

Junction annotation compares detected splice junctions to a reference gene model. Splicing annotation is performed in two levels: splice event level and splice junction level.

RSeQC documentation: [junction_annotation.py](http://rseqc.sourceforge.net/#junction-annotation-py)

![MultiQC - RSeQC junction annotation plot](images/mqc_rseqc_junctionannotation.png)

#### Inner distance

<details markdown="1">
<summary>Output files</summary>

- `<ALIGNER>/rseqc/inner_distance/pdf/`
  - `*.inner_distance_plot.pdf`: PDF file containing inner distance plot.
- `<ALIGNER>/rseqc/inner_distance/rscript/`
  - `*.inner_distance_plot.r`: R script used to generate pdf plot above.
- `<ALIGNER>/rseqc/inner_distance/txt/`
  - `*.inner_distance_freq.txt`: File containing frequency of insert sizes.
  - `*.inner_distance_mean.txt`: File containing mean, median and standard deviation of insert sizes.

</details>

The inner distance script tries to calculate the inner distance between two paired-end reads. It is the distance between the end of read 1 to the start of read 2, and it is sometimes confused with the insert size (see [this blog post](http://thegenomefactory.blogspot.com.au/2013/08/paired-end-read-confusion-library.html) for disambiguation):

This plot will not be generated for single-end data. Very short inner distances are often seen in old or degraded samples (_eg._ FFPE) and values can be negative if the reads overlap consistently.

RSeQC documentation: [inner_distance.py](http://rseqc.sourceforge.net/#inner-distance-py)

![MultiQC - RSeQC inner distance plot](images/mqc_rseqc_innerdistance.png)

#### Junction saturation

<details markdown="1">
<summary>Output files</summary>

- `<ALIGNER>/rseqc/junction_saturation/pdf/`
  - `*.junctionSaturation_plot.pdf`: PDF file containing junction saturation plot.
- `<ALIGNER>/rseqc/junction_saturation/rscript/`
  - `*.junctionSaturation_plot.r`: R script used to generate pdf plot above.

</details>

This script shows the number of splice sites detected within the data at various levels of subsampling. A sample that reaches a plateau before getting to 100% data indicates that all junctions in the library have been detected, and that further sequencing will not yield any more observations. A good sample should approach such a plateau of _Known junctions_, however, very deep sequencing is typically required to saturate all _Novel Junctions_ in a sample.

RSeQC documentation: [junction_saturation.py](http://rseqc.sourceforge.net/#junction-saturation-py)

![MultiQC - RSeQC junction saturation plot](images/mqc_rseqc_junctionsaturation.png)

#### Read duplication

<details markdown="1">
<summary>Output files</summary>

- `<ALIGNER>/rseqc/read_duplication/pdf/`
  - `*.DupRate_plot.pdf`: PDF file containing read duplication plot.
- `<ALIGNER>/rseqc/read_duplication/rscript/`
  - `*.DupRate_plot.r`: R script used to generate pdf plot above.
- `<ALIGNER>/rseqc/read_duplication/xls/`
  - `*.pos.DupRate.xls`: Read duplication rate determined from mapping position of read. First column is “occurrence” or duplication times, second column is number of uniquely mapped reads.
  - `*.seq.DupRate.xls`: Read duplication rate determined from sequence of read. First column is “occurrence” or duplication times, second column is number of uniquely mapped reads.

</details>

This plot shows the number of reads (y-axis) with a given number of exact duplicates (x-axis). Most reads in an RNA-seq library should have a low number of exact duplicates. Samples which have many reads with many duplicates (a large area under the curve) may be suffering excessive technical duplication.

RSeQC documentation: [read_duplication.py](http://rseqc.sourceforge.net/#read-duplication-py)

![MultiQC - RSeQC read duplication plot](images/mqc_rseqc_readduplication.png)

#### BAM stat

<details markdown="1">
<summary>Output files</summary>

- `<ALIGNER>/rseqc/bam_stat/`
  - `*.bam_stat.txt`: Mapping statistics for the BAM file.

</details>

This script gives numerous statistics about the aligned BAM files. A typical output looks as follows:

```txt
#Output (all numbers are read count)
#==================================================
Total records:                                 41465027
QC failed:                                     0
Optical/PCR duplicate:                         0
Non Primary Hits                               8720455
Unmapped reads:                                0

mapq < mapq_cut (non-unique):                  3127757
mapq >= mapq_cut (unique):                     29616815
Read-1:                                        14841738
Read-2:                                        14775077
Reads map to '+':                              14805391
Reads map to '-':                              14811424
Non-splice reads:                              25455360
Splice reads:                                  4161455
Reads mapped in proper pairs:                  21856264
Proper-paired reads map to different chrom:    7648
```

MultiQC plots each of these statistics in a dot plot. Each sample in the project is a dot - hover to see the sample highlighted across all fields.

RSeQC documentation: [bam_stat.py](http://rseqc.sourceforge.net/#bam-stat-py)

#### TIN

<details markdown="1">
<summary>Output files</summary>

- `<ALIGNER>/rseqc/tin/`
  - `*.summary.txt`: File containing TIN results summary.
  - `*.tin.xls`: XLS file containing TIN results.

</details>

This script is designed to evaluate RNA integrity at the transcript level. TIN (transcript integrity number) is named in analogous to RIN (RNA integrity number). RIN (RNA integrity number) is the most widely used metric to evaluate RNA integrity at sample (or transcriptome) level. It is a very useful preventive measure to ensure good RNA quality and robust, reproducible RNA sequencing. This process isn't run by default - please see [this issue](https://github.com/nf-core/rnaseq/issues/769).

RSeQC documentation: [tin.py](http://rseqc.sourceforge.net/#tin-py)

### Qualimap

<details markdown="1">
<summary>Output files</summary>

- `<ALIGNER>/qualimap/<SAMPLE>/`
  - `qualimapReport.html`: Qualimap HTML report that can be viewed in a web browser.
  - `rnaseq_qc_results.txt`: Textual results output.
  - `images_qualimapReport/`: Images required for the HTML report.
  - `raw_data_qualimapReport/`: Raw data required for the HTML report.
  - `css/`: CSS files required for the HTML report.

</details>

[Qualimap](http://qualimap.bioinfo.cipf.es/) is a platform-independent application written in Java and R that provides both a Graphical User Interface (GUI) and a command-line interface to facilitate the quality control of alignment sequencing data. Shortly, Qualimap:

- Examines sequencing alignment data according to the features of the mapped reads and their genomic properties.
- Provides an overall view of the data that helps to to the detect biases in the sequencing and/or mapping of the data and eases decision-making for further analysis.

The [Qualimap RNA-seq QC module](http://qualimap.bioinfo.cipf.es/doc_html/analysis.html#rna-seq-qc) is used within this pipeline to assess the overall mapping and coverage relative to gene features.

![MultiQC - Qualimap gene coverage plot](images/mqc_qualimap_coverage.png)

![MultiQC - Qualimap genomic origin plot](images/mqc_qualimap_features.png)

### dupRadar

<details markdown="1">
<summary>Output files</summary>

- `<ALIGNER>/dupradar/box_plot/`
  - `*_duprateExpBoxplot.pdf`: PDF file containing box plot for duplicate rate relative to mean expression.
- `<ALIGNER>/dupradar/gene_data/`
  - `*_dupMatrix.txt`: Text file containing duplicate metrics per gene.
- `<ALIGNER>/dupradar/histogram/`
  - `*_expressionHist.pdf`: PDF file containing histogram of reads per kilobase values per gene.
- `<ALIGNER>/dupradar/intercepts_slope/`
  - `*_intercept_slope.txt`: Text file containing intercept slope values.
- `<ALIGNER>/dupradar/scatter_plot/`
  - `*_duprateExpDens.pdf`: PDF file containing typical dupRadar 2D density scatter plot.

See [dupRadar docs](https://www.bioconductor.org/packages/devel/bioc/vignettes/dupRadar/inst/doc/dupRadar.html) for further information regarding the content of these files.

</details>

[dupRadar](https://www.bioconductor.org/packages/release/bioc/html/dupRadar.html) is a Bioconductor library written in the R programming language. It generates various QC metrics and plots that relate duplication rate with gene expression levels in order to identify experiments with high technical duplication. A good sample with little technical duplication will only show high numbers of duplicates for highly expressed genes. Samples with technical duplication will have high duplication for all genes, irrespective of transcription level.

![dupRadar - Example good and bad experiment plot](images/dupradar_example_plot.png)

> _Credit: [dupRadar documentation](https://www.bioconductor.org/packages/devel/bioc/vignettes/dupRadar/inst/doc/dupRadar.html)_

### DESeq2

<details markdown="1">
<summary>Output files</summary>

- `<ALIGNER/PSEUDOALIGNER>/deseq2_qc/`
  - `*.plots.pdf`: File containing PCA and hierarchical clustering plots.
  - `*.dds.RData`: File containing R `DESeqDataSet` object generated
    by DESeq2, with either an rlog or vst `assay` storing the
    variance-stabilised data.
  - `*.rds`: Alternative version of the RData file suitable for
    `readRDS` to give user control of the eventual object name.
  - `*pca.vals.txt`: Matrix of values for the first 2 principal components.
  - `*sample.dists.txt`: Sample distance matrix.
  - `R_sessionInfo.log`: File containing information about R, the OS and attached or loaded packages.
- `<ALIGNER/PSEUDOALIGNER>/deseq2_qc/size_factors/`
  - `*.txt`, `*.RData`: Files containing DESeq2 sizeFactors per sample.

</details>

[DESeq2](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) is one of the most commonly used software packages to perform differential expression analysis for RNA-seq datasets.

**This pipeline uses a standardised DESeq2 analysis script to get an idea of the reproducibility across samples within the experiment. Please note that this will not suit every experimental design, and if there are other problems with the experiment then it may not work as well as expected.**

The script included in the pipeline uses DESeq2 to normalise read counts across all of the provided samples in order to create a PCA plot and a clustered heatmap showing pairwise Euclidean distances between the samples in the experiment. These help to show the similarity between groups of samples and can reveal batch effects and other potential issues with the experiment.

By default, the pipeline uses the `vst` transformation which is more suited to larger experiments. You can set the parameter `--deseq2_vst false` if you wish to use the DESeq2 native `rlog` option. See [DESeq2 docs](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#data-transformations-and-visualization) for a more detailed explanation.

Both types of transformation are performed blind, i.e. using across-all-samples variability, without using any prior information on experimental groups (equivalent to using an intercept-only design), as recommended by the [DESeq2 docs](https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#blind-dispersion-estimation).

The PCA plots are generated based alternately on the top five hundred most variable genes, or all genes. The former is the conventional approach that is more likely to pick up strong effects (ie the biological signal) and the latter, when different, is picking up a weaker but consistent effect that is synchronised across many transcripts. We project both of these onto the first two PCs (shown in the top row of the figure below), which is the best two dimensional representation of the variation between samples.

We also explore higher components in terms of experimental factors inferred from sample names. If your sample naming convention follows a strict policy of using underscores to delimit values of experimental factors (for example `WT_UNTREATED_REP1`) and all names have the same number of underscores (so excluding `WT_TREATED_10ml_REP1` from being compatible with the previous label), then any of these factors that are informative (ie label some but not all samples the same) then we individually plot upto the first five PCs, per experimental level, for each of the experimental factors.

The plot on the left hand side shows the standard PC plot - notice the variable number of underscores, meaning that the central plot would not be produced: here we have changed the underscore that is hyphenating the treatment to a '-' character. This allows the central plot to be generated, and we can see that replicate (the 2nd part of the sample name) seems to be affecting the 3rd principal component, but the treatment factor is affecting the more important first two components. The right-most plot shows all pairwise euclidean distances between the samples.

<p align="center"><img src="images/deseq2_qc_plots.png" alt="DESeq2 PCA plots"></p>

![MultiQC - DESeq2 PCA plot](images/mqc_deseq2_pca.png)

<p align="center"><img src="images/mqc_deseq2_clustering.png" alt="MultiQC - DESeq2 sample similarity plot" width="600"></p>

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/<ALIGNER>/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools i.e. FastQC, Cutadapt, SortMeRNA, STAR, RSEM, HISAT2, Salmon, SAMtools, Picard, RSeQC, Qualimap, Preseq and featureCounts. Additionally, various custom content has been added to the report to assess the output of dupRadar, DESeq2 and featureCounts biotypes, and to highlight samples failing a mimimum mapping threshold or those that failed to match the strand-specificity provided in the input samplesheet. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.
