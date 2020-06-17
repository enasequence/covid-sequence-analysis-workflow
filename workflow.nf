#!/usr/bin/env nextflow

params.READS = "/Users/alexey/ebi_projects/covid-sequence-analysis-workflow/SRR11092064_{1,2}.fastq.gz"
params.OUTDIR = "results"

Channel
    .fromFilePairs(params.READS, checkIfExists:true)
    .into {read_pairs_ch; read_pairs2_ch}

process quality_control_pre {
    publishDir params.OUTDIR, mode:'copy'

    cpus 1
    memory '1 GB'
    container 'biocontainers/fastqc:v0.11.8dfsg-2-deb_cv1'

    input:
    tuple val(run_id), path(reads) from read_pairs_ch

    output:
    path("${run_id}_*_fastqc.html")

    script:
    """
    fastqc -t ${task.cpus} -q ${reads[0]} ${reads[1]}
    """
}

process trimming_reads {
    publishDir params.OUTDIR, mode:'copy'

    cpus 1
    memory '1 GB'
    container 'davelabhub/trimmomatic:0.39--1'

    input:
    tuple val(run_id), path(reads) from read_pairs2_ch

    output:
    path "${run_id}*.fq" into trim_reads_ch
    path("${run_id}_trim_summary")

    script:
    """
    trimmomatic PE ${reads} ${run_id}_trim_1.fq \
    ${run_id}_trim_1_un.fq ${run_id}_trim_2.fq ${run_id}_trim_2_un.fq \
    -summary ${run_id}_trim_summary -threads ${task.cpus} \
    SLIDINGWINDOW:5:30 MINLEN:50
    """
}

process quality_control_post {
    publishDir params.OUTDIR, mode:'copy'

    cpus 1
    memory '1 GB'
    container 'biocontainers/fastqc:v0.11.8dfsg-2-deb_cv1'

    input:
    path trimmed_reads from trim_reads_ch

    output:
    path("*.html")

    script:
    """
    fastqc -t ${task.cpus} -q ${trimmed_reads}
    """
}