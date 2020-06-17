#!/usr/bin/env nextflow

params.READS = "/sars-cov2-sequence-analysis/SRR11092064_{1,2}.fastq.gz"
params.OUTDIR = "results"
params.HUMAN_IDX = "/sars-cov2-sequence-analysis/ref/human/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index"
params.RUN_ID = "SRR11092064"

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
    path "${run_id}*.fq" into trim_reads_ch, trim_reads2_ch
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

bt_indices = Channel
    .fromPath("${params.HUMAN_IDX}*", checkIfExists: true)

process align_reads {
    cpus 19
    memory '90 GB'
    container 'alexeyebi/bowtie2_samtools'

    input:
    path trimmed_reads from trim_reads2_ch
    file indices from bt_indices.collect()
    val run_id from params.RUN_ID

    output:
    path "${run_id}_nohuman.bam" into aligned_reads_ch

    script:
    index_base = indices[0].toString() - ~/.rev.\d.bt2?/ - ~/.\d.bt2?/
    """
    bowtie2 --very-sensitive-local -p ${task.cpus} \
    -x $index_base --met-file ${run_id}_bowtie_human_summary \
    -1 ${trimmed_reads[0]} -2 ${trimmed_reads[2]} \
    -U ${trimmed_reads[1]},${trimmed_reads[3]} | \
    samtools view -Sb -f 4 > ${run_id}_nohuman.bam
    """
}

process convert_bam_to_fastq {
    cpus 1
    memory '1 GB'
    container 'alexeyebi/bowtie2_samtools'

    input:
    path bam from aligned_reads_ch
    val run_id from params.RUN_ID

    output:
    path "${run_id}_nohuman*fq" into bam_to_fastq_ch

    script:
    """
    samtools bam2fq -1 ${run_id}_nohuman_1.fq -2 ${run_id}_nohuman_2.fq \
    -s ${run_id}_nohuman_s.fq ${bam} > ${run_id}_nohuman_3.fq
    """
}

process align_reads_to_sars2_genome {
    cpus 19
    memory '90 GB'
    container 'alexeyebi/bowtie2_samtools'

    input:

    output:

    script:
    """
    bowtie2 -p $v_threads --no-mixed --no-discordant --met-file $s \
    -x $v_sars2_idx -1 $f -2 $r | samtools view -bST $v_sars2_fa | \
    samtools sort | samtools view -h -F 4 -b > $bam
    samtools index $bam
    """
}