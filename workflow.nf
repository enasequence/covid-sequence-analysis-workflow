#!/usr/bin/env nextflow

params.OUTDIR = "results"
params.HUMAN_IDX = "/data/ref/human/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index"
params.SARS2_IDX = "/data/ref/sars2/index/NC_045512.2"
params.SARS2_FA = "/data/ref/sars2/fa/NC_045512.2.fa"
params.SARS2_FA_FAI = "/data/ref/sars2/fa/NC_045512.2.fa.fai"

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

    cpus 19
    memory '90 GB'
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

bt_indices_sars2 = Channel
    .fromPath("${params.SARS2_IDX}*", checkIfExists: true)

process align_reads_to_sars2_genome {
    cpus 19
    memory '90 GB'
    container 'alexeyebi/bowtie2_samtools'

    input:
    path fastq from bam_to_fastq_ch
    path sars2_fasta from params.SARS2_FA
    file indices from bt_indices_sars2.collect()
    val run_id from params.RUN_ID

    output:
    path "${run_id}.bam" into sars2_aligned_reads_ch

    script:
    index_base = indices[0].toString() - ~/.rev.\d.bt2?/ - ~/.\d.bt2?/
    """
    bowtie2 -p ${task.cpus} --no-mixed --no-discordant \
    --met-file ${run_id}_bowtie_nohuman_summary -x $index_base \
    -1 ${fastq[0]} -2 ${fastq[1]} | samtools view -bST ${sars2_fasta} | \
    samtools sort | samtools view -h -F 4 -b > ${run_id}.bam
    samtools index ${run_id}.bam
    """
}

process remove_duplicates {
    cpus 19
    memory '90 GB'
    container 'biocontainers/picard:v1.141_cv3'

    input:
    path bam from sars2_aligned_reads_ch
    val run_id from params.RUN_ID

    output:
    path "${run_id}_dep.bam" into remove_duplicates_ch, remove_duplicates2_ch

    script:
    """
    picard MarkDuplicates I=${bam} O=${run_id}_dep.bam REMOVE_DUPLICATES=true \
    M=${run_id}_marked_dup_metrics.txt
    """
}

process check_coverage {
    cpus 1
    memory '1 GB'
    container 'alexeyebi/bowtie2_samtools'

    input:
    path bam from remove_duplicates_ch
    val run_id from params.RUN_ID
    path sars2_fasta from params.SARS2_FA

    output:
    path "${run_id}.pileup" into check_coverage_ch

    script:
    """
    samtools mpileup -A -Q 30 -d 1000000 -f ${sars2_fasta} ${bam} > \
    ${run_id}.pileup
    """
}

process make_small_file_with_coverage {
    publishDir params.OUTDIR, mode:'copy'
    cpus 1
    memory '1 GB'
    container 'alexeyebi/bowtie2_samtools'

    input:
    path pileup from check_coverage_ch
    val run_id from params.RUN_ID

    output:
    path("${run_id}.coverage")

    script:
    """
    cat ${pileup} | awk '{print \$2,","\$3,","\$4}' > ${run_id}.coverage
    """
}

process generate_vcf {
    publishDir params.OUTDIR, mode:'copy'
    cpus 19
    memory '90 GB'
    container 'alexeyebi/bowtie2_samtools'

    input:
    path bam from remove_duplicates2_ch
    path sars2_fasta from params.SARS2_FA
    path sars2_fasta_fai from params.SARS2_FA_FAI
    val run_id from params.RUN_ID

    output:
    path "${run_id}.vcf.gz" into vcf_ch, vcf2_ch
    path("${run_id}.stat")

    script:
    """
    samtools index ${bam}
    lofreq call-parallel --pp-threads ${task.cpus} -f ${sars2_fasta} \
    -o ${run_id}.vcf ${bam}
    bgzip ${run_id}.vcf
    tabix ${run_id}.vcf.gz
    bcftools stats ${run_id}.vcf.gz > ${run_id}.stat
    """
}

process create_consensus_sequence {
    publishDir params.OUTDIR, mode:'copy'
    cpus 19
    memory '90 GB'
    container 'alexeyebi/bowtie2_samtools'

    input:
    path vcf from vcf_ch
    path sars2_fasta from params.SARS2_FA
    path sars2_fasta_fai from params.SARS2_FA_FAI
    val run_id from params.RUN_ID

    output:
    path("${run_id}.cons.fa")

    script:
    """
    bcftools filter -i "DP>50" ${vcf} -o ${run_id}.cfiltered.vcf
    bgzip ${run_id}.cfiltered.vcf
    tabix ${run_id}.cfiltered.vcf.gz
    bcftools filter -i "AF>0.5" ${run_id}.cfiltered.vcf.gz > \
    ${run_id}.cfiltered_freq.vcf
    bgzip -c ${run_id}.cfiltered_freq.vcf > ${run_id}.cfiltered_freq.vcf.gz
    bcftools index ${run_id}.cfiltered_freq.vcf.gz
    bcftools consensus -f ${sars2_fasta} ${run_id}.cfiltered_freq.vcf.gz > \
    ${run_id}.cons.fa
    sed -i "1s/.*/>${run_id}/" ${run_id}.cons.fa
    rm ${run_id}.cfiltered.vcf.gz
    rm ${run_id}.cfiltered.vcf.gz.tbi
    rm ${run_id}.cfiltered_freq.vcf
    rm ${run_id}.cfiltered_freq.vcf.gz.csi
    rm ${run_id}.cfiltered_freq.vcf.gz
    """
}

process filter_snv {
    cpus 1
    memory '1 GB'
    container 'alexeyebi/bowtie2_samtools'

    input:
    path vcf from vcf2_ch
    val run_id from params.RUN_ID

    output:
    path "${run_id}.filtered_freq.vcf" into filtered_freq_vcf_ch

    script:
    """
    bcftools filter -i "DP>50" ${vcf} -o ${run_id}.filtered.vcf
    bgzip ${run_id}.filtered.vcf
    tabix ${run_id}.filtered.vcf.gz
    bcftools filter -i "AF>0.1" ${run_id}.filtered.vcf.gz > \
    ${run_id}.filtered_freq.vcf
    """
}

process annotate_snps {
    publishDir params.OUTDIR, mode:'copy'
    cpus 19
    memory '90 GB'
    container 'alexeyebi/snpeff'

    input:
    path vcf from filtered_freq_vcf_ch
    val run_id from params.RUN_ID

    output:
    path("${run_id}.annot.n.filtered_freq.vcf")

    script:
    """
    cat ${vcf} | sed "s/^NC_045512.2/NC_045512/" > \
    ${run_id}.newchr.filtered_freq.vcf
    snpeff eff -v -s ${run_id}.snpEff_summary.html sars.cov.2 \
    ${run_id}.newchr.filtered_freq.vcf > ${run_id}.annot.n.filtered_freq.vcf
    """
}