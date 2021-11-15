#!/usr/bin/env nextflow

params.OUTDIR = "results"
params.SARS2_FA = "/data/ref/NC_045512.2.fa"
params.SARS2_FA_FAI = "/data/ref/NC_045512.2.fa.fai"

Channel
    .fromFilePairs(params.READS, checkIfExists:true)
    .into {read_pairs_ch; read_pairs2_ch}

process quality_control_pre {
    publishDir params.OUTDIR, mode:'copy'

    cpus 19
    memory '30 GB'
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
    cpus 19
    memory '30 GB'
    container 'davelabhub/trimmomatic:0.39--1'
    publishDir params.OUTDIR, mode:'copy'

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

    cpus 19
    memory '30 GB'
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

process align_reads_to_sars2_genome {
    publishDir params.OUTDIR, mode:'copy'
    cpus 19
    memory '30 GB'
    container 'alexeyebi/bowtie2_samtools'

    input:
    path trimmed_reads from trim_reads2_ch
    path sars2_fasta from params.SARS2_FA
    val run_id from params.RUN_ID

    output:
    path "${run_id}.bam" into sars2_aligned_reads_ch, sars2_aligned_reads_ch2
    path("${run_id}.bam")

    """
    bwa index ${sars2_fasta}
    bwa mem -t ${task.cpus} ${sars2_fasta} ${trimmed_reads[0]} ${trimmed_reads[2]} | samtools view -bF 4 - | samtools sort - > ${run_id}_paired.bam
    bwa mem -t ${task.cpus} ${sars2_fasta} <(cat ${trimmed_reads[1]} ${trimmed_reads[3]}) | samtools view -bF 4 - | samtools sort - > ${run_id}_unpaired.bam
    samtools merge ${run_id}.bam ${run_id}_paired.bam ${run_id}_unpaired.bam
    rm ${run_id}_paired.bam ${run_id}_unpaired.bam
    """
}

process check_coverage {
    publishDir params.OUTDIR, mode:'copy'
    cpus 1
    memory '10 GB'
    container 'alexeyebi/bowtie2_samtools'

    input:
    path bam from sars2_aligned_reads_ch
    val run_id from params.RUN_ID
    path sars2_fasta from params.SARS2_FA

    output:
    path "${run_id}.pileup" into check_coverage_ch
    path("${run_id}.pileup")

    script:
    """
    samtools mpileup -a -A -Q 30 -d 8000 -f ${sars2_fasta} ${bam} > \
    ${run_id}.pileup
    """
}

process make_small_file_with_coverage {
    publishDir params.OUTDIR, mode:'copy'
    cpus 1
    memory '10 GB'
    container 'alexeyebi/bowtie2_samtools'

    input:
    path pileup from check_coverage_ch
    val run_id from params.RUN_ID

    output:
    path("${run_id}.coverage")
    path "${run_id}.coverage" into coverage_ch

    script:
    """
    cat ${pileup} | awk '{print \$2,","\$3,","\$4}' > ${run_id}.coverage
    """
}

process generate_vcf {
    publishDir params.OUTDIR, mode:'copy'
    cpus 10
    memory '30 GB'
    container 'alexeyebi/bowtie2_samtools'

    input:
    path bam from sars2_aligned_reads_ch2
    path sars2_fasta from params.SARS2_FA
    path sars2_fasta_fai from params.SARS2_FA_FAI
    val run_id from params.RUN_ID

    output:
    path "${run_id}.vcf.gz" into vcf_ch, vcf_ch2
    path("${run_id}.stat")
    path("${run_id}_filtered.vcf.gz")

    script:
    """
    samtools index ${bam}
    lofreq indelqual --dindel ${bam} -f ${sars2_fasta} -o ${run_id}_fixed.bam
    samtools index ${run_id}_fixed.bam
    lofreq call-parallel --no-default-filter --call-indels --pp-threads ${task.cpus} -f ${sars2_fasta} -o ${run_id}.vcf ${run_id}_fixed.bam
    lofreq filter --af-min 0.25 -i ${run_id}.vcf -o ${run_id}_filtered.vcf
    bgzip ${run_id}.vcf
    bgzip ${run_id}_filtered.vcf
    tabix ${run_id}.vcf.gz
    bcftools stats ${run_id}.vcf.gz > ${run_id}.stat
    """
}

process annotate_snps {
    publishDir params.OUTDIR, mode:'copy'
    cpus 1
    memory '30 GB'
    container 'alexeyebi/snpeff'

    input:
    path vcf from vcf_ch
    val run_id from params.RUN_ID

    output:
    path("${run_id}.annot.vcf")

    script:
    """
    zcat ${vcf} | sed "s/^NC_045512.2/NC_045512/" > \
    ${run_id}.newchr.vcf
    java -Xmx4g -jar /data/tools/snpEff/snpEff.jar -q -no-downstream -no-upstream -noStats sars.cov.2 ${run_id}.newchr.vcf > ${run_id}.annot.vcf
    """
}

process create_consensus_sequence {
    publishDir params.OUTDIR, mode:'copy'
    cpus 1
    memory '30 GB'
    container 'alexeyebi/vcf_to_consensus'

    input:
    path vcf from vcf_ch2
    val run_id from params.RUN_ID
    path sars2_fasta from params.SARS2_FA
    path sars2_fasta_fai from params.SARS2_FA_FAI
    path coverage from coverage_ch

    output:
    path("${run_id}_consensus.fasta.gz")

    script:
    """
    tabix ${vcf}
    python /data/tools/vcf_to_consensus.py -dp 10 -af 0.25 -v ${vcf} -d ${coverage} -o ${run_id}_consensus.fasta -n ${run_id} -r ${sars2_fasta}
    bgzip ${run_id}_consensus.fasta
    """
}
