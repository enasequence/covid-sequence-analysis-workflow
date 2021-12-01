#!/usr/bin/env nextflow

params.OUTDIR = "gs://prj-int-dev-covid19-nf-gls/illumina-porting-workdir/results"
params.SARS2_FA = "gs://prj-int-dev-covid19-nf-gls/illumina-porting-workdir/data/ref/NC_045512.2.fa"
params.SARS2_FA_FAI = "gs://prj-int-dev-covid19-nf-gls/illumina-porting-workdir/data/ref/NC_045512.2.fa.fai"
params.sampels_idxs = "gs://prj-int-dev-covid19-nf-gls/illumina-porting-workdir/data/illumina.index.tsv"
params.STOREDIR = "gs://prj-int-dev-covid19-nf-gls/illumina-porting-workdir/storeDir"

nextflow.enable.dsl=2

workflow {
    data = Channel
            .fromPath(params.sampels_idxs)
            .splitCsv(header:true, sep:'\t')
            .map{ row-> tuple(row.run_accession, 'ftp://'+row.fastq_ftp.split(';')[0], 'ftp://'+row.fastq_ftp.split(';')[1]) }

    download_fastq(data)
    trimming_reads(download_fastq.out)
}

process download_fastq {
    storeDir params.STOREDIR

    // Use GLS default 1 CPU 1 GB and default quay.io/nextflow/bash
    // cpus 2
    // memory '1 GB'

    input:
    tuple val(sampleId), file(input_file_1), file(input_file_2)
    output:
    tuple val(sampleId), file("${sampleId}_1.fastq.gz"), file("${sampleId}_2.fastq.gz")

    script:
    """
    wget -t 0 -O ${sampleId}_1.fastq.gz \$(cat ${input_file_1})
    wget -t 0 -O ${sampleId}_2.fastq.gz \$(cat ${input_file_2})
    """
}

process trimming_reads {
    storeDir params.STOREDIR
    publishDir params.OUTDIR, mode:'copy'

    cpus 18
    memory '30 GB'
    container 'davelabhub/trimmomatic:0.39--1'

    input:
    tuple val(sampleId), path(reads)    //path(read1), path(read2)

    output:
    tuple val(sampleId), path("${sampleId}*.fq"), emit: trim_reads_ch, trim_reads2_ch
    path("${sampleId}_trim_summary")

    script:
    """
    trimmomatic PE ${reads} ${sampleId}_trim_1.fq ${sampleId}_trim_1_un.fq ${sampleId}_trim_2.fq ${sampleId}_trim_2_un.fq \
        -summary ${sampleId}_trim_summary -threads ${task.cpus} SLIDINGWINDOW:5:30 MINLEN:50
    """
}

process align_reads_to_sars2_genome {
    publishDir params.OUTDIR, mode:'copy'
    cpus 20
    memory '30 GB'
    container 'alexeyebi/bowtie2_samtools'

    input:
    path trimmed_reads from trim_reads2_ch
    path sars2_fasta from params.SARS2_FA
    val sampleId from params.RUN_ID

    output:
    path "${sampleId}.bam" into sars2_aligned_reads_ch, sars2_aligned_reads_ch2
    path("${sampleId}.bam")

    """
    bwa index ${sars2_fasta}
    bwa mem -t ${task.cpus} ${sars2_fasta} ${trimmed_reads[0]} ${trimmed_reads[2]} | samtools view -bF 4 - | samtools sort - > ${sampleId}_paired.bam
    bwa mem -t ${task.cpus} ${sars2_fasta} <(cat ${trimmed_reads[1]} ${trimmed_reads[3]}) | samtools view -bF 4 - | samtools sort - > ${sampleId}_unpaired.bam
    samtools merge ${sampleId}.bam ${sampleId}_paired.bam ${sampleId}_unpaired.bam
    rm ${sampleId}_paired.bam ${sampleId}_unpaired.bam
    """
}

process check_coverage {
    publishDir params.OUTDIR, mode:'copy'
    cpus 2
    memory '10 GB'
    container 'alexeyebi/bowtie2_samtools'

    input:
    path bam from sars2_aligned_reads_ch
    val sampleId from params.RUN_ID
    path sars2_fasta from params.SARS2_FA

    output:
    path "${sampleId}.pileup" into check_coverage_ch
    path("${sampleId}.pileup")

    script:
    """
    samtools mpileup -a -A -Q 30 -d 8000 -f ${sars2_fasta} ${bam} > \
    ${sampleId}.pileup
    """
}

process make_small_file_with_coverage {
    publishDir params.OUTDIR, mode:'copy'
    cpus 2
    memory '10 GB'
    container 'alexeyebi/bowtie2_samtools'

    input:
    path pileup from check_coverage_ch
    val sampleId from params.RUN_ID

    output:
    path("${sampleId}.coverage")
    path "${sampleId}.coverage" into coverage_ch

    script:
    """
    cat ${pileup} | awk '{print \$2,","\$3,","\$4}' > ${sampleId}.coverage
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
    val sampleId from params.RUN_ID

    output:
    path "${sampleId}.vcf.gz" into vcf_ch, vcf_ch2
    path("${sampleId}.stat")
    path("${sampleId}_filtered.vcf.gz")

    script:
    """
    samtools index ${bam}
    lofreq indelqual --dindel ${bam} -f ${sars2_fasta} -o ${sampleId}_fixed.bam
    samtools index ${sampleId}_fixed.bam
    lofreq call-parallel --no-default-filter --call-indels --pp-threads ${task.cpus} -f ${sars2_fasta} -o ${sampleId}.vcf ${sampleId}_fixed.bam
    lofreq filter --af-min 0.25 -i ${sampleId}.vcf -o ${sampleId}_filtered.vcf
    bgzip ${sampleId}.vcf
    bgzip ${sampleId}_filtered.vcf
    tabix ${sampleId}.vcf.gz
    bcftools stats ${sampleId}.vcf.gz > ${sampleId}.stat
    """
}

process annotate_snps {
    publishDir params.OUTDIR, mode:'copy'
    cpus 2
    memory '30 GB'
    container 'alexeyebi/snpeff'

    input:
    path(vcf) from vcf_ch
    val sampleId from params.RUN_ID

    output:
    path("${sampleId}.annot.vcf")

    script:
    """
    zcat ${vcf} | sed "s/^NC_045512.2/NC_045512/" > \
    ${sampleId}.newchr.vcf
    java -Xmx4g -jar /data/tools/snpEff/snpEff.jar -q -no-downstream -no-upstream -noStats sars.cov.2 ${sampleId}.newchr.vcf > ${sampleId}.annot.vcf
    """
}

process create_consensus_sequence {
    publishDir params.OUTDIR, mode:'copy'
    cpus 2
    memory '30 GB'
    container 'alexeyebi/vcf_to_consensus'

    input:
    path(vcf) from vcf_ch2
    val sampleId from params.RUN_ID
    path(sars2_fasta) from params.SARS2_FA
    path(sars2_fasta_fai) from params.SARS2_FA_FAI
    path(coverage) from coverage_ch

    output:
    path("${sampleId}_consensus.fasta.gz")

    script:
    """
    tabix ${vcf}
    python /data/tools/vcf_to_consensus.py -dp 10 -af 0.25 -v ${vcf} -d ${coverage} -o ${sampleId}_consensus.fasta -n ${sampleId} -r ${sars2_fasta}
    bgzip ${sampleId}_consensus.fasta
    """
}
