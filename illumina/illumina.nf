#!/usr/bin/env nextflow

params.INDEX = "gs://prj-int-dev-covid19-nf-gls/prepro/illumina.index.tsv"
params.STOREDIR = "gs://prj-int-dev-covid19-nf-gls/prepro/storeDir"
params.OUTDIR = "gs://prj-int-dev-covid19-nf-gls/prepro/results"

params.STUDY = 'PRJEB45555'
params.TEST_SUBMISSION = 'true'

nextflow.enable.dsl = 2

process map_to_reference {
    storeDir params.STOREDIR

    cpus 8
    memory '8 GB'
    container 'quay.io/yanisasunt/ena-sars-cov2-illumina:1.0'

    input:
    tuple val(run_accession), val(sample_accession), file(input_file_1), file(input_file_2) //from samples_ch
    path(sars2_fasta)
    path(sars2_fasta_fai)
    path(projects_accounts_csv)
    val(study_accession)

    output:
    val(run_accession)
    val(sample_accession)
    file("${run_accession}.bam")
    file("${run_accession}.coverage.gz")
    file("${run_accession}.annot.vcf.gz")
    file("${run_accession}_filtered.vcf.gz")
    file("${run_accession}_consensus.fasta.gz")

    script:
    """
    echo "hello test"
    mamba run download.sh ${run_accession} ${projects_accounts_csv} ${input_file_1} ${input_file_2} ${study_accession} 2>>log.err
    echo "mamba test"
    try:
        mamba run assemblies.sh ${run_accession} ${sars2_fasta} ${task.cpus} 2>>log.err
    catch Exception e:
        mamba run annotation.sh ${run_accession} ${sars2_fasta} ${task.cpus} 2>>log.err
    echo "try test"
    mamba run push_log.sh ${run_accession} log.err
    """
}

include { ena_analysis_submit } from '../nextflow-lib/ena-analysis-submitter.nf'
workflow {
    data = Channel
            .fromPath(params.INDEX)
            .splitCsv(header: true, sep: '\t')
            .map { row -> tuple(row.run_accession, row.sample_accession, 'ftp://' + row.fastq_ftp.split(';').takeRight(2)[0], 'ftp://' + row.fastq_ftp.split(';').takeRight(2)[1]) }
//            .map { row -> tuple(row.run_accession, row.sample_accession, 'ftp://' + row.fastq_ftp.split(';')[0], 'ftp://' + row.fastq_ftp.split(';')[1]) }

    map_to_reference(data, params.SARS2_FA, params.SARS2_FA_FAI, params.SECRETS, params.STUDY)
    ena_analysis_submit(map_to_reference.out, params.SECRETS, params.STUDY, params.TEST_SUBMISSION, params.CONFIG_YAML)
}