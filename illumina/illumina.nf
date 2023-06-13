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
    val(gcloud_project_id)
    path(gcloud_key_file)

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
    gcloud auth activate-service-account --key-file="${gcloud_key_file}"
    gcloud config set project ${gcloud_project_id}
    download.sh ${run_accession} ${projects_accounts_csv} ${input_file_1} ${input_file_2} ${study_accession} 2>&1 >/dev/null | tee -a log.err
    assemblies.sh ${run_accession} ${sars2_fasta} ${task.cpus} 2>&1 >/dev/null | tee -a log.err
    annotation.sh ${run_accession} ${sars2_fasta} ${task.cpus} 2>&1 >/dev/null | tee -a log.err
	push_log.sh ${run_accession} log.err ${PROFILE} ${gcloud_project_id} 
    """
}

include { ena_analysis_submit } from '../nextflow-lib/ena-analysis-submitter.nf'
workflow {
    data = Channel
            .fromPath(params.INDEX)
            .splitCsv(header: true, sep: '\t')
            .map { row -> tuple(row.run_accession, row.sample_accession, 'ftp://' + row.fastq_ftp.split(';').takeRight(2)[0], 'ftp://' + row.fastq_ftp.split(';').takeRight(2)[1]) }
//            .map { row -> tuple(row.run_accession, row.sample_accession, 'ftp://' + row.fastq_ftp.split(';')[0], 'ftp://' + row.fastq_ftp.split(';')[1]) }

    map_to_reference(data, params.SARS2_FA, params.SARS2_FA_FAI, params.SECRETS, params.STUDY, params.GCLOUD_PROJECT_ID, params.GCLOUD_KEY_FILE)
    ena_analysis_submit(map_to_reference.out, params.SECRETS, params.STUDY, params.TEST_SUBMISSION, params.CONFIG_YAML)
}