// Same process duplicated in both pipelines due to a limitation in Nextflow
process ena_analysis_submit {
    publishDir params.OUTDIR, mode: 'copy'
    storeDir params.STOREDIR

    container 'davidyuyuan/ena-analysis-submitter:2.0'

    input:
    val(run_accession)
    val(sample_accession)
    file(output_tgz)
    file(filtered_vcf_gz)
    file(consensus_fasta_gz)
    path(projects_accounts_csv)
    val(study_accession)

    output:
    file("${run_accession}_output/${study_accession}/${run_accession}_output.tar.gz")
    file("${run_accession}_output/${study_accession}/${run_accession}_filtered.vcf.gz")
    file("${run_accession}_output/${study_accession}/${run_accession}_consensus.fasta.gz")
    file("${run_accession}_output/${study_accession}/${run_accession}_submissions.txt")

    script:
    """
    line="\$(grep ${study_accession} ${projects_accounts_csv})"
    webin_id="\$(echo \${line} | cut -d ',' -f 4)"
    webin_password="\$(echo \${line} | cut -d ',' -f 5)"
    
    mkdir -p ${run_accession}_output/${study_accession}
    cp nextflow-bin/config.yaml ${run_accession}_output/${study_accession}
    if [ "${study_accession}" = 'PRJEB45555' ]; then
        analysis_submission.py -t -o ${run_accession}_output/${study_accession} -p PRJEB43947 -s ${sample_accession} -r ${run_accession} -f ${output_tgz} -a PATHOGEN_ANALYSIS -au \${webin_id} -ap \${webin_password}
        analysis_submission.py -t -o ${run_accession}_output/${study_accession} -p PRJEB45554 -s ${sample_accession} -r ${run_accession} -f ${filtered_vcf_gz} -a COVID19_FILTERED_VCF -au \${webin_id} -ap \${webin_password}
        analysis_submission.py -t -o ${run_accession}_output/${study_accession} -p PRJEB45619 -s ${sample_accession} -r ${run_accession} -f ${consensus_fasta_gz} -a COVID19_CONSENSUS -au \${webin_id} -ap \${webin_password}
    else
        analysis_submission.py -t -o ${run_accession}_output/${study_accession} -p ${study_accession} -s ${sample_accession} -r ${run_accession} -f ${output_tgz} -a PATHOGEN_ANALYSIS -au \${webin_id} -ap \${webin_password}
        analysis_submission.py -t -o ${run_accession}_output/${study_accession} -p ${study_accession} -s ${sample_accession} -r ${run_accession} -f ${filtered_vcf_gz} -a COVID19_FILTERED_VCF -au \${webin_id} -ap \${webin_password}
        analysis_submission.py -t -o ${run_accession}_output/${study_accession} -p ${study_accession} -s ${sample_accession} -r ${run_accession} -f ${consensus_fasta_gz} -a COVID19_CONSENSUS -au \${webin_id} -ap \${webin_password}
    fi
    mv ${output_tgz} ${filtered_vcf_gz} ${consensus_fasta_gz} ${run_accession}_output/${study_accession}
    mv ${run_accession}_output/${study_accession}/successful_submissions.txt ${run_accession}_output/${study_accession}/${run_accession}_submissions.txt
    """
}

