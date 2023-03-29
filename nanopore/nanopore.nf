#!/usr/bin/env nextflow

//params.SARS2_FA = "gs://prj-int-dev-covid19-nf-gls/data/NC_045512.2.fa"
//params.SARS2_FA_FAI = "gs://prj-int-dev-covid19-nf-gls/data/NC_045512.2.fa.fai"
//params.SECRETS = "gs://prj-int-dev-covid19-nf-gls/prepro/projects_accounts.csv"

params.INDEX = "gs://prj-int-dev-covid19-nf-gls/prepro/nanopore.index.tsv"
params.STOREDIR = "gs://prj-int-dev-covid19-nf-gls/prepro/storeDir"
params.OUTDIR = "gs://prj-int-dev-covid19-nf-gls/prepro/results"

params.STUDY = 'PRJEB45555'
params.TEST_SUBMISSION = 'true'

//import nextflow.splitter.CsvSplitter
nextflow.enable.dsl = 2

//def fetchRunAccessions(String tsv ) {
//    CsvSplitter splitter = new CsvSplitter().options( header:true, sep:'\t' )
//    BufferedReader reader = new BufferedReader( new FileReader( tsv ) )
//    splitter.parseHeader( reader )
//
//    List<String> run_accessions = [] as List<String>
//    Map<String,String> row
//    while( row = splitter.fetchRecord( reader ) ) {
//        run_accessions.add( row['run_accession'] )
//    }
//    return run_accessions
//}

process map_to_reference {
    storeDir params.STOREDIR

    cpus 4 /* more is better, parallelizes very well*/
    memory '8 GB'
    container 'davidyuyuan/ena-sars-cov2-nanopore:2.0'

    input:
    tuple val(run_accession), val(sample_accession), file(input_file)
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
    line="\$(grep ${study_accession} ${projects_accounts_csv})"
    ftp_id="\$(echo \${line} | cut -d ',' -f 3)"
    ftp_password="\$(echo \${line} | cut -d ',' -f 6)"
    
    if [ "\${ftp_id}" = 'public' ]; then
        wget -t 0 -O ${run_accession}_1.fastq.gz \$(cat ${input_file})
    else
        wget -t 0 -O ${run_accession}_1.fastq.gz \$(cat ${input_file}) --user=\${ftp_id} --password=\${ftp_password}
    fi
    cutadapt -u 30 -u -30 -o ${run_accession}.trimmed.fastq ${run_accession}_1.fastq.gz -m 75 -j ${task.cpus} --quiet

    minimap2 -Y -t ${task.cpus} -x map-ont -a ${sars2_fasta} ${run_accession}.trimmed.fastq | samtools view -bF 4 - | samtools sort -@ ${task.cpus} - > ${run_accession}.bam
    samtools index -@ ${task.cpus} ${run_accession}.bam

    bam_to_vcf.py -b ${run_accession}.bam -r ${sars2_fasta} --mindepth 30 --minAF 0.1 -c ${task.cpus} -o ${run_accession}.vcf
    filtervcf.py -i ${run_accession}.vcf -o ${run_accession}_filtered.vcf
    bgzip ${run_accession}_filtered.vcf

    samtools mpileup -a -A -Q 0 -d 8000 -f ${sars2_fasta} ${run_accession}.bam > ${run_accession}.pileup
    cat ${run_accession}.pileup | awk '{print \$2,","\$3,","\$4}' > ${run_accession}.coverage
    tabix ${run_accession}_filtered.vcf.gz
    # vcf2consensus.py -v ${run_accession}_filtered.vcf.gz -d ${run_accession}.coverage -r ${sars2_fasta} -o ${run_accession}_consensus.fasta -dp 30 -n ${run_accession}
    vcf2consensus.py -v ${run_accession}_filtered.vcf.gz -d ${run_accession}.coverage -r ${sars2_fasta} -o headless_consensus.fasta -dp 30 -n ${run_accession}
    fix_consensus_header.py headless_consensus.fasta > ${run_accession}_consensus.fasta
    bgzip ${run_accession}.coverage

    java -Xmx4g -jar /opt/conda/share/snpeff-5.0-1/snpEff.jar -q -no-downstream -no-upstream -noStats NC_045512.2 ${run_accession}.vcf > ${run_accession}.annot.vcf
    #bgzip ${run_accession}.vcf
    bgzip ${run_accession}.annot.vcf
    bgzip ${run_accession}_consensus.fasta
    """
}

include { ena_analysis_submit } from '../nextflow-lib/ena-analysis-submitter.nf'
workflow {
//    Requires local input.
//    accessions = fetchRunAccessions(params.INDEX)
//    data = Channel.fromSRA( accessions )
//    data.view()
    data = Channel
            .fromPath(params.INDEX)
            .splitCsv(header: true, sep: '\t')
            .map { row -> tuple(row.run_accession, row.sample_accession, 'ftp://' + row.fastq_ftp) }

    map_to_reference(data, params.SARS2_FA, params.SARS2_FA_FAI, params.SECRETS, params.STUDY)
    ena_analysis_submit(map_to_reference.out, params.SECRETS, params.STUDY, params.TEST_SUBMISSION, params.CONFIG_YAML)
}
