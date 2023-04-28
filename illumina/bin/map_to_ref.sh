#!/usr/bin/env bash
set -e
set -o pipefail
run_accession=${1}
projects_accounts_csv=${2}
input_file_1=${3}
input_file_2=${4}
sars2_fasta=${5}
task_cpus=${6}
study_accession=${7}

line=$(grep ${study_accession} ${projects_accounts_csv})
ftp_id=$(echo ${line} | cut -d ',' -f 3)
ftp_password=$(echo ${line} | cut -d ',' -f 6)
if [ "${ftp_id}" = 'public' ]; then
    wget -t 0 -O ${run_accession}_1.fastq.gz ${input_file_1}
    wget -t 0 -O ${run_accession}_2.fastq.gz ${input_file_2}
else
    wget -t 0 -O ${run_accession}_1.fastq.gz ${input_file_1} --user=${ftp_id} --password=${ftp_password}
    wget -t 0 -O ${run_accession}_2.fastq.gz ${input_file_2} --user=${ftp_id} --password=${ftp_password}
fi

trimmomatic PE ${run_accession}_1.fastq.gz ${run_accession}_2.fastq.gz ${run_accession}_trim_1.fq \
${run_accession}_trim_1_un.fq ${run_accession}_trim_2.fq ${run_accession}_trim_2_un.fq \
-summary ${run_accession}_trim_summary -threads ${task_cpus} \
SLIDINGWINDOW:5:30 MINLEN:50

bwa index ${sars2_fasta}
bwa mem -t ${task_cpus} ${sars2_fasta} ${run_accession}_trim_1.fq ${run_accession}_trim_2.fq | samtools view -bF 4 - | samtools sort - > ${run_accession}_paired.bam
bwa mem -t ${task_cpus} ${sars2_fasta} <(cat ${run_accession}_trim_1_un.fq ${run_accession}_trim_2_un.fq) | samtools view -bF 4 - | samtools sort - > ${run_accession}_unpaired.bam
samtools merge ${run_accession}.bam ${run_accession}_paired.bam ${run_accession}_unpaired.bam
rm ${run_accession}_paired.bam ${run_accession}_unpaired.bam

samtools mpileup -a -A -Q 30 -d 8000 -f ${sars2_fasta} ${run_accession}.bam > ${run_accession}.pileup
cat ${run_accession}.pileup | awk '{print \$2,","\$3,","\$4}' > ${run_accession}.coverage

samtools index ${run_accession}.bam
lofreq indelqual --dindel ${run_accession}.bam -f ${sars2_fasta} -o ${run_accession}_fixed.bam
samtools index ${run_accession}_fixed.bam
lofreq call-parallel --no-default-filter --call-indels --pp-threads ${task_cpus} -f ${sars2_fasta} -o ${run_accession}.vcf ${run_accession}_fixed.bam
lofreq filter --af-min 0.25 -i ${run_accession}.vcf -o ${run_accession}_filtered.vcf
bgzip ${run_accession}.vcf
bgzip ${run_accession}_filtered.vcf
tabix ${run_accession}.vcf.gz
bcftools stats ${run_accession}.vcf.gz > ${run_accession}.stat

snpEff -q -no-downstream -no-upstream -noStats NC_045512.2 ${run_accession}.vcf > ${run_accession}.annot.vcf
# vcf_to_consensus.py -dp 10 -af 0.25 -v ${run_accession}.vcf.gz -d ${run_accession}.coverage -o ${run_accession}_consensus.fasta -n ${run_accession} -r ${sars2_fasta}
vcf_to_consensus.py -dp 10 -af 0.25 -v ${run_accession}.vcf.gz -d ${run_accession}.coverage -o headless_consensus.fasta -n ${run_accession} -r ${sars2_fasta}
fix_consensus_header.py headless_consensus.fasta > ${run_accession}_consensus.fasta
bgzip ${run_accession}_consensus.fasta
bgzip ${run_accession}.coverage
bgzip ${run_accession}.annot.vcf