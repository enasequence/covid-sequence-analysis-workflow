#!/usr/bin/env bash
set -e
run_accession=${1}
sars2_fasta=${2}
task_cpus=${3}

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