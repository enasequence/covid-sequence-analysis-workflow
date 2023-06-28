#!/usr/bin/env bash
run_accession=${1}
sars2_fasta=${2}
task_cpus=${3}

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
cat ${run_accession}.pileup | awk '{print $2,",",$3,",",$4}' > ${run_accession}.coverage
