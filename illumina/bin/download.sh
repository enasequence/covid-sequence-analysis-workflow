#!/usr/bin/env bash
run_accession=${1}
projects_accounts_csv=${2}
input_file_1=${3}
input_file_2=${4}
study_accession=${5}

line="$(grep ${study_accession} ${projects_accounts_csv})"
ftp_id="$(echo ${line} | cut -d ',' -f 3)"
ftp_password="$(echo ${line} | cut -d ',' -f 6)"

if [ "${ftp_id}" = 'public' ]; then
    wget -t 0 -O ${run_accession}_1.fastq.gz $(cat ${input_file_1})
    wget -t 0 -O ${run_accession}_2.fastq.gz $(cat ${input_file_2})
else
    wget -t 0 -O ${run_accession}_1.fastq.gz $(cat ${input_file_1}) --user=${ftp_id} --password=${ftp_password}
    wget -t 0 -O ${run_accession}_2.fastq.gz $(cat ${input_file_2}) --user=${ftp_id} --password=${ftp_password}
fi