#!/usr/bin/env bash

batch_index=${1:-'0'}
snapshot_date=${2:-'2022-06-27'}
pipeline=${3:-'illumina'}       # nanopore
profile=${4:-'codon'}   # gls
root_dir=${5:-'/hps/nobackup/tburdett/ena/users/analyser/nextflow'}   # gs://prj-int-dev-covid19-nf-gls
dataset_name=${6:-'sarscov2_metadata'}
project_id=${7:-'prj-int-dev-covid19-nf-gls'}

# DIR where the current script resides
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
#export TMPDIR=/lscratch

pipeline_dir="${root_dir}/${snapshot_date}/${pipeline}_${batch_index}"
output_dir="${DIR}/results/${snapshot_date}/output"; rm -R "${output_dir}"; mkdir -p "${output_dir}"

###########################################################################################
# Concatenate the receipts in a batch
# Analysis archived: ${output_dir}/${snapshot_date}_${pipeline}_${batch_index}_receipts.tsv
###########################################################################################
echo ""
echo "** Getting ${snapshot_date}_${pipeline}_${batch_index}_receipts. **"

printf '%s\t%s\t%s\t%s\n' "analysis_accession" "file_submitted" "time_submitted" "snapshot_date" > \
  "${output_dir}/${snapshot_date}_${pipeline}_${batch_index}_receipts.tsv"
if [ "${profile}" = 'gls' ]; then
#  gsutil -m ls "${pipeline_dir}/publishDir/**/*_submissions.txt" > "${output_dir}/${pipeline}_receipts_${batch_index}.txt"
  gsutil -m cp "${pipeline_dir}/publishDir/**/*_submissions.txt" "${output_dir}"
#  gsutil -m cp "gs://prj-int-dev-covid19-nf-gls/prepro/results/**/*_submissions.txt" "${output_dir}"
  find "${output_dir}" -type f -name '*_submissions.txt' -exec cat {} + >> "${output_dir}/${snapshot_date}_${pipeline}_${batch_index}_receipts.tsv"
else
#  find "${pipeline_dir}/publishDir/" -type f -name '*_submissions.txt' > "${output_dir}/${pipeline}_receipts_${batch_index}.txt"
#  cat "$(cat "${output_dir}/${pipeline}_receipts_${batch_index}.txt")" >> "${output_dir}/${snapshot_date}_${pipeline}_${batch_index}_receipts.tsv"
  find "${pipeline_dir}/publishDir/" -type f -name '*_submissions.txt' -exec cat {} + >> \
    "${output_dir}/${snapshot_date}_${pipeline}_${batch_index}_receipts.tsv"
fi

########################################################
# upload receipts and update the rows with snapshot_date
########################################################
echo ""
echo "** Updating ${dataset_name}.submission_receipts table. **"

gsutil -m cp "${output_dir}/${snapshot_date}_${pipeline}_${batch_index}_receipts.tsv" \
  "gs://${dataset_name}/${snapshot_date}_${pipeline}_${batch_index}_receipts.tsv" && \
  bq --project_id="${project_id}" load --source_format=CSV --replace=false --skip_leading_rows=1 --field_delimiter=tab \
  --autodetect "${dataset_name}.submission_receipts" "gs://${dataset_name}/${snapshot_date}_${pipeline}_${batch_index}_receipts.tsv" \
  "analysis_accession:STRING,file_submitted:STRING,time_submitted,snapshot_date"

sql="UPDATE ${dataset_name}.submission_receipts SET snapshot_date = '""${snapshot_date}""' WHERE snapshot_date is NULL"
bq --project_id="${project_id}" --format=csv query --use_legacy_sql=false "${sql}"
sql="CREATE OR REPLACE TABLE ${dataset_name}.submission_receipts AS SELECT DISTINCT * FROM ${dataset_name}.submission_receipts"
bq --project_id="${project_id}" --format=csv query --use_legacy_sql=false "${sql}"
rm "${output_dir}/${snapshot_date}_${pipeline}_${batch_index}_receipts.tsv"
