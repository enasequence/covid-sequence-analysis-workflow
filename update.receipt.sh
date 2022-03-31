#!/usr/bin/env bash

batch_index=${1:-'0'}
snapshot_date=${2:-'2022-02-25'}
pipeline=${3:-'nanopore'}
profile=${4:-'gls'}
root_dir=${5:-'gs://prj-int-dev-covid19-nf-gls'}
dataset_name=${6:-'sarscov2_metadata'}
project_id=${7:-'prj-int-dev-covid19-nf-gls'}

# DIR where the current script resides
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

output_dir="${DIR}/results/${snapshot_date}/output"
rm -R "${output_dir}"
mkdir -p "${output_dir}"

###########################################################################################
# Concatenate the receipts in a batch
# Analysis archived: ${output_dir}/${snapshot_date}_${pipeline}_${batch_index}_receipts.tsv
###########################################################################################
echo "** Getting ${snapshot_date}_${pipeline}_${batch_index}_receipts. **"

printf '%s\t%s\t%s\t%s\n' "analysis_accession" "file_submitted" "time_submitted" "snapshot_date" > \
  "${output_dir}/${snapshot_date}_${pipeline}_${batch_index}_receipts.tsv"
if [ "${profile}" = 'gls' ]; then
#  gsutil -m ls "${root_dir}/${snapshot_date}/${pipeline}_${batch_index}/publishDir/**/*_submissions.txt" > "${output_dir}/${pipeline}_receipts_${batch_index}.txt"
  gsutil -m cp "${root_dir}/${snapshot_date}/${pipeline}_${batch_index}/publishDir/**/*_submissions.txt" "${output_dir}"
#  gsutil -m cp "gs://prj-int-dev-covid19-nf-gls/prepro/results/**/*_submissions.txt" "${output_dir}"
  find "${output_dir}" -type f -name '*_submissions.txt' -exec cat {} + >> "${output_dir}/${snapshot_date}_${pipeline}_${batch_index}_receipts.tsv"
else
#  find "${root_dir}/${snapshot_date}/${pipeline}_${batch_index}/publishDir/" -type f -name '*_submissions.txt' > "${output_dir}/${pipeline}_receipts_${batch_index}.txt"
#  cat "$(cat "${output_dir}/${pipeline}_receipts_${batch_index}.txt")" >> "${output_dir}/${snapshot_date}_${pipeline}_${batch_index}_receipts.tsv"
  find "${root_dir}/${snapshot_date}/${pipeline}_${batch_index}/publishDir/" -type f -name '*_submissions.txt' -exec cat {} + >> \
    "${output_dir}/${snapshot_date}_${pipeline}_${batch_index}_receipts.tsv"
fi

########################################################
# upload receipts and update the rows with snapshot_date
########################################################
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

##################################
# Update submission_metadata table
##################################
#echo "** Updating ${dataset_name}.submission_metadata table. **"
#
#function gen_metadata {
#  local input_file=$1
#  local index_tsv=$2
#  local metadata=$3
#  local timestamp=$4
#
#  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "run_id"	"platform"	"model"	"first_public"	"first_created"	"country"	"collection_date"	"snapshot_date" > "${metadata}"
#  sed 1d "${input_file}" | while IFS="" read -r receipt || [ -n "$f" ]
#  do
#    run_accession=$(echo "${receipt}" | cut -f 2 | cut -d '_' -f 1 )
#    sample_accession=$(grep "${run_accession}" "${index_tsv}" | cut -f3-4,10-13)
#    printf '%s\t%s\t%s\n' "${run_accession}" "${sample_accession}" "${timestamp}" >> "${metadata}"
#  done
#}
#
#gen_metadata "${output_dir}/${snapshot_date}_${pipeline}_${batch_index}_receipts.tsv" \
#  "${DIR}/results/${snapshot_date}/${pipeline}_to_be_processed_${batch_index}.tsv" \
#  "${output_dir}/${pipeline}_metadata_${batch_index}.tsv" "${snapshot_date}"
#gsutil -m cp "${output_dir}/${pipeline}_metadata_${batch_index}.tsv" "gs://${dataset_name}/${pipeline}_metadata_${batch_index}.tsv" && \
#  bq --project_id="${project_id}" load --source_format=CSV --replace=false --skip_leading_rows=1 --field_delimiter=tab \
#  --autodetect "${dataset_name}.submission_metadata" "gs://${dataset_name}/${pipeline}_metadata_${batch_index}.tsv" \
#  "run_id,platform,model,first_public,first_created,country,collection_date,snapshot_date"
#
#sql="CREATE OR REPLACE TABLE ${dataset_name}.submission_metadata AS SELECT DISTINCT * FROM ${dataset_name}.submission_metadata WHERE platform IS NOT NULL"
#bq --project_id="${project_id}" --format=csv query --use_legacy_sql=false "${sql}"
#
#rm "${output_dir}/${snapshot_date}_${pipeline}_${batch_index}_receipts.tsv"
