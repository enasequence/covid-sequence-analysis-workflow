#!/usr/bin/env bash

# DIR where the current script resides
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

pipeline=${1:-'nanopore'}
profile=${2:-'gls'}
batch_index=${3:-'0'}
snapshot_date=${4:-'2022-03-22'}
batch_size=${5:-'1'}
test_submission=${6:-'true'}
study_accession=${7:-'PRJEB45555'}
root_dir=${8:-'gs://prj-int-dev-covid19-nf-gls'}
dataset_name=${9:-'sarscov2_metadata'}
project_id=${10:-'prj-int-dev-covid19-nf-gls'}

##############################
# Retrieve and reserve a batch
##############################
offset=$((batch_index * batch_size))
echo "** Retrieving and reserving batch ${batch_index} with the size of ${batch_size} from the offset of ${offset}. **"

output_dir="${DIR}/results/${snapshot_date}"
mkdir -p "${output_dir}"

table_name="${pipeline}_to_be_processed"
sql="SELECT * FROM ${project_id}.${dataset_name}.${table_name} LIMIT ${batch_size} OFFSET ${offset}"
bq --project_id="${project_id}" --format=csv query --use_legacy_sql=false --max_rows="${batch_size}" "${sql}" \
  | awk 'BEGIN{ FS=","; OFS="\t" }{$1=$1; print $0 }' > "${output_dir}/${table_name}_${batch_index}.tsv"
gsutil -m cp "${output_dir}/${table_name}_${batch_index}.tsv" "gs://${dataset_name}/${table_name}_${batch_index}.tsv" && \
  bq --project_id="${project_id}" load --source_format=CSV --replace=false --skip_leading_rows=1 --field_delimiter=tab \
  --max_bad_records=0 "${dataset_name}.sra_processing" "gs://${dataset_name}/${table_name}_${batch_index}.tsv"

#################################
# Process the batch with Nextflow
#################################
echo "** Processing samples with ${DIR}/${pipeline}/main.nf. **"

pipeline_dir="${root_dir}/${snapshot_date}/${pipeline}_${batch_index}"
nextflow -C "${DIR}/nextflow-lib/nextflow.config" run "${DIR}/${pipeline}/main.nf" -profile "${profile}" \
      --STUDY "${study_accession}" \
      --SECRETS "${DIR}/data/projects_accounts.csv" \
      --SARS2_FA "${DIR}/data/NC_045512.2.fa" \
      --SARS2_FA_FAI "${DIR}/data/NC_045512.2.fa.fai"\
      --TEST_SUBMISSION "${test_submission}" \
      --INDEX "${output_dir}/${table_name}_${batch_index}.tsv" \
      --OUTDIR "${pipeline_dir}/publishDir" \
      --STOREDIR "${pipeline_dir}/storeDir" \
      -w "${pipeline_dir}/workDir" \
      -with-tower

########################################################################################
# Update submission receipt and submission metadata as well as all the analyses archived
########################################################################################
"${DIR}/update.receipt.sh" "${batch_index}" "${snapshot_date}" "${pipeline}" "${profile}" "${root_dir}" "${dataset_name}" "${project_id}"
"${DIR}/set.archived.sh" "${dataset_name}" "${project_id}"
