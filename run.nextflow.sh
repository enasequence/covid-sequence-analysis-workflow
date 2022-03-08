#!/usr/bin/env bash

# DIR where the current script resides
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

batch_index=${1:-'0'}
snapshot_date=${2:-'2022-02-25'}
batch_size=${3:-'1'}
test_submission=${4:-'true'}
pipeline=${5:-'nanopore'}
profile=${6:-'gls'}
root_dir=${7:-'gs://prj-int-dev-covid19-nf-gls'}
dataset_name=${8:-'sarscov2_metadata'}
project_id=${9:-'prj-int-dev-covid19-nf-gls'}

offset=$((batch_index * batch_size))
table_name="${pipeline}_to_be_processed"
pipeline_dir="${root_dir}/${snapshot_date}/${pipeline}_${batch_index}"
output_dir="${DIR}/results/${snapshot_date}"
mkdir -p "${output_dir}"

##############################
# Retrieve and reserve a batch
##############################
echo "** Retrieving and reserving batch ${batch_index} with the size of ${batch_size} from the offset of ${offset}. **"

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

nextflow -C "${DIR}/nextflow-lib/nextflow.config" run "${DIR}/${pipeline}/main.nf" -profile "${profile}" \
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
