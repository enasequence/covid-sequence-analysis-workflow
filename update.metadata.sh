#!/usr/bin/env bash

j=${1:-'0'}
snapshot_date=${2:-'2022-02-25'}
pipeline=${3:-'nanopore'}
dataset_name=${4:-'sarscov2_metadata'}
project_id=${5:-'prj-int-dev-covid19-nf-gls'}

# DIR where the current script resides
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

output_dir="${DIR}/results/${snapshot_date}/output"
mkdir -p "${output_dir}"

table_name="${pipeline}_to_be_processed"

# Results and metadata
function gen_metadata {
  local input_file=$1
  local index_tsv=$2
  local metadata=$3
  local timestamp=$4

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "run_id"	"platform"	"model"	"first_public"	"first_created"	"country"	"collection_date"	"snapshot_date" > "${metadata}"
  true > "${metadata}.tmp"
  sed 1d "${input_file}" | while IFS="" read -r receipt || [ -n "$f" ]
  do
    run_accession=$(echo "${receipt}" | cut -f 2 | cut -d '_' -f 1 )
    sample_accession=$(grep "${run_accession}" "${index_tsv}" | cut -f3-4,10-13)
    printf '%s\t%s\t%s\n' "${run_accession}" "${sample_accession}" "${timestamp}" >> "${metadata}.tmp"
  done
  sort -u "${metadata}.tmp" >> "${metadata}" && rm "${metadata}.tmp"
}

##################################
# Update submission_metadata table
##################################
gen_metadata "${output_dir}/${snapshot_date}_${pipeline}_${j}_receipts.tsv" "${DIR}/results/${snapshot_date}/${table_name}_${j}.tsv" "${output_dir}/${pipeline}_metadata_${j}.tsv" "${snapshot_date}"
gsutil -m cp "${output_dir}/${pipeline}_metadata_${j}.tsv" "gs://${dataset_name}/${pipeline}_metadata_${j}.tsv"
bq --project_id="${project_id}" load --source_format=CSV --replace=false --skip_leading_rows=1 --field_delimiter=tab \
  --autodetect "${dataset_name}.submission_metadata" "gs://${dataset_name}/${pipeline}_metadata_${j}.tsv" \
  "run_id,platform,model,first_public,first_created,country,collection_date,snapshot_date"
