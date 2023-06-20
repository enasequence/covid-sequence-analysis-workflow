#!/usr/bin/env bash
index=${1:-'1'}
batch_size=${2:-'15000'}
project_id=${3:-'prj-int-dev-covid19-nf-gls'}
dataset_name=${4:-'sarscov2_metadata'}
table_name=${5:-'nanopore_to_be_processed'}
input_dir=${6:-'/data/'}
pipeline=${7:-'nanopore'}
profile=${8:-'gls'}
root_dir=${9:-'gls'}
snapshot_date=${10:-'2023/03/14'}
test_submission=${11:-'true'}

offset=$((index * batch_size))
echo ""
echo "** Retrieving and reserving batch ${index} with the size of ${batch_size} from the offset of ${offset}. **"

sql="SELECT * FROM ${project_id}.${dataset_name}.${table_name} LIMIT ${batch_size} OFFSET ${offset}"
bq --project_id="${project_id}" --format=csv query --use_legacy_sql=false --max_rows="${batch_size}" "${sql}" \
	| awk 'BEGIN{ FS=","; OFS="\t" }{$1=$1; print $0 }' > "${input_dir}/${table_name}_${index}.tsv"
	#Modify the file here with a timestamp field and then load into BQ
	#sql2 = "ALTER TABLE ${dataset_name}.sra_processing ADD COLUMN timestamp DATETIME DEFAULT CURRENT_DATETIME();"
	#bq --project_id="${project_id}" --format=csv query --use_legacy_sql=false "${sql2}"
gsutil -m cp "${input_dir}/${table_name}_${index}.tsv" "gs://${dataset_name}/${table_name}_${index}.tsv" && \
    bq --project_id="${project_id}" load --source_format=CSV --replace=false --skip_leading_rows=1 --field_delimiter=tab \
    --max_bad_records=0 "${dataset_name}.sra_processing" "gs://${dataset_name}/${table_name}_${index}.tsv"
if [ "$profile" == "codon" ]; then
	bsub -n 2 -M 4096 -q production "${DIR}/run.nextflow.sh" "${input_dir}/${table_name}_${index}.tsv" \
    "${pipeline}" "${profile}" "${root_dir}" "${index}" "${snapshot_date}" "${test_submission}"
else
	#TODO: Figure out WHY absolute path is needed
	"/opt/covid-sequence-analysis-workflow/run.nextflow.sh" "${input_dir}/${table_name}_${index}.tsv" "${pipeline}" "${profile}" "${root_dir}" \
	"${index}" "${snapshot_date}" "${test_submission}"  &
fi
