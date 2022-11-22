#!/usr/bin/env bash

# DIR where the current script resides
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
export $(grep -v '^#' .env | xargs)

pipeline=${1:-'illumina'}
profile=${2:-'gls'}
root_dir=${3:-'gs://prj-int-dev-covid19-nf-gls'}
snapshot_date=${4:-'2022-06-26'} #2022-06-27
concurrency=${5:-'120'}   # Maximum concurrency determined by the bottleneck - the submission server at present
batch_size=${6:-'9000'}   # takes 12 hours if 2 jobs, 72 hours if 40 jobs
dataset_name=${7:-'sarscov2_metadata'}
project_id=${8:-'prj-int-dev-covid19-nf-gls'}

# Row count and batches
table_name="${pipeline}_to_be_processed"
sql="SELECT count(*) AS total FROM ${project_id}.${dataset_name}.${table_name}"
row_count=$(bq --project_id="${project_id}" --format=csv query --use_legacy_sql=false "${sql}" | grep -v total)

batch_size=5
row_count=4
############################################
# as defined as queueSize in nextflow.config
############################################
queue_size=4
batches=$(( row_count / batch_size + 1 ))
num_of_jobs=$(( concurrency / queue_size ))
poc_project_id='prj-int-dev-ena-nft-poc-368211'
root_dir='gs://prj-int-dev-covid19-nf-gls-poc'
input_dir="${DIR}/data/${snapshot_date}"; mkdir -p "${input_dir}"
offset=0
test_submission='true'
#for(( i=0; i<batches; i+=num_of_jobs )); do
#  for (( j=i; j<i+num_of_jobs&&j<batches; j++ )); do
for (( j=0; j<num_of_jobs&&j<batches; j++ )); do
  # sql="SELECT * FROM ${project_id}.${dataset_name}.${table_name} LIMIT ${batch_size} OFFSET ${offset}"
  # bq --project_id="${project_id}" --format=csv query --use_legacy_sql=false --max_rows="${batch_size}" "${sql}" \
  #   | awk 'BEGIN{ FS=","; OFS="\t" }{$1=$1; print $0 }' > "${input_dir}/${table_name}_${j}.tsv"
  # gsutil -m cp "${input_dir}/${table_name}_${j}.tsv" "gs://${dataset_name}_poc/${table_name}_${j}.tsv" && \
  #   bq --project_id="${project_id}" load --source_format=CSV --replace=false --skip_leading_rows=1 --field_delimiter=tab \
  #   --max_bad_records=0 "${dataset_name}.sra_processing" "gs://${dataset_name}_poc/${table_name}_${j}.tsv"

#  bsub -n 2 -M 8192 -q production
  # "${DIR}/run.nextflow.sh" "${pipeline}" "${profile}" "${root_dir}" "${j}" "${snapshot_date}" "${batch_size}" &
  
  "${DIR}/test.nextflow.sh" "${input_dir}/${table_name}_${j}.tsv" \
    "${pipeline}" "${profile}" "${root_dir}" "${j}" "${snapshot_date}" "${test_submission}"
done
#done

#max_mem avg_mem swap stat exit_code exec_cwd exec_host
#bjobs -u all -d -o "jobid job_name user submit_time start_time finish_time run_time cpu_used slots min_req_proc max_req_proc nthreads delimiter='^'" > jobs.csv
num_of_snapshots=$(( batches / num_of_jobs + 1 ))
echo "Row count: ${row_count}. Total number of batches: ${batches}, Number of jobs: ${num_of_jobs}, Number of snapshots: ${num_of_snapshots}."
wait
