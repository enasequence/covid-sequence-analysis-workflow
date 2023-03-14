#!/usr/bin/env bash

# DIR where the current script resides
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

pipeline=${1:-'nanopore'}
profile=${2:-'gls'}
root_dir=${3:-'gs://prj-int-dev-covid19-nf-gls'}
snapshot_date=${4:-'2023-03-14'}
concurrency=${5:-'5'}   # Maximum concurrency determined by the bottleneck - the submission server at present
batch_size=${6:-'5'}   # takes 12 hours if 2 jobs, 72 hours if 40 jobs
dataset_name=${7:-'sarscov2_metadata'}
project_id=${8:-'prj-int-dev-covid19-nf-gls'}
test_submission='true'

# Row count and batches
table_name="${pipeline}_to_be_processed"
sql="SELECT count(*) AS total FROM ${project_id}.${dataset_name}.${table_name}"
row_count=$(bq --project_id="${project_id}" --format=csv query --use_legacy_sql=false "${sql}" | grep -v total)

############################################
# as defined as queueSize in nextflow.config
############################################
queue_size=100
batches=$(( (row_count + batch_size - 1) / batch_size ))


batches=$(( (row_count + batch_size - 1) / batch_size ))
num_of_jobs=$(( (concurrency + queue_size - 1) / queue_size ))

input_dir="${DIR}/data/${snapshot_date}"; mkdir -p "${input_dir}"

echo "Pipeline: ${pipeline}"
echo "batches: ${batches}"
echo "num_of_jobs: ${num_of_jobs}"
echo "input_dir: ${input_dir}"
sleep 5
#for(( i=0; i<batches; i+=num_of_jobs )); do
#  for (( j=i; j<i+num_of_jobs&&j<batches; j++ )); do
#sql="SELECT * FROM ${project_id}.${dataset_name}.${table_name} LIMIT ${batch_size} OFFSET ${offset}"

for (( j=0; j<num_of_jobs&&j<batches; j++ )); do
#  bsub -n 2 -M 8192 -q production
  echo "root_dir: ${root_dir}"
  echo "pwd 1: $PWD"
  mkdir -p "${root_dir}/${snapshot_date}/${pipeline}_${j}"; cd "${root_dir}/${snapshot_date}/${pipeline}_${j}" || exit
  echo "pwd 2: $PWD"
  offset=$((j * batch_size))
  "${DIR}/manage_nf.sh" "${j}" "${batch_size}" "${project_id}" "${dataset_name}" "${table_name}" \
  "${input_dir}" "${pipeline}" "${profile}" "${root_dir}" "${snapshot_date}" "${test_submission}"
done

#Move set.archived.sh and update.receipt.sh outside the loop
#Comment out update.receipt.sh
#"${DIR}/update.receipt.sh" "${batch_index}" "${snapshot_date}" "${pipeline}" "${profile}" "${root_dir}" "${dataset_name}" "${project_id}"
"${DIR}/set.archived.sh" "${dataset_name}" "${project_id}"

sql="CREATE OR REPLACE TABLE ${dataset_name}.sra_processing AS SELECT DISTINCT * FROM ${dataset_name}.sra_processing"
bq --project_id="${project_id}" --format=csv query --use_legacy_sql=false "${sql}"

#max_mem avg_mem swap stat exit_code exec_cwd exec_host
#bjobs -u all -d -o "jobid job_name user submit_time start_time finish_time run_time cpu_used slots min_req_proc max_req_proc nthreads delimiter='^'" > jobs.csv
num_of_snapshots=$(( batches / num_of_jobs + 1 ))
echo "Row count: ${row_count}. Total number of batches: ${batches}, Number of jobs: ${num_of_jobs}, Number of snapshots: ${num_of_snapshots}."
wait
