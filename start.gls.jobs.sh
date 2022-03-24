#!/usr/bin/env bash

# DIR where the current script resides
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

pipeline=${1:-'illumina'}
profile=${2:-'gls'}
root_dir=${3:-'gs://prj-int-dev-covid19-nf-gls'}
snapshot_date=${4:-'2022-04-12'}
concurrency=${5:-'120'}   # 120 Maximum concurrency determined by the bottleneck - the submission server at present
batch_size=${6:-'9000'}   # 9000 takes 12 hours if 2 jobs, 72 hours if 40 jobs
dataset_name=${7:-'sarscov2_metadata'}
project_id=${8:-'prj-int-dev-covid19-nf-gls'}

# Row count and batches
table_name="${pipeline}_to_be_processed"
sql="SELECT count(*) AS total FROM ${project_id}.${dataset_name}.${table_name}"
row_count=$(bq --project_id="${project_id}" --format=csv query --use_legacy_sql=false "${sql}" | grep -v total)

############################################
# as defined as queueSize in nextflow.config
############################################
queue_size=5
batches=$(( row_count / batch_size + 1 ))
num_of_jobs=$(( concurrency / queue_size ))

#for(( i=0; i<batches; i+=num_of_jobs )); do
#  for (( j=i; j<i+num_of_jobs&&j<batches; j++ )); do
for (( j=0; j<num_of_jobs&&j<batches; j++ )); do
#  bsub -n 2 -M 4096 -q production
  "${DIR}/run.nextflow.sh" "${pipeline}" "${profile}" "${root_dir}" "${j}" "${snapshot_date}" "${batch_size}" &
done
#sleep 12h
#done

#max_mem avg_mem swap stat exit_code exec_cwd exec_host
#bjobs -u all -d -o "jobid job_name user submit_time start_time finish_time run_time cpu_used slots min_req_proc max_req_proc nthreads delimiter='^'" > jobs.csv
num_of_snapshots=$(( batches / num_of_jobs + 1 ))
echo "Row count: ${row_count}. Total number of batches: ${batches}, Number of jobs: ${num_of_jobs}, Number of snapshots: ${num_of_snapshots}."
