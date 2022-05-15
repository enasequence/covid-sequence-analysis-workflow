#!/usr/bin/env bash

# DIR where the current script resides
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

skip=${1:-'0'}
concurrency=${2:-'400'}   # Maximum concurrency determined by the bottleneck - the submission server and storage space
pipeline=${3:-'illumina'}   # nanopore
root_dir=${4:-'/hps/nobackup/tburdett/ena/users/analyser/nextflow'}  # /hps/nobackup/cochrane/ena/users/analyser/nextflow
batch_size=${5:-'15000'}
profile=${6:-'codon'}
snapshot_date=${7:-'2022-05-23'}  # 2022-03-22 2022-04-12 2022-05-23 2022-06-27
dataset_name=${8:-'sarscov2_metadata'}
project_id=${9:-'prj-int-dev-covid19-nf-gls'}

# Row count and batches
table_name="${pipeline}_to_be_processed"
sql="SELECT count(*) AS total FROM ${project_id}.${dataset_name}.${table_name}"
row_count=$(bq --project_id="${project_id}" --format=csv query --use_legacy_sql=false "${sql}" | grep -v total)

############################################
# as defined as queueSize in nextflow.config
############################################
queue_size=100     # 4
batches=$(( row_count / batch_size + 1 ))
num_of_jobs=$(( concurrency / queue_size ))
#mem_limit=$(( batch_size / 2500 * 2048));mem_limit=$(( mem_limit > 2048 ? mem_limit : 2048 ))

input_dir="${DIR}/data/${snapshot_date}"; mkdir -p "${input_dir}"

for (( batch_index=skip; batch_index<skip+num_of_jobs&&batch_index<batches; batch_index++ )); do
  mkdir -p "${root_dir}/${pipeline}_${batch_index}"; cd "${root_dir}/${pipeline}_${batch_index}" || exit

  offset=$((batch_index * batch_size))
  echo ""
  echo "** Retrieving and reserving batch ${batch_index} with the size of ${batch_size} from the offset of ${offset}. **"

  sql="SELECT * FROM ${project_id}.${dataset_name}.${table_name} LIMIT ${batch_size} OFFSET ${offset}"
  bq --project_id="${project_id}" --format=csv query --use_legacy_sql=false --max_rows="${batch_size}" "${sql}" \
    | awk 'BEGIN{ FS=","; OFS="\t" }{$1=$1; print $0 }' > "${input_dir}/${table_name}_${batch_index}.tsv"
  gsutil -m cp "${input_dir}/${table_name}_${batch_index}.tsv" "gs://${dataset_name}/${table_name}_${batch_index}.tsv" && \
    bq --project_id="${project_id}" load --source_format=CSV --replace=false --skip_leading_rows=1 --field_delimiter=tab \
    --max_bad_records=0 "${dataset_name}.sra_processing" "gs://${dataset_name}/${table_name}_${batch_index}.tsv"

  bsub -n 2 -M 4096 -q production "${DIR}/run.nextflow.sh" "${input_dir}/${table_name}_${batch_index}.tsv" \
    "${pipeline}" "${profile}" "${root_dir}" "${batch_index}" "${snapshot_date}"
done

sql="CREATE OR REPLACE TABLE ${dataset_name}.sra_processing AS SELECT DISTINCT * FROM ${dataset_name}.sra_processing"
bq --project_id="${project_id}" --format=csv query --use_legacy_sql=false "${sql}"

#max_mem avg_mem swap stat exit_code exec_cwd exec_host
#bjobs -u all -d -o "jobid job_name user submit_time start_time finish_time run_time cpu_used slots min_req_proc max_req_proc nthreads delimiter='^'" > jobs.csv
num_of_snapshots=$(( batches / num_of_jobs + 1 ))
echo "Row count: ${row_count}. Total number of batches: ${batches}, Number of jobs: ${num_of_jobs}, Number of snapshots: ${num_of_snapshots}."
