#!/usr/bin/env bash

# DIR where the current script resides
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
export $(grep -v '^#' .env | xargs)

skip=${1:-'0'}
# concurrency=${2:-'500'}   # 5,20 Maximum concurrency determined by the bottleneck - the submission server and storage space
concurrency=${2:-'5000'} # concurrency>=queue_size
pipeline=${3:-'illumina'}   # nanopore
# root_dir=${4:-'/hps/nobackup/tburdett/ena/users/analyser/nextflow'}  # /hps/nobackup/cochrane/ena/users/analyser/nextflow
root_dir=${4:-"/scratch/${PROJECT_ID}/nextflow"} 
# batch_size=${5:-'15000'}
batch_size=${5:-'3750'} # Storage capacity per job (Storage Max ~ 1500 samples in /project/) TODO: Submit array job if cannot handle large file download
profile=${6:-'slurm'}
snapshot_date=${7:-'2022-12-19'}  #2022-09-26 2022-10-24 2022-11-21 2022-12-19
dataset_name=${8:-'sarscov2_metadata'}
project_id=${9:-'prj-int-dev-covid19-nf-gls'}
test_submission=true
# Row count and batches
table_name="${pipeline}_to_be_processed"
# sql="SELECT count(*) AS total FROM ${project_id}.${dataset_name}.${table_name}"
# row_count=$(bq --project_id="${project_id}" --format=csv query --use_legacy_sql=false "${sql}" | grep -v total)
row_count=75000
############################################
# as defined as queueSize in nextflow.config
############################################
queue_size=40     #3,4:100 1-2:20 
batches=$(( row_count / batch_size + 1 ))
num_of_jobs=$(( concurrency / queue_size ))
echo "$batches $num_of_jobs"
#mem_limit=$(( batch_size / 2500 * 2048));mem_limit=$(( mem_limit > 2048 ? mem_limit : 2048 ))

input_dir="${DIR}/data/${snapshot_date}"; mkdir -p "${input_dir}"
parallel_job_limit=5

module use $HOME/.local/easybuild/modules/all
module load Nextflow/22.10.0
nextflow -version

echo "array_size: ${array_size}"
echo "$(( skip+num_of_jobs )) ${batches}"
for (( batch_index=skip;batch_index<skip+num_of_jobs&&batch_index<batches; batch_index++ )); do
  mkdir -p "${root_dir}/${pipeline}_${batch_index}"; cd "${root_dir}/${pipeline}_${batch_index}" || exit
  
  offset=$((batch_index * batch_size))
  echo ""
  echo "** Retrieving and reserving batch ${batch_index} with the size of ${batch_size} from the offset of ${offset}. **"
  sql="SELECT * FROM ${project_id}.${dataset_name}.${table_name} LIMIT ${batch_size} OFFSET ${offset}"
  # bq --project_id="${project_id}" --format=csv query --use_legacy_sql=false --max_rows="${batch_size}" "${sql}" \
  #   | awk 'BEGIN{ FS=","; OFS="\t" }{$1=$1; print $0 }' > "${input_dir}/${table_name}_${batch_index}.tsv"
  if [[ $test_submission = false ]]
  then
    echo "Start production"
    # gsutil -m cp "${input_dir}/${table_name}_${batch_index}.tsv" "gs://${dataset_name}/${table_name}_${batch_index}.tsv" && \
    # bq --project_id="${project_id}" load --source_format=CSV --replace=false --skip_leading_rows=1 --field_delimiter=tab \
    # --max_bad_records=0 "${dataset_name}.sra_processing" "gs://${dataset_name}/${table_name}_${batch_index}.tsv"
  else
    echo "Start test"
    # gsutil -m cp "${input_dir}/${table_name}_${batch_index}.tsv" "gs://${dataset_name}/${table_name}_${batch_index}.tsv" && \
    # bq --project_id="${project_id}" load --source_format=CSV --replace=false --skip_leading_rows=1 --field_delimiter=tab \
    # --max_bad_records=0 "${dataset_name}.sra_processing_test" "gs://${dataset_name}/${table_name}_${batch_index}.tsv"
  fi
  # sbatch --account=$PROJECT_ID -N 2 --ntasks=2 --cpus-per-task=32 -t 72:00:00 -p small --export ALL --wrap="${DIR}/run.nextflow.sh ${input_dir}/${table_name}_${batch_index}.tsv \
  #   ${pipeline} ${profile} ${root_dir} ${batch_index} ${snapshot_date} ${test_submission}"
done
index_range="0-$(($batches - 1 ))%$parallel_job_limit"
echo $index_range
sbatch --account=$PROJECT_ID -N 2 --mem 4096 --cpus-per-task=16 -t 72:00:00 -p small --export ALL --array=$index_range --wrap="${DIR}/array.nextflow.sh ${input_dir}/${table_name} \
    ${pipeline} ${profile} ${root_dir} ${batch_index} ${snapshot_date} ${test_submission}"

# sql="CREATE OR REPLACE TABLE ${dataset_name}.sra_processing AS SELECT DISTINCT * FROM ${dataset_name}.sra_processing"
# bq --project_id="${project_id}" --format=csv query --use_legacy_sql=false "${sql}"

#max_mem avg_mem swap stat exit_code exec_cwd exec_host
#bjobs -u all -d -o "jobid job_name user submit_time start_time finish_time run_time cpu_used slots min_req_proc max_req_proc nthreads delimiter='^'" > jobs.csv
num_of_snapshots=$(( batches / num_of_jobs + 1 ))
echo "Row count: ${row_count}. Total number of batches: ${batches}, Number of jobs: ${num_of_jobs}, Number of snapshots: ${num_of_snapshots}."
