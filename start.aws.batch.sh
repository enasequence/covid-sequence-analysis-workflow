#!/usr/bin/env bash

# DIR where the current script resides
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
export $(grep -v '^#' .env | xargs)
skip=${1:-'0'}
concurrency=${2:-'500'}   # Maximum concurrency determined by the bottleneck - the submission server and storage space
# concurrency=${2:-'5'}
pipeline=${3:-'illumina'}   # nanopore,illumina
root_dir=${4:-'s3://prj-int-dev-ait-eosc-aws-eval/nextflow'} #s3://prj-int-dev-covid19-nf-aws
batch_size=${5:-'15000'}
profile=${6:-'awsbatch'}
snapshot_date=${7:-'2023-03-29'}  #2022-09-26 2022-10-24 2022-11-21 2022-12-19
dataset_name=${8:-'sarscov2_metadata'}
project_id=${9:-'prj-int-dev-covid19-nf-gls'}
project_bucket='prj-int-dev-ait-eosc-aws-eval'
test_submission='false'
input_dir="${DIR}/data/${snapshot_date}"; mkdir -p "${input_dir}"


# Row count and batches
table_name="${pipeline}_to_be_processed"
sql="SELECT count(*) AS total FROM ${project_id}.${dataset_name}.${table_name}"
# row_count=$(bq --project_id="${project_id}" --format=csv query --use_legacy_sql=false "${sql}" | grep -v total)
row_count=75000
############################################
# as defined as queueSize in nextflow.config
############################################
queue_size=100
batches=$(( row_count / batch_size + 1 ))
num_of_jobs=$(( concurrency / queue_size ))

aws --region ${AWS_DEFAULT_REGION} s3 mb ${root_dir}
aws s3 sync "${DIR}/data/" "s3://${project_bucket}/${dataset_name}/" --exclude '${DIR}/data/${snapshot_date}/*'

for (( batch_index=skip; batch_index<skip+num_of_jobs&&batch_index<batches; batch_index++ )); do
  
	offset=$((batch_index * batch_size))
	echo ""
	echo "** Retrieving and reserving batch ${batch_index} with the size of ${batch_size} from the offset of ${offset}. **"

	# sql="SELECT * FROM ${project_id}.${dataset_name}.${table_name} LIMIT ${batch_size} OFFSET ${offset}"
	# bq --project_id="${project_id}" --format=csv query --use_legacy_sql=false --max_rows="${batch_size}" "${sql}" \
	#   | awk 'BEGIN{ FS=","; OFS="\t" }{$1=$1; print $0 }' > "${input_dir}/${table_name}_${batch_index}.tsv"

	# aws s3 cp "${input_dir}/${table_name}_${batch_index}.tsv" "s3://${project_bucket}/${dataset_name}/${snapshot_date}/" 

	# gsutil -m cp "${input_dir}/${table_name}_${batch_index}.tsv" "gs://${dataset_name}/${snapshot_date}/${table_name}_${batch_index}.tsv" && \
    # bq --project_id="${project_id}" load --source_format=CSV --replace=false --skip_leading_rows=1 --field_delimiter=tab \
    # --max_bad_records=0 "${dataset_name}.sra_processing" "gs://${dataset_name}/${snapshot_date}/${table_name}_${batch_index}.tsv"

	###TODO: Insert AWS cli command to submit batch job (head node) here###
	cmd_override=$(cat <<-END
	{
	  "command": ["/usr/local/bin/run.aws.nextflow.sh", \
	  "s3://${project_bucket}/${dataset_name}/${snapshot_date}/${table_name}_${batch_index}.tsv", \
	  "${pipeline}", "awsbatch", "${root_dir}",\
	  "${batch_index}", "${snapshot_date}",\
	  "${test_submission}"], \
	  "environment": [ \
	    {"name": "TOWER_ACCESS_TOKEN", "value": "${TOWER_ACCESS_TOKEN}"}
		]
	}
	END
	)
	
	aws batch submit-job --job-name "submit-job-${snapshot_date}-${pipeline}-${batch_index}" --job-definition "head_node_job" \
	--job-queue "head_queue" --container-overrides "${cmd_override}"
	# break
done
num_of_snapshots=$(( batches / num_of_jobs + 1 ))
echo "Row count: ${row_count}. Total number of batches: ${batches}, Number of jobs: ${num_of_jobs}, Number of snapshots: ${num_of_snapshots}."

