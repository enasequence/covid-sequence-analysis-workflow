#!/usr/bin/env bash

# DIR where the current script resides
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
batch_input=${1}
pipeline=${2:-'nanopore'}
profile=${3:-'codon'}
root_dir=${4:-'gs://prj-int-dev-covid19-nf-gls'}
batch_index=${5:-'0'}
snapshot_date=${6:-'2023-03-13'}
test_submission=${7:-'true'}
study_accession=${8:-'PRJEB45555'}
dataset_name=${9:-'sarscov2_metadata'}
project_id=${10:-'prj-int-dev-covid19-nf-gls'}
source .env
#################################
# Process the batch with Nextflow
#################################
echo ""
echo "** Processing samples with ${DIR}/${pipeline}/${pipeline}.nf. **"

# if [ "$profile" = "awsbatch" ]; then
#       echo "** Retrieve secrets to ${DIR}/${project_id}-sa-credential.json **"
#       aws secretsmanager get-secret-value --secret-id $GOOGLE_APPLICATION_CREDENTIALS_SECRET_ID --query SecretString --output text > "$DIR/${project_id}-sa-credential.json"
#       gcloud auth activate-service-account --key-file="$DIR/${project_id}-sa-credential.json"
#       gcloud config set project ${project_id}
#       project_bucket="prj-int-dev-ait-eosc-aws-eval"
#       aws s3 cp ${batch_input} ${DIR}/data/ # download sample index file from s3 to local dir
#       aws s3 cp "s3://${project_bucket}/${dataset_name}/" "${DIR}/data/" --recursive --exclude "*/*"  # download projects_accounts and .fa files
#       batch_input="${DIR}/data/$(basename -- "$batch_input")" #local path to sample index file
# fi

# pipeline_dir="${root_dir}/${snapshot_date}/${pipeline}_${batch_index}"
# echo "** pipeline_dir: ${pipeline_dir} **"
# nextflow -C "${DIR}/nextflow-lib/nextflow.config" run "${DIR}/${pipeline}/${pipeline}.nf" -profile "${profile}" \
#       --TEST_SUBMISSION "${test_submission}" --STUDY "${study_accession}" \
#       --CONFIG_YAML "${DIR}/${pipeline}/config.yaml" \
#       --SECRETS "${DIR}/data/projects_accounts.csv" \
#       --SARS2_FA "${DIR}/data/NC_045512.2.fa" \
#       --SARS2_FA_FAI "${DIR}/data/NC_045512.2.fa.fai"\
#       --INDEX "${batch_input}" \
#       --OUTDIR "${pipeline_dir}/publishDir" \
#       --STOREDIR "${pipeline_dir}/storeDir" \
#       -w "${pipeline_dir}/workDir" \
#       -with-tower

########################################################################################
# Update submission receipt and submission metadata [as well as all the analyses archived]
########################################################################################
# "${DIR}/update.receipt.sh" "${batch_index}" "${snapshot_date}" "${pipeline}" "${profile}" "${root_dir}" "${dataset_name}" "${project_id}"
# "${DIR}/set.archived.sh" "${dataset_name}" "${project_id}"
# DD_SITE="datadoghq.eu" DD_API_KEY="$DD_API_KEY" python3 "push_log.py"

if ! [[ "$profile" =~ ^(gls|awsbatch)$ ]]; then
      rm -R "${pipeline_dir}/workDir" &
      rm -R "${pipeline_dir}/storeDir" &
      rm -R "${pipeline_dir}/publishDir" &
      wait
      rm -R "${pipeline_dir}"
fi


# Set your S3 bucket and log file path
S3_BUCKET=prj-int-dev-ait-eosc-aws-eval
S3_LOG_FILE=nextflow/2023-04-24/illumina_1/workDir/48/ea615ef8ef0316cebf6cbc147db317/.command.err
# nextflow/2023-04-24/illumina_1/workDir/0a/a9132a0f3fdac34c05c65b60747006/.command.err

# Download the log file from S3 to a temporary file
TMP_FILE=$(mktemp)
aws s3 cp s3://${S3_BUCKET}/${S3_LOG_FILE} ${TMP_FILE}

# Embed the log file contents in the Datadog API request
LOG_DATA=$(cat ${TMP_FILE})
# | jq -R -s 'split("\n")[:-1] | map(split(" ")) | map({"timestamp": .[0], "message": .[1]})')
FTP_URL=$(cat ${TMP_FILE} | grep -oP -m 1 'ftp://ftp\.sra\.ebi\.ac\.uk/.*')
echo ${FTP_URL}
SAMPLE_ID=$(echo ${FTP_URL} | grep -oP '.*\\([^\\]+)\\')

# Print the sample_id
echo ${SAMPLE_ID}

# # Send the log data to Datadog using the API
# echo $(cat << EOF
# [
#   {
#     "ddsource": "nginx",
#     "ddtags": "env:staging,version:5.1",
#     "hostname": "i-012345678",
#     "message": "${LOG_DATA}",
#     "service": "payment"
#   }
# ]
# EOF
# ) | gzip | curl -X POST "https://http-intake.logs.datadoghq.eu/api/v2/logs" \
# -H "Accept: application/json" \
# -H "Content-Type: application/json" \
# -H "Content-Encoding: gzip" \
# -H "DD-API-KEY: ${DD_API_KEY}" \
# --data-binary @-

# if ! [[ "$profile" =~ ^(gls|awsbatch)$ ]]; then
#       rm -R "${pipeline_dir}/workDir" &
#       rm -R "${pipeline_dir}/storeDir" &
#       rm -R "${pipeline_dir}/publishDir" &
#       wait
#       rm -R "${pipeline_dir}"
# fi

# if [ "$profile" = "awsbatch" ]; then
#       aws s3 rm --recursive "${pipeline_dir}/workDir" --quiet &
#       aws s3 rm --recursive "${pipeline_dir}/storeDir" --quiet &
#       aws s3 rm --recursive "${pipeline_dir}/publishDir" --quiet &
#       wait
#       aws s3 rm --recursive "${pipeline_dir}"
# fi

