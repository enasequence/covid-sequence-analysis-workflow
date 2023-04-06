#!/usr/bin/env bash

# DIR where the current script resides
if [ "$profile" = "awsbatch" ]; then
      DIR="/scratch/covid-sequence-analysis-workflow"
else
      DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
fi

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

#################################
# Process the batch with Nextflow
#################################
echo ""
echo "** Processing samples with ${DIR}/${pipeline}/${pipeline}.nf. **"

if [ "$profile" = "awsbatch" ]; then
      echo "** Retrieve secrets: ${SERVICE_ACCOUNT_KEY_FILE} **"
      aws secretsmanager get-secret-value --secret-id $GOOGLE_APPLICATION_CREDENTIALS_SECRET_ARN:project_id:: --query SecretString --output text
      aws secretsmanager get-secret-value --secret-id $GOOGLE_APPLICATION_CREDENTIALS_SECRET_ARN --query SecretString --output text > $SERVICE_ACCOUNT_KEY_FILE 
      if [[ ! -f $SERVICE_ACCOUNT_KEY_FILE  ]] ; then
            echo "File ${SERVICE_ACCOUNT_KEY_FILE} is not there, aborting."
            exit
      fi
      gcloud auth activate-service-account --key-file=$SERVICE_ACCOUNT_KEY_FILE
      gcloud config set project ${project_id}
      project_bucket="prj-int-dev-ait-eosc-aws-eval"
      aws s3 cp ${batch_input} ${DIR}/data/ # download sample index file from s3 to local dir
      aws s3 cp "s3://${project_bucket}/${dataset_name}/" "${DIR}/data/" --recursive --exclude "*/*"  # download projects_accounts and .fa files
      batch_input="${DIR}/data/$(basename -- "$batch_input")" #local path to sample index file
fi

pipeline_dir="${root_dir}/${snapshot_date}/${pipeline}_${batch_index}"
echo "** pipeline_dir: ${pipeline_dir} **"
nextflow -C "${DIR}/nextflow-lib/nextflow.config" run "${DIR}/${pipeline}/${pipeline}.nf" -profile "${profile}" \
      --TEST_SUBMISSION "${test_submission}" --STUDY "${study_accession}" \
      --CONFIG_YAML "${DIR}/${pipeline}/config.yaml" \
      --SECRETS "${DIR}/data/projects_accounts.csv" \
      --SARS2_FA "${DIR}/data/NC_045512.2.fa" \
      --SARS2_FA_FAI "${DIR}/data/NC_045512.2.fa.fai"\
      --INDEX "${batch_input}" \
      --OUTDIR "${pipeline_dir}/publishDir" \
      --STOREDIR "${pipeline_dir}/storeDir" \
      -w "${pipeline_dir}/workDir" \
      -with-tower

########################################################################################
# Update submission receipt and submission metadata [as well as all the analyses archived]
########################################################################################
"${DIR}/update.receipt.sh" "${batch_index}" "${snapshot_date}" "${pipeline}" "${profile}" "${root_dir}" "${dataset_name}" "${project_id}"
"${DIR}/set.archived.sh" "${dataset_name}" "${project_id}"

if [ "$profile" != "gls" ] && [ "$profile" != "awsbatch" ]; then
      rm -R "${pipeline_dir}/workDir" &
      rm -R "${pipeline_dir}/storeDir" &
      rm -R "${pipeline_dir}/publishDir" &
      wait
      rm -R "${pipeline_dir}"
fi

if [ "$profile" = "awsbatch" ]; then
      aws s3 rm --recursive "${pipeline_dir}/workDir" --quiet &
      aws s3 rm --recursive "${pipeline_dir}/storeDir" --quiet &
      aws s3 rm --recursive "${pipeline_dir}/publishDir" --quiet &
      wait
      aws s3 rm --recursive "${pipeline_dir}"
fi

