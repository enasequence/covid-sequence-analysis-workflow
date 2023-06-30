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

echo "$profile"
if [ "$profile" = "codonslurm"]; then
      echo "in"
      DIR="/hps/nobackup/tburdett/ena/users/analyser/covid-sequence-analysis-workflow"
fi
#################################
# Process the batch with Nextflow
#################################
echo ""
echo "** Processing samples with ${DIR}/${pipeline}/${pipeline}.nf. **"
mkdir -p "$DIR/secrets/"
gcloud_key_file="$DIR/secrets/${project_id}-sa-credential.json"


if [ "$profile" = "awsbatch" ]; then
      echo "** Retrieve secrets to $gcloud_key_file **"
      aws secretsmanager get-secret-value --secret-id $GOOGLE_APPLICATION_CREDENTIALS_SECRET_ID --query SecretString --output text > $gcloud_key_file
      gcloud auth activate-service-account --key-file=$gcloud_key_file
      gcloud config set project ${project_id}
      project_bucket="prj-int-dev-ait-eosc-aws-eval"
      aws s3 cp ${batch_input} ${DIR}/data/ # download sample index file from s3 to local dir
      aws s3 cp "s3://${project_bucket}/${dataset_name}/" "${DIR}/data/" --recursive --exclude "*/*"  # download projects_accounts and .fa files
      batch_input="${DIR}/data/$(basename -- "$batch_input")" #local path to sample index file
else
      touch ${gcloud_key_file}
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
      --GCLOUD_PROJECT_ID "${project_id}" \
      --GCLOUD_KEY_FILE "${gcloud_key_file}" \
      -w "${pipeline_dir}/workDir" \
      -with-tower

########################################################################################
# Update submission receipt and submission metadata [as well as all the analyses archived]
########################################################################################
"${DIR}/update.receipt.sh" "${batch_index}" "${snapshot_date}" "${pipeline}" "${profile}" "${root_dir}" "${dataset_name}" "${project_id}"
"${DIR}/set.archived.sh" "${dataset_name}" "${project_id}"

if ! [[ "$profile" =~ ^(gls|awsbatch)$ ]]; then
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

