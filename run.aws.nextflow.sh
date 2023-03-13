#!/usr/bin/env bash

set -ex
s3_input_path=${1:-'s3://prj-int-dev-ait-eosc-aws-eval/sarscov2_metadata/illumina_to_be_processed_0.tsv'}
pipeline=${2:-'illumina'}
profile=${3:-'awsbatch'}
root_dir=${4:-'s3://prj-int-dev-ait-eosc-aws-eval/nextflow'}
batch_index=${5:-'0'}
snapshot_date=${6:-'2022-06-27'}
test_submission=${7:-'true'}
study_accession=${8:-'PRJEB45555'}
dataset_name=${9:-'sarscov2_metadata'}
project_id=${10:-'prj-int-dev-covid19-nf-gls'}

#################################
# Setup headnode
#################################
#!/bin/bash
PIPELINE_URL=${PIPELINE_URL:-https://github.com/enasequence/covid-sequence-analysis-workflow}
BRANCH="aws-batch"
# DIR where the current script resides
DIR="/scratch/covid-sequence-analysis-workflow"
git clone -b ${BRANCH} ${PIPELINE_URL} ${DIR}
cd ${DIR}
# echo ">> Remove container from pipeline config if present."
# sed -i -e '/process.container/d' ./nextflow-lib/nextflow.config
aws s3 ls s3://prj-int-dev-ait-eosc-aws-eval/sarscov2_metadata/
aws s3 cp --region ${AWS_DEFAULT_REGION} ${s3_input_path} ${DIR}/data/
echo "s3_input_path: ${s3_input_path}"
batch_input="${DIR}/data/$(basename -- "$s3_input_path")"
echo "filename: ${batch_input}"

# "$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
#################################
# Process the batch with Nextflow
#################################
echo ""
echo "** Processing samples with ${DIR}/${pipeline}/${pipeline}.nf. **"
echo $profile
pipeline_dir="${root_dir}/${snapshot_date}/${pipeline}_${batch_index}"
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

# ########################################################################################
# # Update submission receipt and submission metadata [as well as all the analyses archived]
# ########################################################################################
# # "${DIR}/update.receipt.sh" "${batch_index}" "${snapshot_date}" "${pipeline}" "${profile}" "${root_dir}" "${dataset_name}" "${project_id}"
# # "${DIR}/set.archived.sh" "${dataset_name}" "${project_id}"

# aws s3 rm --recursive "${pipeline_dir}/workDir" --quiet &
# aws s3 rm --recursive "${pipeline_dir}/storeDir" --quiet &
# aws s3 rm --recursive "${pipeline_dir}/publishDir" --quiet &
# wait
# aws s3 rm --recursive "${pipeline_dir}" 
