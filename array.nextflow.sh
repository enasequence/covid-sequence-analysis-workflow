#!/usr/bin/env bash

# DIR where the current script resides
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
echo $DIR
DIR="/project/project_465000191/covid-sequence-analysis-workflow"
module use $HOME/.local/easybuild/modules/all
module load Nextflow/22.10.0
nextflow -version

batch_input=${1}
pipeline=${2:-'nanopore'}
profile=${3:-'codon'}
root_dir=${4:-'gs://prj-int-dev-covid19-nf-gls'}
batch_index=${5:-'0'}
snapshot_date=${6:-'2022-06-27'}
test_submission=${7:-'false'}
study_accession=${8:-'PRJEB45555'}
dataset_name=${9:-'sarscov2_metadata'}
project_id=${10:-'prj-int-dev-covid19-nf-gls'}

#################################
# Process the batch with Nextflow
#################################
batch_index=$SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID} / Batch index: ${batch_index}"
echo "** Processing samples with ${DIR}/${pipeline}/${pipeline}.nf. **"
echo "INPUT ${batch_input}_${SLURM_ARRAY_TASK_ID}.tsv"

pipeline_dir="${root_dir}/${snapshot_date}/${pipeline}_${batch_index}"
mkdir -p ${pipeline_dir}/{workDir,storeDir,publishDir};

nextflow -C "${DIR}/nextflow-lib/nextflow.config" run "${DIR}/${pipeline}/${pipeline}.nf" \
      -profile "${profile}" \
      --TEST_SUBMISSION "${test_submission}" --STUDY "${study_accession}" \
      --CONFIG_YAML "${DIR}/${pipeline}/config.yaml" \
      --SECRETS "${DIR}/data/projects_accounts.csv" \
      --SARS2_FA "${DIR}/data/NC_045512.2.fa" \
      --SARS2_FA_FAI "${DIR}/data/NC_045512.2.fa.fai" \
      --INDEX "${batch_input}_${SLURM_ARRAY_TASK_ID}.tsv" \
      --OUTDIR "${pipeline_dir}/publishDir" \
      --STOREDIR "${pipeline_dir}/storeDir" \
      -w "${pipeline_dir}/workDir" \
      -with-tower 
      # -stub-run 
########################################################################################
# Update submission receipt and submission metadata [as well as all the analyses archived]
########################################################################################
# "${DIR}/test.update.receipt.sh" "${batch_index}" "${snapshot_date}" "${pipeline}" "${profile}" "${root_dir}" "${dataset_name}" "${project_id}"
# "${DIR}/test.set.archived.sh" "${dataset_name}" "${project_id}"
# rm -rf "${pipeline_dir}/workDir/*" &
# rm -rf "${pipeline_dir}/storeDir/*" &
# rm -rf "${pipeline_dir}/publishDir/*" &
# wait
# rm -R "${pipeline_dir}"
