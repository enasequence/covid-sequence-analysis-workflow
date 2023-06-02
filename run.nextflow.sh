#!/usr/bin/env bash

# DIR where the current script resides
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

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
echo ""
echo "** Processing samples with ${DIR}/${pipeline}/${pipeline}.nf. **"
pipeline_dir="${root_dir}/${snapshot_date}/${pipeline}_${batch_index}"

if [ "$profile" = "lumislurm" ]; then 
      mkdir -p ${pipeline_dir}/{workDir,storeDir,publishDir}
      ${DIR}/hq/start.hq.node.sh "${pipeline_dir}/workDir"
      module use $HOME/.local/easybuild/modules/all
      module load Nextflow/22.10.0
fi

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


if [ "$profile" = "lumislurm" ]; then
      hq worker stop all
      hq server stop
fi
########################################################################################
# Update submission receipt and submission metadata [as well as all the analyses archived]
########################################################################################
# "${DIR}/update.receipt.sh" "${batch_index}" "${snapshot_date}" "${pipeline}" "${profile}" "${root_dir}" "${dataset_name}" "${project_id}"
# "${DIR}/set.archived.sh" "${dataset_name}" "${project_id}"

rm -R "${pipeline_dir}/workDir" &
rm -R "${pipeline_dir}/storeDir" &
rm -R "${pipeline_dir}/publishDir" &
wait
rm -R "${pipeline_dir}"
