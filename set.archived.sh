#!/usr/bin/env bash

dataset_name=${1:-'sarscov2_metadata'}
project_id=${2:-'prj-int-dev-covid19-nf-gls'}

# DIR where the current script resides
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

output_dir="${DIR}/results"; mkdir -p "${output_dir}"

###########################################################
# Get all analyses archived under parent_study = PRJEB45555: PRJEB55347, PRJEB55357, PRJEB57999, PRJEB59442, PRJEB61672
# PRJEB55347: PRJEB43947, PRJEB45554, PRJEB45619
# PRJEB55357: PRJEB55349, PRJEB55355, PRJEB55356
# PRJEB57999: PRJEB57992, PRJEB57993, PRJEB57995
# PRJEB59442: PRJEB59443, PRJEB59444, PRJEB59445
# PRJEB61672: PRJEB61667, PRJEB61668, PRJEB61669
###########################################################
echo ""
echo "** Updating ${dataset_name}.analysis_archived table. **"

# Parent studies <run_ref>
#curl -X POST -H "Content-Type: application/x-www-form-urlencoded" -d 'result=analysis&query=parent_study%3D%22PRJEB61672%22%20OR%20parent_study%3D%22PRJEB55347%22%20OR%20parent_study%3D%22PRJEB55357%22%20OR%20parent_study%3D%22PRJEB57999%22%20OR%20parent_study%3D%22PRJEB59442%22&fields=ALL&format=tsv&limit=10' \
#  "https://www.ebi.ac.uk/ena/portal/api/search" > "${output_dir}/analysis_archived.tsv"
#curl -X POST -H "Content-Type: application/x-www-form-urlencoded" -d 'result=analysis&query=parent_study%3D%22PRJEB61672%22%20OR%20parent_study%3D%22PRJEB55347%22%20OR%20parent_study%3D%22PRJEB55357%22%20OR%20parent_study%3D%22PRJEB57999%22%20OR%20parent_study%3D%22PRJEB59442%22&fields=sample_accession%2Crun_ref%2Canalysis_date%2Csubmitted_bytes%2Canalysis_type%2Cfirst_public&format=tsv&limit=0' \
#  "https://www.ebi.ac.uk/ena/portal/api/search" > "${output_dir}/analysis_archived.tsv"

# Child studies <run_accession>
#curl -X POST -H "Content-Type: application/x-www-form-urlencoded" -d 'result=analysis&query=study_accession%3D%22PRJEB43947%22%20OR%20study_accession%3D%22PRJEB45554%22%20OR%20study_accession%3D%22PRJEB45619%22%20OR%20study_accession%3D%22PRJEB55349%22%20OR%20study_accession%3D%22PRJEB55355%22%20OR%20study_accession%3D%22PRJEB55356%22%20OR%20study_accession%3D%22PRJEB57992%22%20OR%20study_accession%3D%22PRJEB57993%22%20OR%20study_accession%3D%22PRJEB57995%22%20OR%20study_accession%3D%22PRJEB59443%22%20OR%20study_accession%3D%22PRJEB59444%22%20OR%20study_accession%3D%22PRJEB59445%22%20OR%20study_accession%3D%22PRJEB61667%22%20OR%20study_accession%3D%22PRJEB61668%22%20OR%20study_accession%3D%22PRJEB61669%22&fields=ALL&format=tsv&limit=10' \
#  "https://www.ebi.ac.uk/ena/portal/api/search" > "${output_dir}/analysis_archived.tsv"
curl -X POST -H "Content-Type: application/x-www-form-urlencoded" -d 'result=analysis&query=study_accession%3D%22PRJEB43947%22%20OR%20study_accession%3D%22PRJEB45554%22%20OR%20study_accession%3D%22PRJEB45619%22%20OR%20study_accession%3D%22PRJEB55349%22%20OR%20study_accession%3D%22PRJEB55355%22%20OR%20study_accession%3D%22PRJEB55356%22%20OR%20study_accession%3D%22PRJEB57992%22%20OR%20study_accession%3D%22PRJEB57993%22%20OR%20study_accession%3D%22PRJEB57995%22%20OR%20study_accession%3D%22PRJEB59443%22%20OR%20study_accession%3D%22PRJEB59444%22%20OR%20study_accession%3D%22PRJEB59445%22%20OR%20study_accession%3D%22PRJEB61667%22%20OR%20study_accession%3D%22PRJEB61668%22%20OR%20study_accession%3D%22PRJEB61669%22&fields=sample_accession%2Crun_accession%2Canalysis_date%2Csubmitted_bytes%2Canalysis_type%2Cfirst_public&format=tsv&limit=0' \
  "https://www.ebi.ac.uk/ena/portal/api/search" > "${output_dir}/analysis_archived.tsv"
gsutil -m cp "${output_dir}/analysis_archived.tsv" "gs://${dataset_name}/analysis_archived.tsv" && \
  bq --project_id="${project_id}" load --source_format=CSV --replace=true --skip_leading_rows=1 --field_delimiter=tab \
  --autodetect "${dataset_name}.analysis_archived" "gs://${dataset_name}/analysis_archived.tsv" \
  "analysis_accession:STRING,sample_accession:STRING,run_accession:STRING,analysis_date:DATE,submitted_bytes:STRING,analysis_type:STRING,first_public:DATE"

#################################
# delete runs from sra_processing
#################################
echo ""
echo "** Updating ${dataset_name}.sra_processing table. **"

sql="DELETE FROM ${dataset_name}.sra_processing T1 WHERE T1.run_accession IN (SELECT T2.run_accession FROM ${dataset_name}.submission_metadata T2) OR T1.run_accession IN (SELECT T3.run_accession FROM ${dataset_name}.analysis_archived T3)"
bq --project_id="${project_id}" --format=csv query --use_legacy_sql=false "${sql}"
