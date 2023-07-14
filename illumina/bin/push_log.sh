#!/usr/bin/env bash
run_accession=${1}
err_file=${2}
profile=${3}
gcloud_project_id=${4}
CHAR_LIMIT=2000
if grep -qi "error" ${err_file};  then
    err_log_content=$(jq -Rs '.' ${err_file})
    log_msg=$(cat <<-END
        {
            "severity": "ERROR", 
            "message": ${err_log_content:0:CHAR_LIMIT}, 
            "logging.googleapis.com/labels":{ 
                "run_accession" : "$run_accession", 
                "compute_resource" : "$profile"
            }
        }
END
    )
    gcloud logging write covid_pipeline_logs "${log_msg}" --project="${gcloud_project_id}" --payload-type=json --severity=ERROR
fi