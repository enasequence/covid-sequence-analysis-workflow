#!/usr/bin/env bash
run_accession=${1}
err_file=${2}
profile=${3}
gcloud_project_id=${4}
if grep -qi "error" ${err_file};  then
    err_log_content=$(jq -Rs '.' ${err_file})
    log_msg=$(cat <<-END
        {
            "severity": "ERROR", 
            "message": ${err_log_content}, 
            "logging.googleapis.com/labels":{ 
                "run_accession" : "$run_accession", 
                "compute_resource" : "$profile"
            }
        }
END
    )
    gcloud logging write covid_pipeline_logs "${log_msg}" --project="${gcloud_project_id}" --payload-type=json --severity=ERROR
    # ## Send the log data to Datadog using the API
    # err_log_content=$(jq -Rs '.' ${err_file})
    # log_msg=$(jq -n --arg error_status "ERROR" \
    #         --arg ts "$(date '+%Y-%m-%d %H:%M:%S')" \
    #         --arg error_msg "${err_log_content}" \
    #         '{severity:$error_status, message:$error_msg, times:$ts}' )
fi