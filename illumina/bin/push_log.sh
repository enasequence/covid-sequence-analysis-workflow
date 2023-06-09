#!/usr/bin/env bash
run_accession=${1}
err_file=${2}
if grep -qi "error" ${err_file};  then
    ## Send the log data to Datadog using the API
    err_log_content=$(jq -Rs '.' ${err_file})
    log_msg=$(jq -n --arg error_status "ERROR" \
            --arg ts "$(date '+%Y-%m-%d %H:%M:%S')" \
            --arg error_msg "${err_log_content}" \
            '{severity:$error_status, message:$error_msg, times:$ts}' )
fi