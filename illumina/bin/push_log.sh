#!/usr/bin/env bash
run_accession=${1}
err_file=${2}
# if grep -qi "error" ${err_file};  then
## Send the log data to Datadog using the API
err_log_content=$(jq -Rs '.' ${err_file})
echo $(cat << EOF
[
{
    "ddsource": "covid-pipeline",
    "ddtags": "env:aws-batch,version:5.1",
    "hostname": "aws",
    "message": ${err_log_content},
    "status": "error",
    "run_accession": "${run_accession}"
}
]
EOF
  ) | gzip | curl -X POST "https://http-intake.logs.datadoghq.eu/api/v2/logs" \
  -H "Accept: application/json" \
  -H "Content-Type: application/json" \
  -H "Content-Encoding: gzip" \
  -H "DD-API-KEY: ${DD_API_KEY}" \
  --data-binary @-
# fi