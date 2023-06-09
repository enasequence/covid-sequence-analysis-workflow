err_file=log.err
run_accession=SRR24762121
executor=aws
error_status=ERROR
err_log_content=$(jq -Rs '.' ${err_file})
log_msg=$(cat <<-END
	{
        "severity": "${error_status}", 
        "message": ${err_log_content}, 
        "logging.googleapis.com/labels":{ 
            "run_accession" : "$run_accession", 
            "compute_resource" : "$executor"
        },
        "times":"$(date '+%Y-%m-%d %H:%M:%S')"
	}
	END
)
gcloud logging write covid_pipeline_logs "${log_msg}" --project="prj-int-dev-covid19-nf-gls" --payload-type=json --severity=ERROR