#!/usr/bin/env bash

##############################
# Run once per client machine
# Not required for Cloud Shell
##############################

PROJECT=${1:-'prj-int-dev-covid19-nf-gls'}

# Service account key for NF Tower
export SERVICE_ACCOUNT_NAME=nextflow-service-account
export SERVICE_ACCOUNT_ADDRESS=${SERVICE_ACCOUNT_NAME}@${PROJECT}.iam.gserviceaccount.com
export SERVICE_ACCOUNT_KEY=${SERVICE_ACCOUNT_NAME}-private-key.json

gcloud iam service-accounts keys create --project "${PROJECT}" \
  --iam-account="${SERVICE_ACCOUNT_ADDRESS}" \
  --key-file-type=json ${SERVICE_ACCOUNT_KEY}
export SERVICE_ACCOUNT_KEY_FILE=${DIR}/${SERVICE_ACCOUNT_KEY}
export GOOGLE_APPLICATION_CREDENTIALS=${DIR}/${SERVICE_ACCOUNT_KEY}

#https://www.nextflow.io/docs/latest/google.html?highlight=auth
gcloud auth application-default login
