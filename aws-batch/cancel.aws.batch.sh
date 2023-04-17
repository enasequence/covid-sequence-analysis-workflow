#!/usr/bin/env bash
source .env
aws ecr get-login-password | docker login --username AWS --password-stdin ${AWS_ACCOUNT_ID}.dkr.ecr.${AWS_DEFAULT_REGION}.amazonaws.com
for i in $(aws batch list-jobs --job-queue head_queue --job-status running --output text --query "jobSummaryList[*].[jobId]")
do
  echo "Terminate Job: $i"
  aws batch terminate-job --job-id $i --reason "Terminating job."
  echo "Job $i terminated"
done
for i in $(aws batch list-jobs --job-queue worker_queue --job-status runnable --output text --query "jobSummaryList[*].[jobId]")
do
  echo "Cancel Job: $i"
  aws batch cancel-job --job-id $i --reason "Cancelling job."
  echo "Job $i canceled"
done
for i in $(aws batch list-jobs --job-queue worker_queue --job-status running --output text --query "jobSummaryList[*].[jobId]")
do
  echo "Terminate Job: $i"
  aws batch terminate-job --job-id $i --reason "Terminating job."
  echo "Job $i terminated"
done