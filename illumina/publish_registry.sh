#!/usr/bin/env bash
###
# This script will publish the Docker image for headnode which contains `run.aws.nextflow.sh` / `nextflow.config`
# The repository (Branch aws-batch) will be cloned in the docker image.
###
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
### Option 2: Publish docker image to Quay.io
export REPO_URI="quay.io/yanisasunt/ena-sars-cov2-illumina"
export IMG_TAG="1.0" # $(date +%F)
echo $REPO_URI $DIR
# Log in to your Quay.io account
echo "$QUAY_PASSWORD" | docker login quay.io --username $QUAY_USERNAME --password-stdin
docker rmi ${REPO_URI}:${IMG_TAG}
docker build -t ${REPO_URI}:${IMG_TAG} ${DIR}
echo "IMAGE NAME = ${REPO_URI}:${IMG_TAG}"
docker tag ${REPO_URI}:${IMG_TAG} ${REPO_URI}:${IMG_TAG}
docker push ${REPO_URI}:${IMG_TAG}