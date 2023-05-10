#!/usr/bin/env bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
# sudo yum install jq -y
# sudo service docker start
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