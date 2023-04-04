###
# This script will publish the Docker image which contains `run.aws.nextflow.sh` / `nextflow.config`
# This reposiroty (Branch aws-batch) will then be cloned in run.aws.nextflow.sh
###
# sudo yum install jq -y
export $(grep -v '^#' .env | xargs)
echo $AWS_DEFAULT_REGION
# Create ECR for head node
# aws ecr get-login-password | docker login --username AWS --password-stdin ${AWS_ACCOUNT_ID}.dkr.ecr.${AWS_DEFAULT_REGION}.amazonaws.com
# aws ecr create-repository --tags Key=covid-pipeline,Value=true  --repository-name nextflow-head
# export REPO_URI="${AWS_ACCOUNT_ID}.dkr.ecr.${AWS_DEFAULT_REGION}.amazonaws.com/nextflow-head"
export REPO_URI="quay.io/mingyanisa/covid-analysis-nf-head"
echo $REPO_URI
export IMG_TAG="latest" # $(date +%F)
# Log in to your Quay.io account
docker login quay.io
# docker build -t ${REPO_URI}:${IMG_TAG} .
docker build -t ${REPO_URI}:${IMG_TAG} .
docker tag ${REPO_URI}:${IMG_TAG} quay.io/username/repository:${IMG_TAG}
# docker push ${REPO_URI}:${IMG_TAG}
echo "IMAGE NAME = ${REPO_URI}:${IMG_TAG}"