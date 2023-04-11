###
# This script will publish the Docker image for headnode which contains `run.aws.nextflow.sh` / `nextflow.config`
# The repository (Branch aws-batch) will be cloned in the docker image.
###
# sudo yum install jq -y
export $(grep -v '^#' .env | xargs)
echo $AWS_DEFAULT_REGION
### Option 1: Create ECR for head node
# aws ecr get-login-password | docker login --username AWS --password-stdin ${AWS_ACCOUNT_ID}.dkr.ecr.${AWS_DEFAULT_REGION}.amazonaws.com
# aws ecr create-repository --tags Key=covid-pipeline,Value=true  --repository-name nextflow-head
# export REPO_URI="${AWS_ACCOUNT_ID}.dkr.ecr.${AWS_DEFAULT_REGION}.amazonaws.com/nextflow-head"

### Option 2: Publish docker image to Quay.io
export REPO_URI="quay.io/mingyanisa/covid-analysis-nf-head"
echo $REPO_URI
export IMG_TAG="latest" # $(date +%F)
# Log in to your Quay.io account
echo "$QUAY_PASSWORD" | docker login quay.io --username $QUAY_USERNAME --password-stdin
docker rmi ${REPO_URI}:${IMG_TAG}
docker build -t ${REPO_URI}:${IMG_TAG} .
echo "IMAGE NAME = ${REPO_URI}:${IMG_TAG}"
docker tag ${REPO_URI}:${IMG_TAG} ${REPO_URI}:${IMG_TAG}
docker push ${REPO_URI}:${IMG_TAG}