###
# This script will publish the Docker image which contains `run.aws.nextflow.sh` / `nextflow.config`
# This reposiroty (Branch aws-batch) will then be cloned in run.aws.nextflow.sh
###
# sudo yum install jq -y
export $(grep -v '^#' .env | xargs)
echo $AWS_DEFAULT_REGION
# Create ECR for head node
aws ecr get-login-password | docker login --username AWS --password-stdin ${AWS_ACCOUNT_ID}.dkr.ecr.${AWS_DEFAULT_REGION}.amazonaws.com
aws ecr create-repository --tags Key=covid-pipeline,Value=true  --repository-name nextflow-head
export REPO_URI=$(aws ecr describe-repositories --repository-names=nextflow-head |jq -r '.repositories[0].repositoryUri')
echo $REPO_URI
# chmod +x ${HOME}/run.aws.nextflow.sh
export IMG_TAG="latest"
# $(date +%F)
docker build -t ${REPO_URI}:${IMG_TAG} .
docker push ${REPO_URI}:${IMG_TAG}
echo "IMAGE NAME = ${REPO_URI}:${IMG_TAG}"

# Create Job CE / Job definition / Job Queue -> TF
