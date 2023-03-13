export $(grep -v '^#' .env | xargs)
aws ecr create-repository --tags Key=covid-pipeline,Value=true  --repository-name covid-sequence-analysis-workflow-head
export REPO_URI=$(aws ecr describe-repositories --repository-names=nextflow-head |jq -r '.repositories[0].repositoryUri')
echo $REPO_URI
chmod +x ${HOME}/projects/covid-sequence-analysis-workflow/run.aws.nextflow.sh
export IMG_TAG=$(date +%F).1
docker build -t ${REPO_URI}:${IMG_TAG} .
docker push ${REPO_URI}:${IMG_TAG}
echo "IMAGE NAME = ${REPO_URI}:${IMG_TAG}"
# echo "BUCKET_NAME_RESULTS = ${BUCKET_NAME_RESULTS}"