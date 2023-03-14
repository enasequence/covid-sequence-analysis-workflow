# sudo yum install jq -y
export $(grep -v '^#' .env | xargs)
echo $AWS_DEFAULT_REGION
# Create ECR for head node
aws ecr get-login-password | docker login --username AWS --password-stdin ${AWS_ACCOUNT_ID}.dkr.ecr.${AWS_DEFAULT_REGION}.amazonaws.com
aws ecr create-repository --tags Key=covid-pipeline,Value=true  --repository-name nextflow-head
export REPO_URI=$(aws ecr describe-repositories --repository-names=nextflow-head |jq -r '.repositories[0].repositoryUri')
echo $REPO_URI
chmod +x ${HOME}/projects/covid-sequence-analysis-workflow/run.aws.nextflow.sh
export IMG_TAG="latest"
# $(date +%F)
docker build -t ${REPO_URI}:${IMG_TAG} .
# docker push ${REPO_URI}:${IMG_TAG}
echo "IMAGE NAME = ${REPO_URI}:${IMG_TAG}"
# TF: Create Job CE / Job Queue
# Create new CE for head node

#Create new job def for head node
head_node_def_name="nextflow-parallel-head-node"
latest_revision=$(aws batch describe-job-definitions  --job-definition-name "${head_node_def_name}"  \
   --query='jobDefinitions[?status==`ACTIVE`].revision' --output=json | jq '.[0]')
aws batch deregister-job-definition --job-definition "${head_node_def_name}:${latest_revision}"
container_prop=$(cat <<-END
	{
        "image":"${REPO_URI}:${IMG_TAG}",
        "vcpus": 1,
        "memory": 1032,
        "command": ["/usr/local/bin/run.aws.nextflow.sh"],
        "jobRoleArn": "arn:aws:iam::${AWS_ACCOUNT_ID}:role/ecsTaskExecutionRole",
        "user": "root"
	}
	END
)
echo "${container_prop}"
# aws batch register-job-definition --job-definition-name "${job_def_name}" --type "container" \
#     --container-properties "${container_prop}"