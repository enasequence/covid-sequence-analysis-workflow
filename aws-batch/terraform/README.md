# Terraform for Nextflow on AWS
This repository is for deploying AWS components and defining roles/policies needed to run the covid19-analysis pipeline on AWS Batch.

There are two approaches to setup the infrastructure for running Nextflow on AWS Batch. The first approach is to use `Tower Forge` which will automatically create all AWS components for you but only required IAM user to grant permission on AWS resource creation. 

Another approach is to create your AWS batch components by yourself. This repository is the Terraform code for creating AWS Batch components for the pipeline that use hybrid cloud (GCP&AWS).

## Modules
- AWS Batch
## Pre-script
- Run `aws configure sso` to create aws session.
- Run `./aws-batch/publish_head_node_ecr.sh` in [covid-sequence-analysis-workflow](https://github.com/enasequence/covid-sequence-analysis-workflow) to publish ECR repository for head node
- Run `packer build .` in /aws-batch directory to deploy custom AMI for worker nodes.
## Apply infrastructure
Run `./tf-apply.sh` 
## Post-script
### Pipeline repository
- Get job queue / job definition of head node and substitute in `start.aws.batch.sh` script.
- Get job queue of worker node substitute in `nextflow-lib/nextflow.config` script.
- Run `start.aws.batch.sh` script to submit job to head node and start the pipeline.


## Resource

[Packer build custom AMI](https://github.com/aws/amazon-ecs-ami)