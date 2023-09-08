provider "aws" {
    profile = var.PROFILE
    region = var.REGION
}

module "aws_batch" {
  source = "./modules"
  nf-token = var.NF_TOKEN
}

terraform {
  required_providers {
      aws = {
        source  = "hashicorp/aws"
        version = "~> 4.60.0"
      }
  }
}