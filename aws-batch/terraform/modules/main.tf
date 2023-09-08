# Security Setup
# Retrieves the default vpc for this region
data "aws_vpc" "default" {
  default = true
}

# Security Group for batch processing - default
data "aws_security_group" "default" {
  vpc_id = data.aws_vpc.default.id
  filter {
    name   = "group-name"
    values = ["default"]
  }
}
# Retrieves the subnet ids in the default vpc
data "aws_subnets" "all_default_subnets" {
  filter {
    name   = "vpc-id"
    values = [data.aws_vpc.default.id]
  }
}
data "aws_kms_key" "secretsmanager" {
  key_id = "alias/aws/secretsmanager"
}
resource "aws_secretsmanager_secret" "gcloud_credential_secret" {
    name = "gcloud_credential_secret"
    tags = {
        created-by = "terraform"
    }
}

resource "aws_secretsmanager_secret_version" "gcloud_credential_secret" {
    secret_id = aws_secretsmanager_secret.gcloud_credential_secret.id
    secret_string = file("../secrets/prj-int-dev-covid19-nf-gls-worker-node.json")
}

resource "aws_iam_role" "batch_ec2_role" {    
    name = "batch_ec2_role"
    assume_role_policy = jsonencode({
        Version= "2012-10-17",
        Statement = [
            {
                Action = "sts:AssumeRole",
                Principal = {
                    "Service":"ec2.amazonaws.com"
                },
                Effect = "Allow",
            }
        ]
    })
    tags = {
        created-by = "terraform"
    }
}

resource "aws_iam_policy" "read_gcloud_secret_policy" {
    name        = "GetGCPSecret"
    description = "Policy to get gcloud secrets"
    policy = jsonencode({
        Version = "2012-10-17"
        Statement = [
            {
                Effect = "Allow",
                Action = [
                    "secretsmanager:GetResourcePolicy",
                    "secretsmanager:GetSecretValue",
                    "secretsmanager:DescribeSecret",
                    "secretsmanager:ListSecretVersionIds"
                ],
                Resource = aws_secretsmanager_secret_version.gcloud_credential_secret.arn
            },
            {
                Effect= "Allow",
                Action= [
                    "secretsmanager:GetRandomPassword",
                    "secretsmanager:ListSecrets"
                ],
                Resource= "*"
            },
            {
                Effect= "Allow",
                Action= [
                    "kms:Decrypt"
                ],
                Resource= data.aws_kms_key.secretsmanager.arn
            }
        ]
    })
    depends_on = [
        aws_secretsmanager_secret_version.gcloud_credential_secret
    ]
}

# Service role policy which will be attached to CE
resource "aws_iam_role" "batch_role" {
    name               = "batch_role"
    assume_role_policy = jsonencode(
    {
        Version = "2012-10-17",
        Statement = [
        {
            Action = "sts:AssumeRole",
            Effect = "Allow",
            Principal = {
                "Service": "batch.amazonaws.com"
            }
        }
        ]
    }
    )
    tags = {
        created-by = "terraform"
    }
}

# Attach the Batch policy to the Batch role
resource "aws_iam_role_policy_attachment" "policy_attachment" {
    role       = aws_iam_role.batch_role.name
    for_each = toset([
        "arn:aws:iam::aws:policy/service-role/AWSBatchServiceRole",
        "arn:aws:iam::aws:policy/AmazonECS_FullAccess",
        aws_iam_policy.read_gcloud_secret_policy.arn,
    ])
    policy_arn = each.value
    depends_on = [
        aws_iam_policy.read_gcloud_secret_policy
    ]
}

# Assign the EC2 role to the EC2 profile
resource "aws_iam_instance_profile" "batch_ec2_profile" {
    name = "batch_ec2_profile"
    role = aws_iam_role.batch_ec2_role.name
}

# Attach the EC2 container service policy to the EC2 role 
resource "aws_iam_role_policy_attachment" "batch_ec2_policy_attachment" {
    role       = aws_iam_role.batch_ec2_role.name
    for_each = toset([
        "arn:aws:iam::aws:policy/service-role/AmazonEC2ContainerServiceforEC2Role",
        "arn:aws:iam::aws:policy/AmazonS3FullAccess",
        "arn:aws:iam::aws:policy/AWSBatchFullAccess",
        aws_iam_policy.read_gcloud_secret_policy.arn,
        # "arn:aws:iam::aws:policy/AmazonEC2ContainerRegistryFullAccess",
    ])
    policy_arn = each.value
    # policy_arn = "arn:aws:iam::aws:policy/service-role/AmazonEC2ContainerServiceforEC2Role"
}

# IAM Role for jobs
resource "aws_iam_role" "batch_job_role" {
    name               = "batch_job_role"
    assume_role_policy = jsonencode(
    {
        Version = "2012-10-17",
        Statement = [
        {
            Action = "sts:AssumeRole",
            Effect = "Allow",
            Principal = {
                "Service": "ecs-tasks.amazonaws.com"
            }
        }
        ]
    })
    tags = {
        created-by = "terraform"
    }
}

resource "aws_iam_role_policy_attachment" "batch_job_policy_attachment" {
    role       = aws_iam_role.batch_job_role.name
    for_each = toset([
        "arn:aws:iam::aws:policy/AmazonS3FullAccess",
        "arn:aws:iam::aws:policy/AWSBatchFullAccess",
        "arn:aws:iam::aws:policy/AmazonECS_FullAccess",
        aws_iam_policy.read_gcloud_secret_policy.arn,
    ])
    policy_arn = each.value
}

# Create head node CE
resource "aws_batch_compute_environment" "head_node_ce" {
    compute_environment_name = "head-node-ce"
    compute_resources {
        instance_role = aws_iam_instance_profile.batch_ec2_profile.arn
        instance_type = [
            "optimal"
        ]
        max_vcpus = 15
        min_vcpus = 0
        security_group_ids = [
            data.aws_security_group.default.id
        ]
        subnets = data.aws_subnets.all_default_subnets.ids
        type    = "EC2"
    }
    service_role = aws_iam_role.batch_role.arn
    type = "MANAGED"
    tags = {
        created-by = "terraform"
    }
}

# Create head node job queue
resource "aws_batch_job_queue" "head_queue" {
    name     = "head_queue"
    state    = "ENABLED"
    priority = 1
    compute_environments = [
        aws_batch_compute_environment.head_node_ce.arn
    ]
    depends_on = [aws_batch_compute_environment.head_node_ce]
    tags = {
        created-by = "terraform"
    }
}

resource "aws_batch_job_definition" "head_node_def" {
    name = "head_node_job"
    type = "container"
    parameters = {}
    container_properties = jsonencode(
    {
        image = "quay.io/mingyanisa/covid-analysis-nf-head",
        jobRoleArn = "${aws_iam_role.batch_job_role.arn}",
        vcpus = 1,
        memory = 8192,
        environment = [
            { 
                name = "GOOGLE_APPLICATION_CREDENTIALS_SECRET_ARN",
                value = aws_secretsmanager_secret_version.gcloud_credential_secret.arn
            }
        ],
        command = ["/usr/local/bin/run.aws.nextflow.sh"],
        user = "root",
    }
    )
    tags = {
        created-by = "terraform"
    }
}
# AMI for worker nodes
data "aws_ami" "worker_custom_ami" {
  filter {
    name   = "name"
    values = ["worker-aws-linux2-custom-ami"]
  }
}

resource "aws_batch_job_definition" "worker_node_def" {
    for_each    = var.image_ids
    name = each.key
    type = "container"
    parameters = {
        nf-token = var.nf-token
    }
    container_properties = jsonencode(
    {
        image = each.value,
        vcpus = 1, 
        memory = 1024, 
        environment = [],
        volumes = [
            {
                "host": {
                    "sourcePath": "/home/ec2-user/miniconda"
                },
                "name": "aws-cli"
            }
        ],
        mountPoints = [
            {
                "containerPath": "/home/ec2-user/miniconda",
                "readOnly": true,
                "sourceVolume": "aws-cli"
            }
        ],
        command = ["true"]
    })
    # timeout {
    #     attempt_duration_seconds=3600
    # }
    tags = {
        created-by = "terraform"
    }
}

# Create worker node CE
resource "aws_batch_compute_environment" "worker_node_ce" {
    compute_environment_name = "worker-node-ce"
    compute_resources {
        instance_role = aws_iam_instance_profile.batch_ec2_profile.arn
        instance_type = [
            "optimal"
        ]
        max_vcpus = 5120 //1024*5
        min_vcpus = 0
        security_group_ids = [
            data.aws_security_group.default.id
        ]
        subnets = data.aws_subnets.all_default_subnets.ids
        type    = "SPOT" # EC2
        allocation_strategy = "BEST_FIT_PROGRESSIVE"
        bid_percentage = 70
        ec2_configuration {
            image_id_override = data.aws_ami.worker_custom_ami.image_id 
            image_type = "ECS_AL2"
        }
    }
    service_role = aws_iam_role.batch_role.arn
    type = "MANAGED"
    tags = {
        created-by = "terraform"
    }
}

# Create worker node job queue
resource "aws_batch_job_queue" "worker_queue" {
    name     = "worker_queue"
    state    = "ENABLED"
    priority = 1
    compute_environments = [
        aws_batch_compute_environment.worker_node_ce.arn
    ]
    depends_on = [aws_batch_compute_environment.worker_node_ce]
    tags = {
        created-by = "terraform"
    }
}