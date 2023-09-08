//export $(grep -v '^#' .env | xargs)
packer {
  required_plugins {
    amazon = {
      version = ">= 0.0.2"
      source  = "github.com/hashicorp/amazon"
    }
  }
}

source "amazon-ebs" "worker-al2" {
  ami_name      = "worker-aws-linux2-custom-ami"
  instance_type = "m5.large"
  region        = "eu-west-2"
  ami_description = "Amazon Linux AMI 2.0.20230109 x86_64 ECS HVM GP2"
  launch_block_device_mappings {
    volume_size           = var.block_device_size_gb
    delete_on_termination = true
    volume_type           = "gp2"
    device_name           = "/dev/xvda" 
  }
  source_ami_filter {
    filters = {
      name = "amzn2-ami-ecs-hvm-2.0.20230109-x86_64-ebs"
    }
    most_recent = true
    owners      = ["amazon"]
  }
  ssh_interface = "public_ip"
  ssh_username = "ec2-user"
  tags = {
    os_version          = "Amazon Linux 2"
    ami_type            = "al2"
    // ami_version         = "2.0.${var.ami_version}"
    ami_version         = "2.0.20230109"
  }
}

build {
  name    = "learn-packer"
  sources = [
    "source.amazon-ebs.worker-al2"
  ]

  provisioner "shell" {
    inline_shebang = "/bin/sh -ex"
    inline = [
        "cd $HOME",
        "sudo yum install -y bzip2 wget",
        "wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh",
        "bash Miniconda3-latest-Linux-x86_64.sh -b -f -p $HOME/miniconda",
        "$HOME/miniconda/bin/conda install -c conda-forge -y awscli",
        "rm Miniconda3-latest-Linux-x86_64.sh",
        "./miniconda/bin/aws --version"
    ]
  }
}
