export $(grep -v '^#' .env | xargs)
terraform init -upgrade
terraform plan
terraform apply