variable "image_ids" {
  type = map
  default = {
    "nf-davidyuyuan-ena-sars-cov2-illumina-2-0" = "davidyuyuan/ena-sars-cov2-illumina:2.0"
    "nf-davidyuyuan-ena-sars-cov2-nanopore-2-0" = "davidyuyuan/ena-sars-cov2-nanopore:2.0"
    "nf-sands0-ena-analysis-submitter-2-3"  = "sands0/ena-analysis-submitter:2.3"
  }
}

variable "nf-token" {
  type = string
}