# AWS Batch setup
provider "aws" {
  region  = "us-east-1"
}

variable "default_tags" {
  description = "Default resource tags"
  type        = map(string)
  default     = {
    team = "science"
    project = "scpca"
    purpose = "nextflow-batch"
    config = "https://github.com/AlexsLemonade/alsf-scpca/tree/main/aws"
  }

}

resource "aws_batch_job_queue" "nf_default_queue" {
  name     = "nextflow-batch-default-queue"
  tags     = var.default_tags
  state    = "ENABLED"
  priority = 1
  compute_environments = [
    aws_batch_compute_environment.nf_spot.arn,
  ]
}

resource "aws_batch_job_queue" "nf_bigdisk_queue" {
  name     = "nextflow-batch-bigdisk-queue"
  tags     = var.default_tags
  state    = "ENABLED"
  priority = 1
  compute_environments = [
    aws_batch_compute_environment.nf_spot_bigdisk.arn,
  ]
}

resource "aws_batch_job_queue" "nf_autoscale_queue" {
  name = "nextflow-batch-autoscale-queue"
  tags = var.default_tags
  state = "ENABLED"
  priority = 1
  compute_environments = [
    aws_batch_compute_environment.nf_spot_auto_scaled_ebs.arn,
  ]
}


resource "aws_batch_job_queue" "nf_priority_queue" {
  name     = "nextflow-batch-priority-queue"
  state    = "ENABLED"
  priority = 1
  compute_environments = [
    aws_batch_compute_environment.nf_ondemand.arn,
  ]
}
