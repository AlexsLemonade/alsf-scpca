# AWS Batch setup
provider "aws" {
  profile = "default"
  region  = "us-east-1"
}

resource "aws_batch_job_queue" "nf_default_queue" {
  name     = "nextflow-batch-default-queue"
  state    = "ENABLED"
  priority = 1
  compute_environments = [
    aws_batch_compute_environment.nf_spot.arn,
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