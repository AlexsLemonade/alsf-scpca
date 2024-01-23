terraform {
  backend "s3" {
    bucket = "ccdl-nextflow-terraform-state"
    key    = "terraform.tfstate"
    region = "us-east-1"
  }
}
