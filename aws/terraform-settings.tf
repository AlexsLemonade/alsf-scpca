terraform {
  backend "s3" {
    bucket = "ccdl-nextflow-terraform-state"
    key    = "terraform.tfstate"
    encrypt = true
    region = "us-east-1"
  }
}
