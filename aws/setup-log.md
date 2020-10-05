# AWS setup log

This document is to track setup steps for AWS resources that will be used for the ScPCA project

## S3 buckets

### ScPCA buckets
All of the raw data for the ScPCA data is being uploaded to S3 in the `ccdl-scpca-data` bucket.
Most of this upload is happening through the the MASV portal, which can be found at <https://ccdl-scpca.portal.massive.app>

The files uploaded that way are organized by MASV and the submitter, and will require moving around for better organization.

Files for input to workflows (fastq, etc.) will be stored in the `s3://ccdl-scpca-data/raw` path.

### Nextflow buckets

Working files for Nextflow are stored in `s3://nextflow-ccdl-data/work/`.
The `work/` bucket subdirectory is set so that files expire after 90 days.
For more permanent storage, output files should be saved to a different prefix within this bucket, or better, placed in the `s3://nextflow-ccdl-results/` bucket.
This results bucket has versioning enabled to protect against accidental overwrites.

## AMI creation

The AMI for batch jobs was set up manually, following the instructions at
https://www.nextflow.io/docs/latest/awscloud.html#custom-ami

Briefly:
1. We used the base AMI _Amazon ECS-Optimized Amazon Linux AMI_. (note, not Amazon Linux 2).
This was launched with a 50GB EBS volume for data (compared to the default 22) in addition to the base 8GB boot volume.
All volumes are encrypted.

2. All packages were updated with `sudo yum update`
3. `/etc/sysconfig/docker-storage` edited to include `--storage-opt dm.basesize=50GB`
4. Install AWS CLI with:
```
sudo yum install -y bzip2 wget
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -f -p $HOME/miniconda
$HOME/miniconda/bin/conda install -c conda-forge -y awscli
rm Miniconda3-latest-Linux-x86_64.sh
```
5. Saved AMI from AWS web interface

**The current stored AMI is** `ccdl-nextflow-base-v1.1`: `ami-0a17541ba17115761`


## Terraform setup

The setup of the AWS batch queue and related resources for Nextflow usage are handled by the Terraform (`.tf`) files in this directory.


The general setup currently consists of two Batch queues, `nextflow-batch-default-queue`, which includes up to 100 vCPUs and `nextflow-batch-priority-queue`, which includes up to 20 vCPUs.
`nextflow-batch-default-queue` should be used in all but exceptional circumstances, as it uses spot instances for its processing, at a significant cost savings.
The spot price threshold is set high (currently 100%!) so it should never not run.

To use these files, you will need to install Terraform on your local machine, either from https://www.terraform.io/downloads.html or via Conda, Homebrew, or the like.
(You will also probably want the `aws` cli tools, and to have run `aws configure` to set your access key ID and  secret access key.)

Once you have installed terraform, you can run `terraform init` from this directory to get started.
To implement any changes, you will want a copy of the current state file, which is not included in this repository for security reasons.
With that file in the directory, you can test any changes with `terraform plan` and then implement them on AWS with `terraform apply`.

Changes to the compute environments do not always seem to go smoothly, as Terraform does not always properly shut down the job queues.
If you encounter an error with apply, you may need to run `terraform taint aws_batch_job_queue.nf_default_queue && terraform taint aws_batch_job_queue.nf_priority_queue` to force recreation of the job queues.