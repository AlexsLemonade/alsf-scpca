# AWS setup log

This document is to track setup steps for AWS resources that will be used for the ScPCA project

## S3 bucket

All of the raw data for the ScPCA data is being uploaded to S3 in the `ccdl-scpca-data` bucket.
Most of this upload is happening through the the MASV portal, which can be found at <https://ccdl-scpca.portal.massive.app>

The files uploaded that way are organized by MASV and the submitter, and will require moving around for better organization.

Files for input to workflows (fastq, etc.) will be stored in the `s3://ccdl-scpca-data/raw` path.

## AMI creation

The AMI for batch jobs was set up manually, following the instructions at
https://www.nextflow.io/docs/latest/awscloud.html#custom-ami

Briefly:
1. We used the base AMI _Amazon ECS-Optimized Amazon Linux AMI_. (note, not Amazon Linux 2).
This was launched with an 80GB EBS volume for data (compared to teh default 22) in addition to the base 8GB boot volume.

2. All packages were updated with `sudo yum update`
3. `/etc/sysconfig/docker-storage` edited to include `--storage-opt dm.basesize=40GB`
4. Install AWS CLI with:
```
sudo yum install -y bzip2 wget
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -f -p $HOME/miniconda
$HOME/miniconda/bin/conda install -c conda-forge -y awscli
rm Miniconda3-latest-Linux-x86_64.sh
```
5. Saved AMI from AWS web interface

**The current stored AMI is** `ccdl-nextflow-base`: `ami-0a8857ac38c35157f`
