This folder contains a Dockerfile for the cellranger analysis.


## Building the image

In order to build this image, the cellranger software and bcl2fastq source code components must be downloaded separately to comply with licensing, and should be placed in this folder (`images/cellranger`).
- The current version of cellranger is `cellranger-6.1.2.tar.gz` and can be downloaded from [10X Genomics Website](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/6.0) after agreeing to their license terms.
- The current version of bcl2fastq is 2.20.0 and can be downloaded from [Illumina](https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html) after logging into an Illumina account.
  - As we are using an Centos image, you will need to download the "Linux rpm" version of the software.
  - Note that the file Illumina provides to download is a `.zip` file that should contain the `bcl2fastq2-v2.20.0.422-Linux-x86_64.rpm` file specified in the `Dockerfile`. 

Following download of cellranger and bcl2fastq, you can build the image running the following command from this `images/cellranger` working directory:

```
docker build . -t scpca-cellranger:6.1.2
```

At this point, the image should be ready for use on the local machine.

## Using the AWS ECR

To use the AWS ECR docker repository, you will need to login  with:
```
aws ecr get-login-password --region us-east-1 | docker login --username AWS --password-stdin 589864003899.dkr.ecr.us-east-1.amazonaws.com
```

Note that for the code above to work, you must have set up `aws` command line tools on your machine, and have run [`aws configure`](https://docs.aws.amazon.com/cli/latest/userguide/cli-configure-quickstart.html) to set up your credentials for AWS access.
This is a private repository; access is only available to Alex's Lemonade (CCDL) users at this time.

### Pulling the current image

To pull the current version of the container (skipping the build steps), you can use the following command (once you are logged in):

```
docker pull 589864003899.dkr.ecr.us-east-1.amazonaws.com/scpca-cellranger:latest
```

### Pushing to AWS ECR

If the updated image needs to be pushed to AWS ECR, you can follow the outline steps below (updating version numbers as needed).
*Do not push an image unless you are sure it is working!*

The current image was pushed with the following commands.

```
docker tag scpca-cellranger:6.1.2 589864003899.dkr.ecr.us-east-1.amazonaws.com/scpca-cellranger:latest
docker tag scpca-cellranger:6.1.2 589864003899.dkr.ecr.us-east-1.amazonaws.com/scpca-cellranger:6.1.2
docker push 589864003899.dkr.ecr.us-east-1.amazonaws.com/scpca-cellranger:6.1.2
docker push 589864003899.dkr.ecr.us-east-1.amazonaws.com/scpca-cellranger:latest
```
