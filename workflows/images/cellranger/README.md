This folder contains a Dockerfile for the cellranger analysis.

The cellranger archive component must be downloaded separately to comply with licensing, and should be placed in this folder.
The current version of this file is `cellranger-4.0.0.tar.gz` and can be downloaded from [10X Genomics Website](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest) after agreeing to their license terms.


```
docker build . -t scpca-cellranger:4.0.0
```

This is then tagged and pushed to or pulled from AWS ECR with the following steps:

First, login to the AWS ECR docker repository with:

```
aws ecr get-login-password --region us-east-1 | docker login --username AWS --password-stdin 589864003899.dkr.ecr.us-east-1.amazonaws.com
```
Note that for the code above to work, you must have set up `aws` command line tools on your machine, and have run [`aws configure`](https://docs.aws.amazon.com/cli/latest/userguide/cli-configure-quickstart.html) to set up your credentials for AWS access.
This is a private repository; access is only available to Alexslemonade users at this time.


Following login (which only needs to be done once), tag and push the docker image with: 

```
docker tag scpca-cellranger:4.0.0 589864003899.dkr.ecr.us-east-1.amazonaws.com/scpca-cellranger:latest
docker tag scpca-cellranger:4.0.0 589864003899.dkr.ecr.us-east-1.amazonaws.com/scpca-cellranger:4.0.0
docker push 589864003899.dkr.ecr.us-east-1.amazonaws.com/scpca-cellranger:4.0.0
docker push 589864003899.dkr.ecr.us-east-1.amazonaws.com/scpca-cellranger:latest
```

To pull the current version of the container, you can use:

```
docker pull 589864003899.dkr.ecr.us-east-1.amazonaws.com/scpca-cellranger:latest
```

