This folder contains a Dockerfile for the cellranger analysis.


## Building the image

In order to build this image, the spaceranger software source code components must be downloaded separately to comply with licensing, and should be placed in this folder (`images/spaceranger`).
- The current version of cellranger is `spaceranger-1.3.1.tar.gz` and can be downloaded from [10X Genomics Website](https://support.10xgenomics.com/spatial-gene-expression/software/downloads/latest) after agreeing to their license terms.

Following download of spaceranger, you can build the image running the following command from this `images/spaceranger` working directory:

```
docker buildx build . -t scpca-spaceranger:1.3.1 --platform linux/amd64
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
docker pull 589864003899.dkr.ecr.us-east-1.amazonaws.com/scpca-spaceranger:latest
```

### Pushing to AWS ECR

If the updated image needs to be pushed to AWS ECR, you can follow the outline steps below (updating version numbers as needed).
*Do not push an image unless you are sure it is working!*

The current image was pushed with the following commands.

```
docker tag scpca-spaceranger:1.3.1 589864003899.dkr.ecr.us-east-1.amazonaws.com/scpca-spaceranger:latest
docker tag scpca-spaceranger:1.3.1 589864003899.dkr.ecr.us-east-1.amazonaws.com/scpca-spaceranger:1.3.1
docker push 589864003899.dkr.ecr.us-east-1.amazonaws.com/scpca-spaceranger:1.3.1
docker push 589864003899.dkr.ecr.us-east-1.amazonaws.com/scpca-spaceranger:latest
```
